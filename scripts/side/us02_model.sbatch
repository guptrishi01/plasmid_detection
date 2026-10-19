#!/bin/bash
#SBATCH --job-name=pd_us_model
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# us02_model.sbatch -- side analysis (not in the spec): summarise the 92 Nonpareil curves
# from us01, then model the odds of under-sequencing (logistic vs beta regression, compared
# by leave-one-out CV; scripts/R/undersequencing.R). Also tabulates, for the pilot runs,
# raw-read vs host-depleted Nonpareil coverage (what using raw reads in us01 costs) next
# to the stage 4 assembled fraction.
# Output: $paths.undersequencing/ (nonpareil_all.tsv, model_comparison_loocv.tsv,
#         odds_ratios.tsv, per_run_predictions.tsv, *_summary.txt, *.png, pilot_calibration.tsv)
# Submit:  sbatch scripts/side/us02_model.sbatch      (after us01; not an array)
# Tune:    params.undersequencing.coverage_threshold, threshold_sweep

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
OUT="$PATHS_UNDERSEQUENCING"
mapfile -t NPO < <(grep -v '^\s*$' "$PATHS_RUN_ACCESSIONS" | sed "s|^|$OUT/nonpareil/|; s|$|.npo|")
require "${NPO[@]}"

activate "$ENVS_ASSEMBLY"     # Nonpareil R package
Rscript "$PD_REPO/scripts/R/nonpareil_summary.R" "$OUT/nonpareil_all.tsv" "${NPO[@]}"
mapfile -t PILOT < <(run_samples | sed "s|.*|$PATHS_ASSEMBLY/&/nonpareil/&.npo|")
require "${PILOT[@]}" "$PATHS_TABLES/04_assembly.tsv"
Rscript "$PD_REPO/scripts/R/nonpareil_summary.R" "$OUT/pilot_hostfree.tsv" "${PILOT[@]}"
deactivate

# pilot: run, coverage on raw R1 (us01), coverage on host-depleted R1 (stage 4d), assembled %
awk -F'\t' 'FNR==1{f++; next}
    f==1{raw[$1]=$3} f==2{hf[$1]=$3} f==3{asm[$1]=$5}
    END{print "run\tcoverage_raw\tcoverage_hostfree\tassembled_fraction_pct";
        for (r in hf) printf "%s\t%.4f\t%.4f\t%s\n", r, raw[r], hf[r], asm[r]}' \
    "$OUT/nonpareil_all.tsv" "$OUT/pilot_hostfree.tsv" "$PATHS_TABLES/04_assembly.tsv" \
    | sort > "$OUT/pilot.tsv"

activate "$ENVS_MODEL"
Rscript "$PD_REPO/scripts/R/undersequencing.R" "$OUT/nonpareil_all.tsv" "$PATHS_SAMPLES" \
    "$PARAMS_STATS_SITE_GROUPS" "$PARAMS_UNDERSEQUENCING_COVERAGE_THRESHOLD" \
    "$PARAMS_UNDERSEQUENCING_THRESHOLD_SWEEP" "$OUT/pilot.tsv" "$OUT"
column -t "$OUT/pilot_calibration.tsv"
