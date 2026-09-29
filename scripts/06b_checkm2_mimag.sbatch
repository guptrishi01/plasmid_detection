#!/bin/bash
#SBATCH --job-name=pd_checkm2
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 06b_checkm2_mimag.sbatch -- stage 6, one run per array task: CheckM2 on the MetaBAT2 bins,
# then MIMAG tiers from completeness/contamination (params.binning.mimag_*):
#   high   >= 90% complete, < 5% contamination
#   medium >= 50% complete, < 10% contamination
#   low    everything else
# MIMAG "high-quality draft" also requires 23S/16S/5S rRNA and >= 18 tRNAs; CheckM2 does
# not report those, so "high" here is completeness/contamination only (see README).
# Medium and up are copied to $paths.binning/pass/ for stage 7.
# Output: $paths.binning/<run>/checkm2/quality_report.tsv, $paths.binning/<run>/mimag.tsv
# Submit:  sbatch --array=1-$(grep -c . metadata/pilot_samples.txt) scripts/06b_checkm2_mimag.sbatch

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_CHECKM2"
SRR=$(task_sample)
OUT="$PATHS_BINNING/$SRR"
[[ -n "$(ls "$OUT"/bins/*.fa 2>/dev/null)" ]] || { echo "no bins in $OUT/bins (run 06a)" >&2; exit 1; }
echo "[$(date)] $SRR: CheckM2 on $(ls "$OUT"/bins/*.fa | wc -l) bins"

checkm2 predict --input "$OUT/bins" -x fa --output-directory "$OUT/checkm2" \
    --threads "$PARAMS_BINNING_THREADS" --force
deactivate

pdpy plasmid_detection.tables mimag --report "$OUT/checkm2/quality_report.tsv" --run "$SRR" \
    --high-comp "$PARAMS_BINNING_MIMAG_HIGH_COMPLETENESS" --high-cont "$PARAMS_BINNING_MIMAG_HIGH_CONTAMINATION" \
    --medium-comp "$PARAMS_BINNING_MIMAG_MEDIUM_COMPLETENESS" --medium-cont "$PARAMS_BINNING_MIMAG_MEDIUM_CONTAMINATION" \
    --out "$OUT/mimag.tsv"

mkdir -p "$PATHS_BINNING/pass"
rm -f "$PATHS_BINNING/pass/${SRR}".bin.*.fa
awk -F'\t' 'NR>1 && ($5=="high" || $5=="medium") {print $2}' "$OUT/mimag.tsv" | while read -r bin; do
    cp "$OUT/bins/${bin}.fa" "$PATHS_BINNING/pass/"
done
echo "tiers:"; cut -f5 "$OUT/mimag.tsv" | tail -n +2 | sort | uniq -c
echo "[$(date)] $SRR: done"
