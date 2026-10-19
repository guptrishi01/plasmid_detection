#!/bin/bash
#SBATCH --job-name=pd_us_nonpareil
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# us01_nonpareil_all.sbatch -- side analysis (not in the spec), one run per array task over
# ALL 92 runs (paths.run_accessions, not paths.run_samples): Nonpareil (same mode as stage
# 4d, params.assembly.nonpareil_args) on the raw R1 reads. Raw, because stages 2-3 have only
# run on the pilot; us02 compares raw vs host-depleted Nonpareil on the 4 pilot runs to show
# what that costs.
# Output: $paths.undersequencing/nonpareil/<run>.{npo,npa,npc,npl}
# Submit:  sbatch --array=1-$(grep -c . metadata/srp121672_wgs_run_accessions.txt) scripts/side/us01_nonpareil_all.sbatch

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_ASSEMBLY"
SRR=$(grep -v '^\s*$' "$PATHS_RUN_ACCESSIONS" | sed -n "${SLURM_ARRAY_TASK_ID:?submit as an array}p")
[[ -n "$SRR" ]] || { echo "no run for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
IN="$PATHS_RAW_READS/${SRR}_1.fastq.gz"
require "$IN"
OUT="$PATHS_UNDERSEQUENCING/nonpareil"
mkdir -p "$OUT" "$PATHS_WORK"
rm -f "$OUT/$SRR".np?
TMP=$(mktemp -d "$PATHS_WORK/us_nonpareil.XXXX"); trap 'rm -rf "$TMP"' EXIT
echo "[$(date)] $SRR: Nonpareil on raw R1"
zcat "$IN" > "$TMP/${SRR}_1.fastq"     # Nonpareil reads one uncompressed file
# shellcheck disable=SC2086
nonpareil -s "$TMP/${SRR}_1.fastq" -f fastq $PARAMS_ASSEMBLY_NONPAREIL_ARGS \
    -b "$OUT/$SRR" -t "$PARAMS_UNDERSEQUENCING_THREADS"
require "$OUT/$SRR.npo"
echo "[$(date)] $SRR: done"
