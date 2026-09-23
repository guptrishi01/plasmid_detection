#!/bin/bash
#SBATCH --job-name=pd_manifest
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 01b_collect_manifest.sbatch -- combine the per-run stage-1 manifests into
# paths.raw_manifest and check all 92 runs are present, PASS, and read-only.
#     sbatch --dependency=afterany:<01 array job id> scripts/01b_collect_manifest.sbatch

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"

expected=$(grep -c . "$PATHS_RUN_ACCESSIONS")
problems=0
{
    head -n 1 "$(ls "$PATHS_MANIFEST_DIR"/*.tsv | head -1)"
    while read -r srr; do
        f="$PATHS_MANIFEST_DIR/$srr.tsv"
        if [[ -f "$f" ]]; then tail -n +2 "$f"; fi
    done < "$PATHS_RUN_ACCESSIONS"
} > "$PATHS_RAW_MANIFEST"

while read -r srr; do
    f="$PATHS_MANIFEST_DIR/$srr.tsv"
    if [[ ! -f "$f" ]]; then echo "MISSING manifest: $srr"; problems=$((problems + 1)); continue; fi
    grep -q $'\tPASS$' "$f" || { echo "NOT PASS: $srr"; problems=$((problems + 1)); }
    for r in 1 2; do
        fq="$PATHS_RAW_READS/${srr}_$r.fastq.gz"
        [[ -f "$fq" ]] || { echo "MISSING file: $fq"; problems=$((problems + 1)); continue; }
        [[ -w "$fq" ]] && { echo "WRITABLE: $fq"; problems=$((problems + 1)); }
    done
done < "$PATHS_RUN_ACCESSIONS"

rows=$(( $(wc -l < "$PATHS_RAW_MANIFEST") - 1 ))
echo "manifest: $rows/$expected runs -> $PATHS_RAW_MANIFEST"
echo "problems: $problems"
(( problems == 0 && rows == expected ))
