#!/bin/bash
# sync_results.sh -- copy the shareable part of $paths.results (scratch) into results/ in
# the repo, so reports and tables can be committed and shown. Run on a login node:
#     bash scripts/sync_results.sh [stage_dir ...]      (default: every stage dir)
# e.g. bash scripts/sync_results.sh 02_qc 03_host tables
#
# Copies reports, logs, tables, QC summaries. Never copies reads, contigs, graphs,
# alignments, or tool working dirs, and skips any file over 20 MB (GitHub rejects files
# over 100 MB; the assembly outputs themselves stay on scratch). Skipped large files are
# listed so nothing is dropped silently.

set -euo pipefail
PD_REPO=${PD_REPO:-$HOME/gibas/plasmid_detection}
eval "$(PYTHONPATH="$PD_REPO/src" "$PD_REPO/conda/envs/pd-core/bin/python" -m plasmid_detection.config shell)"
SRC=$PATHS_RESULTS
DEST=$PD_REPO/results
MAX=20M

if (( $# )); then STAGES=("$@"); else mapfile -t STAGES < <(ls "$SRC"); fi
mkdir -p "$DEST"
for stage in "${STAGES[@]}"; do
    [[ -d "$SRC/$stage" ]] || { echo "no $SRC/$stage" >&2; exit 1; }
    rsync -a --prune-empty-dirs --max-size="$MAX" \
        --exclude 'megahit/' --exclude 'bt2_index/' --exclude 'intermediate_files/' \
        --exclude '*.bam' --exclude '*.bai' --exclude '*.sam' \
        --exclude '*.fa' --exclude '*.fna' --exclude '*.fasta' --exclude '*.fastg' \
        --exclude '*.fastq' --exclude '*.fastq.gz' --exclude '*.fq' --exclude '*.fq.gz' \
        --exclude '*.bt2' --exclude '*.bt2l' --exclude '*.npa' \
        "$SRC/$stage/" "$DEST/$stage/"
    big=$(find "$SRC/$stage" -type f -size +"$MAX" ! -name '*.bam' ! -name '*.fa' ! -name '*.fasta' \
          ! -name '*.fastg' ! -name '*.fastq.gz' ! -path '*/megahit/*' ! -path '*/intermediate_files/*' | wc -l)
    echo "$stage: $(find "$DEST/$stage" -type f | wc -l) files, $(du -sh "$DEST/$stage" | cut -f1)" \
         "(skipped $big other files > $MAX)"
done
