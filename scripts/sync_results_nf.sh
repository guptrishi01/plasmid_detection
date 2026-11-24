#!/bin/bash
# sync_results_nf.sh -- copy the shareable part of a Nextflow --outdir (scratch) into the
# repo for lab meetings. A 92-run --outdir is ~20 GB; scripts/sync_results_nf.filter keeps
# tables, summaries and small per-run reports (~130 MB) and leaves sequences, BAMs, depth
# tables and raw tool dumps on scratch. Each run mirrors the destination (--delete).
#     bash scripts/sync_results_nf.sh [outdir] [dest]
# Defaults: /scratch/rgupta25/plasmid_detection/results_nf -> results/full_92
set -euo pipefail
PD_REPO=${PD_REPO:-$HOME/gibas/plasmid_detection}
SRC=${1:-/scratch/rgupta25/plasmid_detection/results_nf}
DEST=${2:-$PD_REPO/results/full_92}
[[ -d "$SRC/tables" ]] || { echo "no $SRC/tables: not a pipeline --outdir" >&2; exit 1; }
mkdir -p "$DEST"
rsync -a --delete --prune-empty-dirs --max-size=20M --filter="merge $PD_REPO/scripts/sync_results_nf.filter" \
    "$SRC/" "$DEST/"
echo "$DEST: $(find "$DEST" -type f | wc -l) files, $(du -sh "$DEST" | cut -f1)"
big=$(find "$SRC/tables" $([[ -d "$SRC/09_rollup" ]] && echo "$SRC/09_rollup") -maxdepth 1 -type f -size +20M | wc -l)
(( big == 0 )) || echo "note: $big table(s) > 20 MB not copied" >&2
