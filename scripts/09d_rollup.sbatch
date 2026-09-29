#!/bin/bash
#SBATCH --job-name=pd_rollup
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 09d_rollup.sbatch -- stage 9 rollup:
#   plasmid summary per sample (geNomad, mob_typer mobility, SCAPP, AMR on plasmids)
#   resistome matrices: gene x sample read counts over AMR hits (own reads on own
#     assembly), all / mobile (plasmid+virus) / chromosome -- for AMRFinderPlus and RGI
#   comparison with the paper's read-level approach: RGI bwt (09a) vs RGI contig level
# Output: $paths.tables/09_plasmid_summary.tsv, $paths.rollup/resistome_*.tsv,
#         $paths.tables/09_rgi_read_vs_contig.tsv

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
mkdir -p "$PATHS_TABLES" "$PATHS_ROLLUP"
pdpy plasmid_detection.tables plasmid-summary --samples "$PATHS_RUN_SAMPLES" \
    --contigs-dir "$PATHS_CONTIGS" --assembly-dir "$PATHS_ASSEMBLY" --out "$PATHS_TABLES/09_plasmid_summary.tsv"
for src in amrfinder rgi; do
    pdpy plasmid_detection.tables resistome --samples "$PATHS_RUN_SAMPLES" \
        --contigs-dir "$PATHS_CONTIGS" --source "$src" --out-dir "$PATHS_ROLLUP"
done
pdpy plasmid_detection.tables rgi-compare --samples "$PATHS_RUN_SAMPLES" \
    --contigs-dir "$PATHS_CONTIGS" --bwt-dir "$PATHS_ROLLUP/rgi_bwt" \
    --out "$PATHS_TABLES/09_rgi_read_vs_contig.tsv"
column -t "$PATHS_TABLES/09_plasmid_summary.tsv"
echo "RGI read vs contig:"; cut -f1,6 "$PATHS_TABLES/09_rgi_read_vs_contig.tsv" | tail -n +2 | sort | uniq -c
