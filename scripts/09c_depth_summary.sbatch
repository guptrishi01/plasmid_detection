#!/bin/bash
#SBATCH --job-name=pd_depth_summary
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 09c_depth_summary.sbatch -- stage 9: add the "full" point for each depth run (stage 4
# MEGAHIT + stage 5 geNomad, same params), collect all points, plot plasmid yield vs depth.
# Output: $paths.tables/09_depth_experiment.tsv, $paths.rollup/depth/plasmid_yield_vs_depth.{pdf,png}

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
for SRR in $PARAMS_ROLLUP_DEPTH_SAMPLES; do
    OUT="$PATHS_ROLLUP/depth/${SRR}_full"; mkdir -p "$OUT"
    TOTAL=$(awk -F'\t' 'NR==2{print $4}' "$PATHS_HOST/${SRR}.counts.tsv")
    pdpy plasmid_detection.tables depth-yield --run "$SRR" --level full --pairs "$TOTAL" \
        --contigs "$PATHS_ASSEMBLY/$SRR/${SRR}.contigs.fa" --genomad-dir "$PATHS_CONTIGS/$SRR/genomad" \
        --out "$OUT/yield.tsv"
done
mkdir -p "$PATHS_TABLES"
pdpy plasmid_detection.tables depth --depth-dir "$PATHS_ROLLUP/depth" --out "$PATHS_TABLES/09_depth_experiment.tsv"
column -t "$PATHS_TABLES/09_depth_experiment.tsv"

activate "$ENVS_STATS"
Rscript "$PD_REPO/scripts/R/depth_plot.R" "$PATHS_TABLES/09_depth_experiment.tsv" \
    "$PATHS_ROLLUP/depth/plasmid_yield_vs_depth"
