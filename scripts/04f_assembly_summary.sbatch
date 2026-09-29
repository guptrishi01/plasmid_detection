#!/bin/bash
#SBATCH --job-name=pd_asm_summary
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 04f_assembly_summary.sbatch -- stage 4 table: contigs, total bp, N50 (QUAST, contigs >=
# params.assembly.min_contig_len) and assembled fraction (share of the run's own reads that
# Bowtie2 aligned back to its own assembly) per sample.
# Output: $paths.tables/04_assembly.tsv

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
pdpy plasmid_detection.tables assembly --samples "$PATHS_RUN_SAMPLES" \
    --assembly-dir "$PATHS_ASSEMBLY" --out "$PATHS_TABLES/04_assembly.tsv"
column -t "$PATHS_TABLES/04_assembly.tsv"
