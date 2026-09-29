#!/bin/bash
#SBATCH --job-name=pd_host_summary
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 03b_host_summary.sbatch -- stage 3 summary: MultiQC of the Bowtie2 logs, and the table
# of raw, trimmed, human, retained pairs per sample.
# Output: $paths.host/multiqc/, $paths.tables/03_host_depletion.tsv

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_QC"   # MultiQC lives in pd-qc
mkdir -p "$PATHS_TABLES"

multiqc --force --outdir "$PATHS_HOST/multiqc" --filename host_depletion "$PATHS_HOST/bowtie2_logs"

pdpy plasmid_detection.tables host --samples "$PATHS_RUN_SAMPLES" \
    --fastp-dir "$PATHS_QC/fastp" --counts-dir "$PATHS_HOST" --out "$PATHS_TABLES/03_host_depletion.tsv"
column -t "$PATHS_TABLES/03_host_depletion.tsv"
