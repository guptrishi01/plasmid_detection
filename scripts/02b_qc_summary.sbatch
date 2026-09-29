#!/bin/bash
#SBATCH --job-name=pd_qc_summary
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 02b_qc_summary.sbatch -- stage 2 summary: MultiQC over raw FastQC + fastp + trimmed
# FastQC, and the per-sample table (raw pairs, trimmed pairs, Q30, duplication, insert size).
# Output: $paths.qc/multiqc/, $paths.tables/02_read_qc.tsv

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_QC"
mkdir -p "$PATHS_TABLES"

multiqc --force --outdir "$PATHS_QC/multiqc" --filename read_qc \
    "$PATHS_QC/fastqc_raw" "$PATHS_QC/fastp" "$PATHS_QC/fastqc_trimmed"

pdpy plasmid_detection.tables qc --samples "$PATHS_RUN_SAMPLES" \
    --fastp-dir "$PATHS_QC/fastp" --out "$PATHS_TABLES/02_read_qc.tsv"
column -t "$PATHS_TABLES/02_read_qc.tsv"
