#!/bin/bash
#SBATCH --job-name=pd_qc
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 02_read_qc.sbatch -- stage 2, one run per array task:
#   FastQC on raw reads -> fastp -> FastQC on trimmed reads
# Submit:  sbatch --array=1-$(grep -c . metadata/pilot_samples.txt) scripts/02_read_qc.sbatch
# Then:    sbatch scripts/02b_qc_summary.sbatch   (MultiQC over everything + table)
# Tune:    params.qc.fastp_args

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_QC"
SRR=$(task_sample)
T=$PARAMS_QC_THREADS

RAW1="$PATHS_RAW_READS/${SRR}_1.fastq.gz"; RAW2="$PATHS_RAW_READS/${SRR}_2.fastq.gz"
TRIM1="$PATHS_TRIMMED/${SRR}_1.fastq.gz"; TRIM2="$PATHS_TRIMMED/${SRR}_2.fastq.gz"
require "$RAW1" "$RAW2"
mkdir -p "$PATHS_TRIMMED" "$PATHS_QC"/{fastqc_raw,fastp,fastqc_trimmed}
echo "[$(date)] $SRR: stage 2 read QC"

fastqc --quiet --threads "$T" -o "$PATHS_QC/fastqc_raw" "$RAW1" "$RAW2"

# shellcheck disable=SC2086  # fastp_args is a list of options
fastp -i "$RAW1" -I "$RAW2" -o "$TRIM1.tmp.gz" -O "$TRIM2.tmp.gz" \
    --thread "$T" $PARAMS_QC_FASTP_ARGS \
    --json "$PATHS_QC/fastp/${SRR}.fastp.json" --html "$PATHS_QC/fastp/${SRR}.fastp.html" \
    --report_title "$SRR" 2> "$PATHS_QC/fastp/${SRR}.fastp.log"
mv "$TRIM1.tmp.gz" "$TRIM1"; mv "$TRIM2.tmp.gz" "$TRIM2"

fastqc --quiet --threads "$T" -o "$PATHS_QC/fastqc_trimmed" "$TRIM1" "$TRIM2"
echo "[$(date)] $SRR: done"
