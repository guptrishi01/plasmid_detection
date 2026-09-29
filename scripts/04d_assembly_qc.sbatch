#!/bin/bash
#SBATCH --job-name=pd_asm_qc
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 04d_assembly_qc.sbatch -- stage 4, one run per array task: QUAST on the MEGAHIT contigs,
# Nonpareil (read redundancy / coverage estimate) on the host-depleted R1.
# Output: $paths.assembly/<run>/quast/, $paths.assembly/<run>/nonpareil/<run>.npo
# Submit:  sbatch --array=1-$(grep -c . metadata/pilot_samples.txt) scripts/04d_assembly_qc.sbatch
# Tune:    params.assembly.nonpareil_args

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_ASSEMBLY"
SRR=$(task_sample)
T=$PARAMS_ASSEMBLY_QC_THREADS
DIR="$PATHS_ASSEMBLY/$SRR"
require "$DIR/${SRR}.contigs.fa" "$PATHS_HOSTFREE/${SRR}_1.fastq.gz"
echo "[$(date)] $SRR: QUAST + Nonpareil"

quast.py "$DIR/${SRR}.contigs.fa" -o "$DIR/quast" -t "$T" \
    --min-contig "$PARAMS_ASSEMBLY_MIN_CONTIG_LEN" --labels "$SRR"

# Nonpareil reads one uncompressed file
mkdir -p "$DIR/nonpareil"
mkdir -p "$PATHS_WORK"
mkdir -p "$PATHS_WORK"
TMP=$(mktemp -d "$PATHS_WORK/nonpareil.XXXX"); trap 'rm -rf "$TMP"' EXIT
zcat "$PATHS_HOSTFREE/${SRR}_1.fastq.gz" > "$TMP/${SRR}_1.fastq"
# shellcheck disable=SC2086
nonpareil -s "$TMP/${SRR}_1.fastq" -f fastq $PARAMS_ASSEMBLY_NONPAREIL_ARGS \
    -b "$DIR/nonpareil/$SRR" -t "$T"
echo "[$(date)] $SRR: done"
