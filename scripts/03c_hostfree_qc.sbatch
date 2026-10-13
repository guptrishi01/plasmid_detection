#!/bin/bash
#SBATCH --job-name=pd_hostfree_qc
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 03c_hostfree_qc.sbatch -- stage 3 add-on: FastQC on the host-depleted reads (the stage 4
# input), then one MultiQC over those + the Bowtie2 logs. Not in the original spec; added
# 2026-10-12 so the reads that go into assembly are QC'd too.
# Output: $paths.host/fastqc_hostfree/, $paths.host/multiqc/hostfree_qc.html
# Submit:  sbatch scripts/03c_hostfree_qc.sbatch      (not an array; loops over runs)

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_QC"
OUT="$PATHS_HOST/fastqc_hostfree"
mkdir -p "$OUT"
FILES=()
while read -r s; do
    for r in 1 2; do
        f="$PATHS_HOSTFREE/${s}_${r}.fastq.gz"
        require "$f"
        FILES+=("$f")
    done
done < <(run_samples)
echo "[$(date)] FastQC on ${#FILES[@]} host-free files"
fastqc --quiet --threads "$PARAMS_QC_THREADS" -o "$OUT" "${FILES[@]}"

# --dirs: keeps sample names unique if more FastQC sets are added to this report later
multiqc --force --dirs --dirs-depth 1 --outdir "$PATHS_HOST/multiqc" --filename hostfree_qc \
    "$PATHS_HOST/bowtie2_logs" "$OUT"
echo "[$(date)] done: $PATHS_HOST/multiqc/hostfree_qc.html"
