#!/bin/bash
#SBATCH --job-name=pd_binning
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 06a_binning.sbatch -- stage 6, one ASSEMBLY per array task: jgi_summarize_bam_contig_depths
# over ALL stage-4 BAMs on the assembly (its own reads + every system mate), then MetaBAT2.
# Output: $paths.binning/<run>/depth.txt, $paths.binning/<run>/bins/<run>.bin.<n>.fa
# Submit:  sbatch --array=1-$(grep -c . metadata/pilot_samples.txt) scripts/06a_binning.sbatch
# Tune:    params.binning.metabat2_min_contig (>= 1500), params.binning.metabat2_seed

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_BINNING"
SRR=$(task_sample)
CONTIGS="$PATHS_ASSEMBLY/$SRR/${SRR}.contigs.fa"
BAMS=("$PATHS_ASSEMBLY/$SRR"/map/*.bam)
OUT="$PATHS_BINNING/$SRR"
require "$CONTIGS" "${BAMS[@]}"
rm -rf "$OUT"; mkdir -p "$OUT/bins"
echo "[$(date)] $SRR: depths over ${#BAMS[@]} BAMs + MetaBAT2"

jgi_summarize_bam_contig_depths --outputDepth "$OUT/depth.txt" "${BAMS[@]}"
metabat2 --inFile "$CONTIGS" --abdFile "$OUT/depth.txt" --outFile "$OUT/bins/${SRR}.bin" \
    --minContig "$PARAMS_BINNING_METABAT2_MIN_CONTIG" --seed "$PARAMS_BINNING_METABAT2_SEED" \
    --numThreads "$PARAMS_BINNING_THREADS"
echo "bins: $(ls "$OUT/bins" | grep -c '\.fa$')"
echo "[$(date)] $SRR: done"
