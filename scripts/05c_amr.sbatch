#!/bin/bash
#SBATCH --job-name=pd_amr
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 05c_amr.sbatch -- stage 5, one run per array task: AMRFinderPlus (--plus) and RGI with
# CARD on the WHOLE assembly (every contig, not just plasmid/virus ones).
# RGI --local reads ./localDB, so it runs in its output dir with localDB linked in.
# Output: $paths.contigs/<run>/amrfinder.tsv, $paths.contigs/<run>/rgi/<run>.rgi.txt
# Submit:  sbatch --array=1-$(grep -c . metadata/pilot_samples.txt) scripts/05c_amr.sbatch
# Tune:    params.contigs.amrfinder_args

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
SRR=$(task_sample)
T=$PARAMS_CONTIGS_THREADS
CONTIGS="$PATHS_ASSEMBLY/$SRR/${SRR}.contigs.fa"
DIR="$PATHS_CONTIGS/$SRR"
require "$CONTIGS"
mkdir -p "$DIR/rgi"
echo "[$(date)] $SRR: AMRFinderPlus + RGI"

activate "$ENVS_AMRFINDERPLUS"
# shellcheck disable=SC2086
amrfinder --nucleotide "$CONTIGS" $PARAMS_CONTIGS_AMRFINDER_ARGS --database "$DB_AMRFINDERPLUS" \
    --threads "$T" --name "$SRR" --output "$DIR/amrfinder.tsv"
deactivate
echo "AMRFinderPlus hits: $(( $(wc -l < "$DIR/amrfinder.tsv") - 1 ))"

activate "$ENVS_RGI"
cd "$DIR/rgi"
ln -sfn "$DB_CARD_DIR/localDB" localDB
rgi main --input_sequence "$CONTIGS" --output_file "$DIR/rgi/$SRR.rgi" --input_type contig \
    --alignment_tool DIAMOND --threads "$T" --local --clean
deactivate
echo "RGI hits: $(( $(wc -l < "$DIR/rgi/$SRR.rgi.txt") - 1 ))"
echo "[$(date)] $SRR: done"
