#!/bin/bash
#SBATCH --job-name=pd_genomad
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 05a_genomad.sbatch -- stage 5, one run per array task: geNomad end-to-end on the whole
# assembly -> plasmid and virus contigs.
# Output: $paths.contigs/<run>/genomad/*_summary/{*_plasmid.fna,*_plasmid_summary.tsv,
#         *_virus.fna,*_virus_summary.tsv}
# PlasLLM (requested for the plasmid calls) is not installed yet: waiting on its source.
# It belongs here, next to geNomad, on the same contigs.
# Submit:  sbatch --array=1-$(grep -c . metadata/pilot_samples.txt) scripts/05a_genomad.sbatch
# Tune:    params.contigs.genomad_args

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_GENOMAD"
SRR=$(task_sample)
CONTIGS="$PATHS_ASSEMBLY/$SRR/${SRR}.contigs.fa"
OUT="$PATHS_CONTIGS/$SRR/genomad"
require "$CONTIGS"
rm -rf "$OUT"; mkdir -p "$OUT"
echo "[$(date)] $SRR: geNomad end-to-end"

# shellcheck disable=SC2086
genomad end-to-end $PARAMS_CONTIGS_GENOMAD_ARGS --threads "$PARAMS_CONTIGS_THREADS" \
    "$CONTIGS" "$OUT" "$DB_GENOMAD"

for kind in plasmid virus; do
    f=$(ls "$OUT"/*_summary/*_${kind}_summary.tsv)
    echo "$kind contigs: $(( $(wc -l < "$f") - 1 ))"
done
echo "[$(date)] $SRR: done"
