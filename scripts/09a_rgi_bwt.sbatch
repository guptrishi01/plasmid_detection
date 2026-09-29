#!/bin/bash
#SBATCH --job-name=pd_rgi_bwt
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 09a_rgi_bwt.sbatch -- stage 9, one run per array task: RGI bwt (read level, Bowtie2
# aligner, CARD canonical references loaded by 00_setup) on the host-depleted pairs, for
# the read-level vs contig-level comparison (the paper quantified ARGs from reads).
# Output: $paths.rollup/rgi_bwt/<run>/<run>.gene_mapping_data.txt
# Submit:  sbatch --array=1-$(grep -c . metadata/pilot_samples.txt) scripts/09a_rgi_bwt.sbatch

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_RGI"
SRR=$(task_sample)
IN1="$PATHS_HOSTFREE/${SRR}_1.fastq.gz"; IN2="$PATHS_HOSTFREE/${SRR}_2.fastq.gz"
require "$IN1" "$IN2"
OUT="$PATHS_ROLLUP/rgi_bwt/$SRR"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
ln -sfn "$DB_CARD_DIR/localDB" localDB
echo "[$(date)] $SRR: RGI bwt"
rgi bwt --read_one "$IN1" --read_two "$IN2" --aligner bowtie2 --threads "$PARAMS_ROLLUP_THREADS" \
    --output_file "$OUT/$SRR" --local --clean
echo "ARO terms with mapped reads: $(( $(wc -l < "$OUT/$SRR.gene_mapping_data.txt") - 1 ))"
echo "[$(date)] $SRR: done"
