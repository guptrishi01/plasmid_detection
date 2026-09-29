#!/bin/bash
#SBATCH --job-name=pd_megahit
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 04a_megahit.sbatch -- stage 4, one run per array task: MEGAHIT on the host-depleted
# pairs. Contig headers are prefixed with the run id (">SRR6231181_k127_42"). The assembly
# graph is kept: intermediate_contigs/k<K>.contigs.fa for the largest K MEGAHIT reached,
# converted with megahit_toolkit contig2fastg (final.contigs.fa merges contigs from all k,
# so it is not one graph). SCAPP (04e) reads that FASTG.
# Submit:  sbatch --array=1-$(grep -c . metadata/pilot_samples.txt) scripts/04a_megahit.sbatch
# Tune:    params.assembly.megahit_args, params.assembly.min_contig_len

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_ASSEMBLY"
SRR=$(task_sample)
T=$PARAMS_ASSEMBLY_THREADS

IN1="$PATHS_HOSTFREE/${SRR}_1.fastq.gz"; IN2="$PATHS_HOSTFREE/${SRR}_2.fastq.gz"
require "$IN1" "$IN2"
OUT="$PATHS_ASSEMBLY/$SRR"
mkdir -p "$OUT"
rm -rf "$OUT/megahit"   # MEGAHIT refuses to write into an existing directory
echo "[$(date)] $SRR: MEGAHIT"

# MEGAHIT sizes itself from the node's total RAM, not the Slurm allocation: pass 90% of
# the allocation in bytes instead
MEM_BYTES=$(( ${SLURM_MEM_PER_NODE:-65536} * 1024 * 1024 * 9 / 10 ))
# shellcheck disable=SC2086
megahit -1 "$IN1" -2 "$IN2" -o "$OUT/megahit" -t "$T" -m "$MEM_BYTES" \
    --min-contig-len "$PARAMS_ASSEMBLY_MIN_CONTIG_LEN" $PARAMS_ASSEMBLY_MEGAHIT_ARGS

# Contigs with sample-prefixed headers (first word only; MEGAHIT's flag/multi/len kept
# after a space so tools that read the first word see "<SRR>_k<K>_<n>")
awk -v s="$SRR" '/^>/{sub(/^>/, ">" s "_")} {print}' "$OUT/megahit/final.contigs.fa" \
    > "$OUT/${SRR}.contigs.fa"

# Assembly graph from the largest-k contig set
K=$(ls "$OUT/megahit/intermediate_contigs" | sed -n 's/^k\([0-9]*\)\.contigs\.fa$/\1/p' | sort -n | tail -1)
[[ -n "$K" ]] || { echo "no intermediate_contigs/k*.contigs.fa to build the graph from" >&2; exit 1; }
megahit_toolkit contig2fastg "$K" "$OUT/megahit/intermediate_contigs/k${K}.contigs.fa" \
    > "$OUT/${SRR}.k${K}.fastg"
echo "$K" > "$OUT/graph_k.txt"

grep -c '^>' "$OUT/${SRR}.contigs.fa" | xargs echo "contigs >= $PARAMS_ASSEMBLY_MIN_CONTIG_LEN bp:"
echo "graph: $OUT/${SRR}.k${K}.fastg"
echo "[$(date)] $SRR: done"
