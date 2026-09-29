#!/bin/bash
#SBATCH --job-name=pd_depth
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 09b_depth_experiment.sbatch -- stage 9 depth experiment, one (run, depth) per array task:
# subsample the host-depleted pairs with seqkit to N million pairs, assemble (MEGAHIT, same
# params as stage 4), geNomad, record plasmid yield. Tasks enumerate
#   params.rollup.depth_samples x params.rollup.depth_levels_million_pairs
# The "full" point reuses the stage 4/5 assembly + geNomad (same params) -- see 09c.
# A level above a run's available pairs is skipped and logged.
# Submit:  sbatch --array=1-8 scripts/09b_depth_experiment.sbatch   (2 runs x 4 levels)

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
read -ra RUNS <<< "$PARAMS_ROLLUP_DEPTH_SAMPLES"
read -ra LEVELS <<< "$PARAMS_ROLLUP_DEPTH_LEVELS_MILLION_PAIRS"
i=$(( SLURM_ARRAY_TASK_ID - 1 ))
(( i < ${#RUNS[@]} * ${#LEVELS[@]} )) || { echo "task $SLURM_ARRAY_TASK_ID > ${#RUNS[@]}x${#LEVELS[@]}" >&2; exit 1; }
SRR=${RUNS[i / ${#LEVELS[@]}]}
LEVEL=${LEVELS[i % ${#LEVELS[@]}]}
N=$(( LEVEL * 1000000 ))
IN1="$PATHS_HOSTFREE/${SRR}_1.fastq.gz"; IN2="$PATHS_HOSTFREE/${SRR}_2.fastq.gz"
require "$IN1" "$IN2" "$PATHS_HOST/${SRR}.counts.tsv"
TOTAL=$(awk -F'\t' 'NR==2{print $4}' "$PATHS_HOST/${SRR}.counts.tsv")
OUT="$PATHS_ROLLUP/depth/${SRR}_${LEVEL}M"
rm -rf "$OUT"; mkdir -p "$OUT"
echo "[$(date)] $SRR at ${LEVEL}M pairs (available: $TOTAL)"
if (( N > TOTAL )); then
    echo "skipped: ${LEVEL}M > $TOTAL available pairs" | tee "$OUT/SKIPPED"
    exit 0
fi

# Same seed + same record count in R1 and R2 -> seqkit keeps the same records in both;
# sample a little over N, then take exactly N
activate "$ENVS_ASSEMBLY"
P=$(awk -v n="$N" -v t="$TOTAL" 'BEGIN{p=n*1.02/t; if(p>1)p=1; printf "%.6f", p}')
for r in 1 2; do
    in="$PATHS_HOSTFREE/${SRR}_${r}.fastq.gz"
    seqkit sample -p "$P" -s "$PARAMS_ROLLUP_DEPTH_SEED" "$in" 2>/dev/null \
        | seqkit head -n "$N" -o "$OUT/sub_${r}.fastq.gz"
done
n1=$(seqkit stats -T "$OUT/sub_1.fastq.gz" | awk 'NR==2{print $4}')
n2=$(seqkit stats -T "$OUT/sub_2.fastq.gz" | awk 'NR==2{print $4}')
last1=$(seqkit seq -n -i "$OUT/sub_1.fastq.gz" | tail -1)
last2=$(seqkit seq -n -i "$OUT/sub_2.fastq.gz" | tail -1)
[[ "$n1" == "$N" && "$n2" == "$N" && "$last1" == "$last2" ]] \
    || { echo "subsample out of sync: R1=$n1 ($last1) R2=$n2 ($last2), want $N" >&2; exit 1; }

MEM_BYTES=$(( ${SLURM_MEM_PER_NODE:-65536} * 1024 * 1024 * 9 / 10 ))
# shellcheck disable=SC2086
megahit -1 "$OUT/sub_1.fastq.gz" -2 "$OUT/sub_2.fastq.gz" -o "$OUT/megahit" \
    -t "$PARAMS_ROLLUP_THREADS" -m "$MEM_BYTES" \
    --min-contig-len "$PARAMS_ASSEMBLY_MIN_CONTIG_LEN" $PARAMS_ASSEMBLY_MEGAHIT_ARGS
awk -v s="$SRR" '/^>/{sub(/^>/, ">" s "_")} {print}' "$OUT/megahit/final.contigs.fa" > "$OUT/contigs.fa"
deactivate

activate "$ENVS_GENOMAD"
# shellcheck disable=SC2086
genomad end-to-end $PARAMS_CONTIGS_GENOMAD_ARGS --threads "$PARAMS_ROLLUP_THREADS" \
    "$OUT/contigs.fa" "$OUT/genomad" "$DB_GENOMAD"
deactivate

pdpy plasmid_detection.tables depth-yield --run "$SRR" --level "${LEVEL}M" --pairs "$N" \
    --contigs "$OUT/contigs.fa" --genomad-dir "$OUT/genomad" --out "$OUT/yield.tsv"
rm -rf "$OUT/sub_1.fastq.gz" "$OUT/sub_2.fastq.gz" "$OUT/megahit/intermediate_contigs"
cat "$OUT/yield.tsv"
echo "[$(date)] done"
