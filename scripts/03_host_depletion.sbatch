#!/bin/bash
#SBATCH --job-name=pd_host
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 03_host_depletion.sbatch -- stage 3, one run per array task: Bowtie2 of the trimmed
# pairs against GRCh38.p14; keep only pairs where NEITHER mate maps (SAM flags 4+8 both
# set, i.e. samtools -f 12). Note bowtie2 --un-conc would also keep pairs where one mate
# maps, which is not what the spec asks for.
# Submit:  sbatch --array=1-$(grep -c . metadata/pilot_samples.txt) scripts/03_host_depletion.sbatch
# Then:    sbatch scripts/03b_host_summary.sbatch
# Tune:    params.host.bowtie2_args

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_HOST"
SRR=$(task_sample)
T=$PARAMS_HOST_THREADS

IN1="$PATHS_TRIMMED/${SRR}_1.fastq.gz"; IN2="$PATHS_TRIMMED/${SRR}_2.fastq.gz"
OUT1="$PATHS_HOSTFREE/${SRR}_1.fastq.gz"; OUT2="$PATHS_HOSTFREE/${SRR}_2.fastq.gz"
require "$IN1" "$IN2" "$REFERENCE_HUMAN_BOWTIE2_INDEX.1.bt2"
mkdir -p "$PATHS_HOSTFREE" "$PATHS_HOST/bowtie2_logs"
LOG="$PATHS_HOST/bowtie2_logs/${SRR}.bowtie2.log"
echo "[$(date)] $SRR: stage 3 host depletion"

# shellcheck disable=SC2086
bowtie2 -p "$T" $PARAMS_HOST_BOWTIE2_ARGS -x "$REFERENCE_HUMAN_BOWTIE2_INDEX" \
        -1 "$IN1" -2 "$IN2" 2> "$LOG" \
    | samtools view -u -f 12 -F 256 - \
    | samtools fastq -@ 4 -n -c 6 -1 "$OUT1.tmp.gz" -2 "$OUT2.tmp.gz" -0 /dev/null -s /dev/null -
mv "$OUT1.tmp.gz" "$OUT1"; mv "$OUT2.tmp.gz" "$OUT2"
cat "$LOG"

TRIMMED=$(awk '/reads; of these:/{print $1; exit}' "$LOG")
R1=$(zcat "$OUT1" | awk 'END{print NR/4}')
R2=$(zcat "$OUT2" | awk 'END{print NR/4}')
[[ "$R1" == "$R2" ]] || { echo "R1/R2 pair counts differ: $R1 vs $R2" >&2; exit 1; }
printf 'run\ttrimmed_pairs\thuman_pairs\tretained_pairs\n%s\t%s\t%s\t%s\n' \
    "$SRR" "$TRIMMED" "$((TRIMMED - R1))" "$R1" > "$PATHS_HOST/${SRR}.counts.tsv"
cat "$PATHS_HOST/${SRR}.counts.tsv"
echo "[$(date)] $SRR: done"
