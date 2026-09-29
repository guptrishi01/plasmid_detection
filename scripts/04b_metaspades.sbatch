#!/bin/bash
#SBATCH --job-name=pd_metaspades
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=3-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 04b_metaspades.sbatch -- stage 4: one metaSPAdes attempt on host-depleted reads of
# params.assembly.metaspades_run, or (if empty) the smallest run in paths.run_samples. "Record what happens": the
# outcome is written to $paths.assembly/metaspades/attempt.tsv whether it succeeds or not
# (exit code, wall time, peak memory, contigs if any). A failure is a result, not an error.
# Submit:  sbatch scripts/04b_metaspades.sbatch      (not an array)
# Tune:    params.assembly.metaspades_memory_gb (keep below --mem above)

set -uo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_ASSEMBLY"
T=$PARAMS_ASSEMBLY_THREADS

SRR=$PARAMS_ASSEMBLY_METASPADES_RUN
if [[ -z "$SRR" ]]; then
    SRR=$(while read -r s; do printf '%s\t%s\n' "$s" "$(sample_field "$s" spots)"; done < <(run_samples) \
          | sort -t$'\t' -k2,2n | head -1 | cut -f1)
fi
IN1="$PATHS_HOSTFREE/${SRR}_1.fastq.gz"; IN2="$PATHS_HOSTFREE/${SRR}_2.fastq.gz"
require "$IN1" "$IN2"
OUT="$PATHS_ASSEMBLY/metaspades"
rm -rf "$OUT/$SRR"; mkdir -p "$OUT"
echo "[$(date)] metaSPAdes on $SRR ($(sample_field "$SRR" spots) spots)"

start=$(date +%s)
/usr/bin/time -v -o "$OUT/${SRR}.time.txt" \
    spades.py --meta -1 "$IN1" -2 "$IN2" -o "$OUT/$SRR" -t "$T" -m "$PARAMS_ASSEMBLY_METASPADES_MEMORY_GB"
rc=$?
elapsed=$(( $(date +%s) - start ))
maxrss_gb=$(awk -F': ' '/Maximum resident set size/{printf "%.1f", $2/1024/1024}' "$OUT/${SRR}.time.txt")
contigs=NA
[[ -s "$OUT/$SRR/contigs.fasta" ]] && contigs=$(grep -c '^>' "$OUT/$SRR/contigs.fasta")
last_error=$(grep -iE "error|exception|killed|not enough memory" "$OUT/$SRR/spades.log" 2>/dev/null | tail -1)

printf 'run\tspots\texit_code\telapsed_h\tmax_rss_gb\tcontigs\tlast_error_line\n' > "$OUT/attempt.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$SRR" "$(sample_field "$SRR" spots)" "$rc" \
    "$(awk -v e="$elapsed" 'BEGIN{printf "%.2f", e/3600}')" "$maxrss_gb" "$contigs" "${last_error:-none}" \
    >> "$OUT/attempt.tsv"
cat "$OUT/attempt.tsv"
echo "[$(date)] metaSPAdes attempt recorded (exit $rc)"
exit 0
