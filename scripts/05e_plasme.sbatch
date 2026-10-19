#!/bin/bash
#SBATCH --job-name=pd_plasme
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 05e_plasme.sbatch -- stage 5, one run per array task: PLASMe on the whole assembly, a
# second plasmid caller next to geNomad (in place of "PlasLLM", which doesn't exist).
# BLASTn against PLSDB for close relatives, order-specific Transformer models for diverged
# plasmids. Runs on CPU (pd-plasme has CPU-only torch). Independent of 05a-05d.
# Output: $paths.contigs/<run>/plasme/<run>.plasme.fna (+ _report.csv)
# Submit:  sbatch --array=1-$(grep -c . metadata/pilot_samples.txt) scripts/05e_plasme.sbatch [MODE]
# Tune:    params.contigs.plasme_args (-m balance | high-sensitivity | high-precision)
# MODE (optional, for caller comparisons): runs `-m MODE` instead of plasme_args and writes
# to plasme_MODE/, so it never overwrites the main run in plasme/.

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_PLASME"
SRR=$(task_sample)
CONTIGS="$PATHS_ASSEMBLY/$SRR/${SRR}.contigs.fa"
MODE=${1:-}
if [[ -n "$MODE" ]]; then ARGS="-m $MODE"; OUT="$PATHS_CONTIGS/$SRR/plasme_$MODE"
else ARGS=$PARAMS_CONTIGS_PLASME_ARGS; OUT="$PATHS_CONTIGS/$SRR/plasme"; fi
require "$CONTIGS" "$DB_PLASME_DIR/plsdb_Mar30.dmnd"   # .done is empty; require needs -s
# PLASMe exits 0 without running if the output exists: start clean
rm -rf "$OUT"; mkdir -p "$OUT"
TMP=$(mktemp -d "$PATHS_WORK/plasme.XXXX"); trap 'rm -rf "$TMP"' EXIT
echo "[$(date)] $SRR: PLASMe ($ARGS) -> $OUT"

# shellcheck disable=SC2086
python "$DB_PLASME_CODE/PLASMe.py" "$CONTIGS" "$OUT/${SRR}.plasme.fna" \
    -d "$DB_PLASME_DIR" -t "$PARAMS_CONTIGS_THREADS" --temp "$TMP" $ARGS
# PLASMe runs BLAST/DIAMOND/Prodigal without checking their exit codes: the report is
# the proof it finished
require "$OUT/${SRR}.plasme.fna_report.csv"

n=$(grep -c '^>' "$OUT/${SRR}.plasme.fna" || true)
echo "PLASMe plasmid contigs: $n"
G=$(ls "$PATHS_CONTIGS/$SRR"/genomad/*_summary/*_plasmid_summary.tsv 2>/dev/null || true)
if [[ -n "$G" ]]; then
    both=$(comm -12 <(grep '^>' "$OUT/${SRR}.plasme.fna" | cut -c2- | cut -d' ' -f1 | sort -u) \
                    <(tail -n +2 "$G" | cut -f1 | sort -u) | wc -l)
    echo "also called plasmid by geNomad: $both of $n (geNomad: $(( $(wc -l < "$G") - 1 )))"
fi
echo "[$(date)] $SRR: done"
