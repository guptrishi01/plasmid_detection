#!/bin/bash
#SBATCH --job-name=pd_checkv_mob
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 05b_checkv_mobtyper.sbatch -- stage 5, one run per array task: CheckV on geNomad's virus
# contigs; MOB-suite mob_typer (--multi: each contig typed as its own plasmid) on geNomad's
# plasmid contigs. An empty contig set is recorded, not an error.
# Output: $paths.contigs/<run>/checkv/, $paths.contigs/<run>/mob_typer.tsv
# Submit:  sbatch --array=1-$(grep -c . metadata/pilot_samples.txt) scripts/05b_checkv_mobtyper.sbatch

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
SRR=$(task_sample)
T=$PARAMS_CONTIGS_THREADS
DIR="$PATHS_CONTIGS/$SRR"
VIRUS=$(ls "$DIR"/genomad/*_summary/*_virus.fna 2>/dev/null || true)
PLASMID=$(ls "$DIR"/genomad/*_summary/*_plasmid.fna 2>/dev/null || true)
[[ -n "$VIRUS" && -n "$PLASMID" ]] || { echo "no geNomad summary in $DIR/genomad; run 05a first" >&2; exit 1; }
echo "[$(date)] $SRR: CheckV + mob_typer"

rm -rf "$DIR/checkv"
if [[ -s "$VIRUS" ]]; then
    activate "$ENVS_CHECKV"
    checkv end_to_end "$VIRUS" "$DIR/checkv" -d "$DB_CHECKV" -t "$T"
    deactivate
else
    mkdir -p "$DIR/checkv"; echo "no virus contigs" > "$DIR/checkv/EMPTY"
fi

rm -f "$DIR/mob_typer.tsv"
if [[ -s "$PLASMID" ]]; then
    activate "$ENVS_MOBSUITE"
    mob_typer --multi --infile "$PLASMID" --out_file "$DIR/mob_typer.tsv" \
        --database_directory "$DB_MOB_SUITE" --num_threads "$T"
    deactivate
else
    echo "no plasmid contigs" > "$DIR/mob_typer.EMPTY"
fi
echo "[$(date)] $SRR: done"
