#!/bin/bash
#SBATCH --job-name=pd_drep
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 07a_drep.sbatch -- stage 7 (all runs together): dRep at params.mags.drep_ani (95% ANI)
# over every medium-and-up MAG in $paths.binning/pass/, using CheckM2's completeness and
# contamination (--genomeInfo). dRep's own quality filter is set to the MIMAG medium
# thresholds so it doesn't silently drop MAGs stage 6 passed (its defaults are 75/25).
# Output: $paths.mags/drep/ (dereplicated_genomes/, data_tables/Cdb.csv, Wdb.csv)
# Submit:  sbatch scripts/07a_drep.sbatch      (not an array)

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
PASS=("$PATHS_BINNING"/pass/*.fa)
[[ -e "${PASS[0]}" ]] || { echo "no medium+ MAGs in $PATHS_BINNING/pass (run 06b)" >&2; exit 1; }
mkdir -p "$PATHS_MAGS"
echo "[$(date)] dRep over ${#PASS[@]} MAGs"

pdpy plasmid_detection.tables drep-info --binning-dir "$PATHS_BINNING" --samples "$PATHS_RUN_SAMPLES" \
    --out "$PATHS_MAGS/genome_info.csv"

rm -rf "$PATHS_MAGS/drep"
if (( ${#PASS[@]} == 1 )); then
    # One MAG: nothing to dereplicate; lay out what dRep would have produced
    mkdir -p "$PATHS_MAGS/drep/dereplicated_genomes" "$PATHS_MAGS/drep/data_tables"
    cp "${PASS[0]}" "$PATHS_MAGS/drep/dereplicated_genomes/"
    g=$(basename "${PASS[0]}")
    printf 'genome,secondary_cluster\n%s,1_1\n' "$g" > "$PATHS_MAGS/drep/data_tables/Cdb.csv"
    printf 'genome,cluster\n%s,1_1\n' "$g" > "$PATHS_MAGS/drep/data_tables/Wdb.csv"
    echo "only one MAG passed; dRep skipped"
    exit 0
fi
activate "$ENVS_DREP"
dRep dereplicate "$PATHS_MAGS/drep" -g "${PASS[@]}" -p "$PARAMS_MAGS_THREADS" \
    -sa "$PARAMS_MAGS_DREP_ANI" --genomeInfo "$PATHS_MAGS/genome_info.csv" \
    -comp "$PARAMS_BINNING_MIMAG_MEDIUM_COMPLETENESS" -con "$PARAMS_BINNING_MIMAG_MEDIUM_CONTAMINATION"
echo "dereplicated genomes: $(ls "$PATHS_MAGS/drep/dereplicated_genomes" | wc -l)"
echo "[$(date)] done"
