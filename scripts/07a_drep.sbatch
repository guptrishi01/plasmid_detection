#!/bin/bash
#SBATCH --job-name=pd_drep
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 07a_drep.sbatch -- stage 7 (all runs together), dRep at params.mags.drep_ani (95% ANI)
# in two levels, using CheckM2's completeness and contamination (--genomeInfo):
#   1. per binning group (params.binning.groups): the medium-and-up MAGs in
#      $paths.binning/pass/ from the group's runs -> one MAG per species per group
#      ($paths.mags/groups/<group>/drep/)
#   2. catalog: every group's dereplicated MAGs together -> one MAG per species overall
#      ($paths.mags/drep/), which GTDB-Tk (7b), CoverM over every run (7c) and stage 8 use
# dRep's own quality filter is set to the MIMAG medium thresholds so it doesn't silently
# drop MAGs stage 6 passed (its defaults are 75/25). A level with one genome skips dRep and
# gets the same layout by hand.
# Output: $paths.mags/groups/<group>/drep/, $paths.mags/drep/ (dereplicated_genomes/,
#         data_tables/Cdb.csv, Wdb.csv)
# Submit:  sbatch scripts/07a_drep.sbatch      (not an array)

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
# Medium+ MAGs of the runs in paths.run_samples (pass/ can hold runs from an earlier run list)
PASS=()
while read -r run; do
    for f in "$PATHS_BINNING/pass/${run}".bin.*.fa; do [[ -e "$f" ]] && PASS+=("$f"); done
done < <(run_samples)
(( ${#PASS[@]} )) || { echo "no medium+ MAGs in $PATHS_BINNING/pass (run 06b)" >&2; exit 1; }
mkdir -p "$PATHS_MAGS"
pdpy plasmid_detection.tables drep-info --binning-dir "$PATHS_BINNING" --samples "$PATHS_RUN_SAMPLES" \
    --out "$PATHS_MAGS/genome_info.csv"

# dereplicate <out_dir> <genome.fa>...
dereplicate(){
    local out=$1; shift
    rm -rf "$out"
    if (( $# == 1 )); then
        # One MAG: nothing to dereplicate; lay out what dRep would have produced
        mkdir -p "$out/dereplicated_genomes" "$out/data_tables"
        cp "$1" "$out/dereplicated_genomes/"
        local g; g=$(basename "$1")
        printf 'genome,secondary_cluster\n%s,1_1\n' "$g" > "$out/data_tables/Cdb.csv"
        printf 'genome,cluster\n%s,1_1\n' "$g" > "$out/data_tables/Wdb.csv"
        return 0
    fi
    activate "$ENVS_DREP"
    dRep dereplicate "$out" -g "$@" -p "$PARAMS_MAGS_THREADS" \
        -sa "$PARAMS_MAGS_DREP_ANI" --genomeInfo "$PATHS_MAGS/genome_info.csv" \
        -comp "$PARAMS_BINNING_MIMAG_MEDIUM_COMPLETENESS" -con "$PARAMS_BINNING_MIMAG_MEDIUM_CONTAMINATION"
    deactivate
}

# Every run must be in a group, or its MAGs would silently miss both levels
while read -r run; do require_group "$run" > /dev/null; done < <(run_samples)

# 1. per group
rm -rf "$PATHS_MAGS/groups"
echo "[$(date)] level 1: dRep per binning group (${#PASS[@]} MAGs)"
for spec in $PARAMS_BINNING_GROUPS; do
    g=${spec%%:*}
    GENOMES=()
    while read -r run; do
        [[ "$(run_group "$run")" == "$g" ]] || continue
        for f in "$PATHS_BINNING/pass/${run}".bin.*.fa; do [[ -e "$f" ]] && GENOMES+=("$f"); done
    done < <(run_samples)
    (( ${#GENOMES[@]} )) || continue
    dereplicate "$PATHS_MAGS/groups/$g/drep" "${GENOMES[@]}"
    echo "  $g: ${#GENOMES[@]} MAGs -> $(ls "$PATHS_MAGS/groups/$g/drep/dereplicated_genomes" | wc -l)"
done

# 2. catalog across groups
GROUP_REPS=("$PATHS_MAGS"/groups/*/drep/dereplicated_genomes/*.fa)
echo "[$(date)] level 2: catalog dRep over ${#GROUP_REPS[@]} group MAGs"
dereplicate "$PATHS_MAGS/drep" "${GROUP_REPS[@]}"
echo "catalog MAGs: $(ls "$PATHS_MAGS/drep/dereplicated_genomes" | wc -l)"
echo "[$(date)] done"
