#!/bin/bash
#SBATCH --job-name=pd_plasmid_callers
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 05f_plasmid_callers.sbatch -- stage 5 (all runs): compare plasmid call sets per run --
# geNomad, PLASMe (main run, plasme/), PLASMe high-precision (plasme_high-precision/, from
# `05e_plasme.sbatch high-precision`), each PLASMe set also with calls that are mostly
# chromosome-shared removed (params.contigs.plasme_max_ambiguous_fraction). Reports size,
# overlap with geNomad, contigs with AMR-type genes, and geNomad's plasmid score of the calls.
# Output: $paths.tables/05_plasmid_callers.tsv
# Submit:  sbatch scripts/05f_plasmid_callers.sbatch     (after 05a, 05c, 05e)

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
mkdir -p "$PATHS_TABLES"
pdpy plasmid_detection.tables plasmid-callers --samples "$PATHS_RUN_SAMPLES" \
    --contigs-dir "$PATHS_CONTIGS" --max-amb "$PARAMS_CONTIGS_PLASME_MAX_AMBIGUOUS_FRACTION" \
    --out "$PATHS_TABLES/05_plasmid_callers.tsv"
column -t "$PATHS_TABLES/05_plasmid_callers.tsv"
