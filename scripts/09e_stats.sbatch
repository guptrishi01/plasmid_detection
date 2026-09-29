#!/bin/bash
#SBATCH --job-name=pd_stats
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 09e_stats.sbatch -- stage 9 statistics on the AMRFinderPlus resistome (all and mobile):
#   sample metadata + BioSample antibiotic concentrations (covariate; ND/<LOQ encoded per
#   params.stats.nd_encoding / loq_encoding, raw values kept) -> metadata/stats_samples.tsv
#   CLR; ALDEx2 or ANCOM-BC2 (params.stats.method) for each params.stats.comparisons pair
#   within each system and within each time point; PERMANOVA (Aitchison) as exploratory;
#   Benjamini-Hochberg. Controls are a reference only: never subtracted, never tested.
#   Tests with fewer than params.stats.min_group_size samples per group are skipped and
#   logged (the 4-run pilot skips all of them; that is expected).
# Output: $paths.rollup/stats/

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
# shellcheck disable=SC2086
pdpy plasmid_detection.tables stats-samples --samples "$PATHS_RUN_SAMPLES" --sheet "$PATHS_SAMPLES" \
    --site-groups $PARAMS_STATS_SITE_GROUPS --nd-encoding "$PARAMS_STATS_ND_ENCODING" \
    --loq-encoding "$PARAMS_STATS_LOQ_ENCODING" --out "$PATHS_METADATA/stats_samples.tsv"

activate "$ENVS_STATS"
mkdir -p "$PATHS_ROLLUP/stats"
for part in all mobile; do
    Rscript "$PD_REPO/scripts/R/stats.R" "$PATHS_METADATA/stats_samples.tsv" \
        "$PATHS_ROLLUP/resistome_amrfinder_${part}.tsv" "$PATHS_ROLLUP/stats/$part"
done
