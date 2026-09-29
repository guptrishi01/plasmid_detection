#!/bin/bash
#SBATCH --job-name=pd_plasmid_hosts
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 08b_plasmid_hosts.sbatch -- stage 8: which geNomad plasmid contigs sit in which MAG, and
# three plasmid-to-host inferences, reported side by side (all three are inferences):
#   1. co-binning: the plasmid contig is in a medium+ MetaBAT2 bin (-> its dRep MAG)
#   2. MOB-suite host range (mob_typer predicted_host_range_overall_rank/name)
#   3. coverage correlation: Spearman between the plasmid contig's coverage profile and
#      each medium+ bin's mean profile across the stage-4 BAMs on that assembly. Only
#      computed with >= params.mag_links.min_samples_for_correlation BAMs (else blank).
# Output: $paths.mag_links/{plasmids_in_mags.tsv,plasmid_host_inference.tsv}, copied to tables/
# Submit:  sbatch scripts/08b_plasmid_hosts.sbatch      (not an array)

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
mkdir -p "$PATHS_MAG_LINKS" "$PATHS_TABLES"
pdpy plasmid_detection.tables plasmid-hosts --samples "$PATHS_RUN_SAMPLES" \
    --binning-dir "$PATHS_BINNING" --mags-dir "$PATHS_MAGS" --contigs-dir "$PATHS_CONTIGS" \
    --min-samples "$PARAMS_MAG_LINKS_MIN_SAMPLES_FOR_CORRELATION" --out-dir "$PATHS_MAG_LINKS"
cp "$PATHS_MAG_LINKS/plasmids_in_mags.tsv" "$PATHS_TABLES/08_plasmids_in_mags.tsv"
cp "$PATHS_MAG_LINKS/plasmid_host_inference.tsv" "$PATHS_TABLES/08_plasmid_host_inference.tsv"
echo "plasmid contigs in medium+ MAGs: $(( $(wc -l < "$PATHS_TABLES/08_plasmids_in_mags.tsv") - 1 ))"
echo "plasmid contigs with any inference row: $(( $(wc -l < "$PATHS_TABLES/08_plasmid_host_inference.tsv") - 1 ))"
