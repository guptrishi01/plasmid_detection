#!/bin/bash
#SBATCH --job-name=pd_cov_tag
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 05d_coverage_tag.sbatch -- stage 5, one run per array task:
#   CoverM contig over every stage-4 BAM on this assembly (coverage + TPM per contig per
#   reads-run), then tag each AMR hit (AMRFinderPlus and RGI) as plasmid, virus, or
#   chromosome by the contig it sits on (geNomad calls; a hit inside a geNomad provirus
#   region is "virus"; "chromosome" = neither plasmid nor virus), then count reads over
#   each hit's coordinates in every BAM (samtools bedcov -c).
# Output: $paths.contigs/<run>/{coverm_contig.tsv,amr_tagged.tsv,amr_counts.tsv}
# Submit:  sbatch --array=1-$(grep -c . metadata/pilot_samples.txt) scripts/05d_coverage_tag.sbatch
# Tune:    params.contigs.coverm_methods

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
SRR=$(task_sample)
DIR="$PATHS_CONTIGS/$SRR"
BAMS=("$PATHS_ASSEMBLY/$SRR"/map/*.bam)
require "${BAMS[@]}" "$DIR/amrfinder.tsv" "$DIR/rgi/$SRR.rgi.txt"
echo "[$(date)] $SRR: CoverM contig over ${#BAMS[@]} BAMs + AMR tagging"

activate "$ENVS_COVERM"
# shellcheck disable=SC2086
coverm contig --bam-files "${BAMS[@]}" --methods $PARAMS_CONTIGS_COVERM_METHODS \
    --threads "$PARAMS_CONTIGS_THREADS" --output-file "$DIR/coverm_contig.tsv"
deactivate

pdpy plasmid_detection.tables tag-amr --run "$SRR" \
    --genomad-dir "$DIR/genomad" --amrfinder "$DIR/amrfinder.tsv" --rgi "$DIR/rgi/$SRR.rgi.txt" \
    --out "$DIR/amr_tagged.tsv" --bed "$DIR/amr_hits.bed"
echo "hits by source x tag:"; cut -f2,3 "$DIR/amr_tagged.tsv" | tail -n +2 | sort | uniq -c

# Reads over each AMR hit's own coordinates, per BAM (integer counts for the stage 9
# resistome matrix and statistics)
activate "$ENVS_ASSEMBLY"
samtools bedcov -c "$DIR/amr_hits.bed" "${BAMS[@]}" > "$DIR/amr_bedcov.txt"
deactivate
pdpy plasmid_detection.tables amr-counts --bedcov "$DIR/amr_bedcov.txt" --bams "${BAMS[@]}" \
    --out "$DIR/amr_counts.tsv"
echo "[$(date)] $SRR: done"
