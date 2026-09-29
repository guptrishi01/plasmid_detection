#!/bin/bash
#SBATCH --job-name=pd_mag_abund
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 07c_mag_abundance.sbatch -- stage 7: CoverM genome over ALL runs in paths.run_samples
# (host-depleted pairs) against the dereplicated MAGs, then the per-MAG table: quality
# (CheckM2 + MIMAG tier), taxonomy (GTDB-Tk), abundance per sample.
# CoverM genome does its own read mapping (default minimap2-sr); Bowtie2 is the stage-4
# mapper only.
# Output: $paths.mags/coverm_genome.tsv, $paths.tables/07_mags.tsv
# Submit:  sbatch scripts/07c_mag_abundance.sbatch     (not an array)
# Tune:    params.mags.coverm_genome_methods

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_COVERM"
GENOMES="$PATHS_MAGS/drep/dereplicated_genomes"
[[ -n "$(ls "$GENOMES"/*.fa 2>/dev/null)" ]] || { echo "no genomes in $GENOMES (run 07a)" >&2; exit 1; }
COUPLED=()
while read -r s; do
    require "$PATHS_HOSTFREE/${s}_1.fastq.gz" "$PATHS_HOSTFREE/${s}_2.fastq.gz"
    COUPLED+=("$PATHS_HOSTFREE/${s}_1.fastq.gz" "$PATHS_HOSTFREE/${s}_2.fastq.gz")
done < <(run_samples)
echo "[$(date)] CoverM genome: $(( ${#COUPLED[@]} / 2 )) runs vs $(ls "$GENOMES"/*.fa | wc -l) MAGs"

# shellcheck disable=SC2086
coverm genome --coupled "${COUPLED[@]}" --genome-fasta-directory "$GENOMES" -x fa \
    --methods $PARAMS_MAGS_COVERM_GENOME_METHODS --threads "$PARAMS_MAGS_THREADS" \
    --output-file "$PATHS_MAGS/coverm_genome.tsv"
deactivate

pdpy plasmid_detection.tables mags --samples "$PATHS_RUN_SAMPLES" --binning-dir "$PATHS_BINNING" \
    --mags-dir "$PATHS_MAGS" --out "$PATHS_TABLES/07_mags.tsv"
column -t "$PATHS_TABLES/07_mags.tsv" | head -30
echo "[$(date)] done"
