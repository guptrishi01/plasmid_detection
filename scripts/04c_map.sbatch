#!/bin/bash
#SBATCH --job-name=pd_map
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 04c_map.sbatch -- stage 4, one ASSEMBLY per array task: Bowtie2 maps the run's own reads
# back to its assembly, and the reads of every other run of the same system (location in
# the sample sheet, restricted to paths.run_samples) to it too. One sorted, indexed BAM
# per (assembly, reads) pair gives each contig a coverage profile across samples, which
# binning (stage 6) and plasmid-host coverage correlation (stage 8) use.
# Output: $paths.assembly/<asm>/map/<reads_run>.bam (+ .bai, .bowtie2.log)
# Submit:  sbatch --array=1-$(grep -c . metadata/pilot_samples.txt) scripts/04c_map.sbatch
# Tune:    params.assembly.bowtie2_args, params.assembly.map_threads

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_ASSEMBLY"
ASM=$(task_sample)
T=$PARAMS_ASSEMBLY_MAP_THREADS
DIR="$PATHS_ASSEMBLY/$ASM"
CONTIGS="$DIR/${ASM}.contigs.fa"
require "$CONTIGS"
mkdir -p "$DIR/map" "$DIR/bt2_index"

mapfile -t MATES < <(system_mates "$ASM")
echo "[$(date)] $ASM ($(sample_field "$ASM" location)): mapping ${#MATES[@]} runs: ${MATES[*]}"

bowtie2-build --threads "$T" "$CONTIGS" "$DIR/bt2_index/$ASM" > "$DIR/bt2_index/build.log" 2>&1

for R in "${MATES[@]}"; do
    IN1="$PATHS_HOSTFREE/${R}_1.fastq.gz"; IN2="$PATHS_HOSTFREE/${R}_2.fastq.gz"
    require "$IN1" "$IN2"
    BAM="$DIR/map/${R}.bam"
    # shellcheck disable=SC2086
    bowtie2 -p "$T" $PARAMS_ASSEMBLY_BOWTIE2_ARGS -x "$DIR/bt2_index/$ASM" -1 "$IN1" -2 "$IN2" \
            2> "$DIR/map/${R}.bowtie2.log" \
        | samtools sort -@ 4 -m 2G -o "$BAM.tmp" -
    mv "$BAM.tmp" "$BAM"
    samtools index "$BAM"
    echo "  $R -> $ASM: $(tail -1 "$DIR/map/${R}.bowtie2.log")"
done
rm -rf "$DIR/bt2_index"
echo "[$(date)] $ASM: done"
