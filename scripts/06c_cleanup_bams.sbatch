#!/bin/bash
#SBATCH --job-name=pd_bam_cleanup
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 06c_cleanup_bams.sbatch -- stage 6 (all runs together): delete the cross-sample BAMs from
# 04c (another group mate's reads on this run's assembly) once everything that reads them
# is done: 05d (coverm_contig.tsv, amr_counts.tsv) and 06a (depth.txt). Their only later
# use, the stage 8 coverage correlation, reads coverm_contig.tsv. Before deleting, checks
# that coverm_contig.tsv and depth.txt have one column set per group mate, i.e. that they
# were built from this mapping set. Kept: the run's own BAM (+ .bai) and every
# .bowtie2.log. Deleted files are listed in <asm>/map/DELETED_BAMS.txt.
# Re-running 05d or 06a afterwards needs 04c again.
# Submit:  sbatch scripts/06c_cleanup_bams.sbatch      (not an array; after 05d and 06a)

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"

# Exact header match (CoverM: "<run> Mean"; jgi depth: "<run>", "<run>-var")
has_column(){ awk -F'\t' -v c="$2" 'NR==1{for(i=1;i<=NF;i++) if($i==c) f=1; exit !f}' "$1"; }

# Check every run first, delete only if all pass
while read -r ASM; do
    require_group "$ASM" > /dev/null
    mapfile -t MATES < <(group_mates "$ASM")
    COV="$PATHS_CONTIGS/$ASM/coverm_contig.tsv"; DEPTH="$PATHS_BINNING/$ASM/depth.txt"
    require "$COV" "$PATHS_CONTIGS/$ASM/amr_counts.tsv" "$DEPTH"
    for R in "${MATES[@]}"; do
        has_column "$COV" "$R Mean" || { echo "$COV has no column for $R: re-run 05d" >&2; exit 1; }
        has_column "$DEPTH" "$R" || { echo "$DEPTH has no column for $R: re-run 06a" >&2; exit 1; }
    done
    echo "$ASM: 05d + 06a cover all ${#MATES[@]} group mates"
done < <(run_samples)

while read -r ASM; do
    MAP="$PATHS_ASSEMBLY/$ASM/map"
    n=0
    for BAM in "$MAP"/*.bam; do
        [[ -e "$BAM" ]] || continue
        [[ "$(basename "$BAM" .bam)" == "$ASM" ]] && continue
        echo "$(date +%F)	$(basename "$BAM")	$(stat -c %s "$BAM")" >> "$MAP/DELETED_BAMS.txt"
        rm -f "$BAM" "$BAM.bai"
        n=$(( n + 1 ))
    done
    echo "$ASM: deleted $n cross-sample BAMs; kept $(ls "$MAP"/*.bam | xargs -n1 basename)"
done < <(run_samples)
echo "[$(date)] done"
