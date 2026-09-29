#!/bin/bash
#SBATCH --job-name=pd_mag_amr
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 08a_mag_amr.sbatch -- stage 8: AMRFinderPlus (params.contigs.amrfinder_args, i.e. --plus)
# on each dereplicated MAG; per-MAG outputs combined into one table.
# Output: $paths.mag_links/amrfinder/<mag>.tsv, $paths.tables/08_mag_amr.tsv
# Submit:  sbatch scripts/08a_mag_amr.sbatch      (not an array)

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_AMRFINDERPLUS"
GENOMES=("$PATHS_MAGS"/drep/dereplicated_genomes/*.fa)
[[ -e "${GENOMES[0]}" ]] || { echo "no dereplicated MAGs (run 07a)" >&2; exit 1; }
mkdir -p "$PATHS_MAG_LINKS/amrfinder" "$PATHS_TABLES"
echo "[$(date)] AMRFinderPlus on ${#GENOMES[@]} MAGs"

for g in "${GENOMES[@]}"; do
    mag=$(basename "$g" .fa)
    # shellcheck disable=SC2086
    amrfinder --nucleotide "$g" $PARAMS_CONTIGS_AMRFINDER_ARGS --database "$DB_AMRFINDERPLUS" \
        --threads "$PARAMS_MAG_LINKS_THREADS" --name "$mag" --output "$PATHS_MAG_LINKS/amrfinder/$mag.tsv"
done
first=1
for f in "$PATHS_MAG_LINKS"/amrfinder/*.tsv; do
    if (( first )); then cat "$f"; first=0; else tail -n +2 "$f"; fi
done > "$PATHS_TABLES/08_mag_amr.tsv"
echo "hits: $(( $(wc -l < "$PATHS_TABLES/08_mag_amr.tsv") - 1 ))"
echo "[$(date)] done"
