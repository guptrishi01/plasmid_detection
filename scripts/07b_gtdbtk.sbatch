#!/bin/bash
#SBATCH --job-name=pd_gtdbtk
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 07b_gtdbtk.sbatch -- stage 7: GTDB-Tk classify_wf (GTDB r226, via GTDBTK_DATA_PATH set in
# pd-gtdbtk) on the dereplicated MAGs.
# Output: $paths.mags/gtdbtk/ (gtdbtk.bac120.summary.tsv, gtdbtk.ar53.summary.tsv)
# Submit:  sbatch scripts/07b_gtdbtk.sbatch     (not an array)
# Tune:    params.mags.gtdbtk_args

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_GTDBTK"
GENOMES="$PATHS_MAGS/drep/dereplicated_genomes"
[[ -n "$(ls "$GENOMES"/*.fa 2>/dev/null)" ]] || { echo "no genomes in $GENOMES (run 07a)" >&2; exit 1; }
echo "[$(date)] GTDB-Tk on $(ls "$GENOMES"/*.fa | wc -l) MAGs; GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH"

rm -rf "$PATHS_MAGS/gtdbtk"
# shellcheck disable=SC2086
gtdbtk classify_wf --genome_dir "$GENOMES" -x fa --out_dir "$PATHS_MAGS/gtdbtk" \
    --cpus "$PARAMS_MAGS_THREADS" $PARAMS_MAGS_GTDBTK_ARGS
echo "[$(date)] done"
