#!/bin/bash
#SBATCH --job-name=pd_scapp
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 04e_scapp.sbatch -- stage 4, one run per array task: SCAPP on the MEGAHIT assembly graph
# (k<K>.fastg from 04a, with -k K) plus the run's host-depleted reads, for circular
# plasmids. Output: $paths.assembly/<run>/scapp/
# Submit:  sbatch --array=1-$(grep -c . metadata/pilot_samples.txt) scripts/04e_scapp.sbatch
# Tune:    params.assembly.scapp_args (e.g. "-l 1000 -m 0.5")

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
activate "$ENVS_SCAPP"
SRR=$(task_sample)
DIR="$PATHS_ASSEMBLY/$SRR"
require "$DIR/graph_k.txt"
K=$(cat "$DIR/graph_k.txt")
GRAPH="$DIR/${SRR}.k${K}.fastg"
IN1="$PATHS_HOSTFREE/${SRR}_1.fastq.gz"; IN2="$PATHS_HOSTFREE/${SRR}_2.fastq.gz"
require "$GRAPH" "$IN1" "$IN2"
rm -rf "$DIR/scapp"
echo "[$(date)] $SRR: SCAPP on k=$K graph"

# shellcheck disable=SC2086
if ! scapp -g "$GRAPH" -k "$K" -r1 "$IN1" -r2 "$IN2" -o "$DIR/scapp" \
        -p "$PARAMS_ASSEMBLY_SCAPP_THREADS" $PARAMS_ASSEMBLY_SCAPP_ARGS; then
    # SCAPP 0.1.4 crashes in makeblastdb when the graph yields no candidate cycles (empty
    # *.cycs.fasta). That is a result (no circular plasmids), not a failure; anything else is.
    cycs=$(ls "$DIR"/scapp/intermediate_files/*.cycs.fasta 2>/dev/null | head -1)
    if [[ -n "$cycs" && ! -s "$cycs" ]]; then
        echo "SCAPP: no candidate cycles in the graph -> 0 circular plasmids" | tee "$DIR/scapp/NO_CYCLES"
    else
        echo "SCAPP failed (see log above)" >&2; exit 1
    fi
fi
ls -la "$DIR/scapp"
echo "[$(date)] $SRR: done"
