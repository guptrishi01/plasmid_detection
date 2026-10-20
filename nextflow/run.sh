#!/bin/bash
# run.sh -- launch nextflow/main.nf with config/config.yaml as the params-file.
#   bash nextflow/run.sh -profile slurm,pilot -resume     # 4-run pilot
#   bash nextflow/run.sh -profile slurm -resume           # all 92 WGS runs of PRJNA415974
#   bash nextflow/run.sh -profile slurm --id SRP121672 --runs SRR6231181,SRR6231191
# Extra arguments go straight to `nextflow run`. On the cluster, submit nextflow/submit.sbatch.
set -euo pipefail
PD_REPO=${PD_REPO:-$HOME/gibas/plasmid_detection}
export NXF_HOME=${NXF_HOME:-/scratch/rgupta25/plasmid_detection/nxf_home}
export NXF_CONDA_CACHEDIR=${NXF_CONDA_CACHEDIR:-$NXF_HOME/conda}
mkdir -p "$NXF_HOME"
PARAMS="$NXF_HOME/params.$$.json"
trap 'rm -f "$PARAMS"' EXIT
PYTHONPATH="$PD_REPO/src" "$PD_REPO/conda/envs/pd-core/bin/python" -m plasmid_detection.config nextflow > "$PARAMS"
cd "$PD_REPO"
"$PD_REPO/conda/envs/pd-nextflow/bin/nextflow" run "$PD_REPO/nextflow/main.nf" -params-file "$PARAMS" "$@"
