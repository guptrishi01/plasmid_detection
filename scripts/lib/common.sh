# Sourced by every scripts/*.sbatch: loads conda and exports config/config.yaml as env vars
# (paths.raw_reads -> PATHS_RAW_READS, ...). Batch scripts run from a Slurm spool copy, so
# the repo location comes from PD_REPO, not from this file's path.

PD_REPO=${PD_REPO:-$HOME/gibas/plasmid_detection}
# The core env is found before the config is read, so its location is fixed here
PD_CORE_ENV=${PD_CORE_ENV:-$PD_REPO/conda/envs/pd-core}

module load mambaforge/23.11
source "$(conda info --base)/etc/profile.d/conda.sh"

if [[ ! -x "$PD_CORE_ENV/bin/python" ]]; then
    echo "pd-core env not found at $PD_CORE_ENV; run scripts/00_setup.sbatch first" >&2
    exit 1
fi
eval "$(PYTHONPATH="$PD_REPO/src" "$PD_CORE_ENV/bin/python" -m plasmid_detection.config shell)"

export CONDA_PKGS_DIRS=$PATHS_CONDA_PKGS
