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

# ---- helpers for stages 2-9 ------------------------------------------------------------

# Runs processed by stages 2-9 (paths.run_samples), one per line
run_samples(){ grep -v '^\s*$' "$PATHS_RUN_SAMPLES"; }

# The run for this array task: line $SLURM_ARRAY_TASK_ID of paths.run_samples.
# Per-sample stages are submitted as  sbatch --array=1-$(grep -c . <run_samples file>) ...
task_sample(){
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        echo "not an array task; submit with --array=1-N (N = runs in $PATHS_RUN_SAMPLES)" >&2
        return 1
    fi
    local srr
    srr=$(run_samples | sed -n "${SLURM_ARRAY_TASK_ID}p")
    [[ -n "$srr" ]] || { echo "no run on line $SLURM_ARRAY_TASK_ID of $PATHS_RUN_SAMPLES" >&2; return 1; }
    echo "$srr"
}

# Sample-sheet column for a run: sample_field SRR6231181 location  (columns of paths.samples)
sample_field(){
    awk -F'\t' -v r="$1" -v c="$2" 'NR==1{for(i=1;i<=NF;i++) if($i==c) k=i; next}
                                    $1==r{print $k}' "$PATHS_SAMPLES"
}

# System mates: runs in paths.run_samples from the same location (system) as $1, itself
# included. Stage 4 maps every mate's reads to each assembly.
system_mates(){
    local loc
    loc=$(sample_field "$1" location)
    while read -r s; do
        [[ "$(sample_field "$s" location)" == "$loc" ]] && echo "$s"
    done < <(run_samples)
}

# Fail early with a clear message if an input from an earlier stage is missing
require(){
    local f
    for f in "$@"; do
        [[ -s "$f" ]] || { echo "missing input: $f (run the earlier stage first)" >&2; exit 1; }
    done
}

# Python helpers (tables, tagging, summaries) from src/plasmid_detection, run in pd-core
pdpy(){ PYTHONPATH="$PD_REPO/src" "$PD_CORE_ENV/bin/python" -m "$@"; }

# conda (de)activation with `set -u` relaxed: some envs' activate.d scripts (R, CheckM2,
# GTDB-Tk) read unset variables and would abort a `set -euo pipefail` script
activate(){ set +u; conda activate "$1"; set -u; }
deactivate(){ set +u; conda deactivate; set -u; }
