#!/bin/bash
#SBATCH --job-name=pd_nextflow
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# submit.sbatch -- the Nextflow head job: stays running and submits every pipeline task to
# Slurm itself (at most 20 queued at once). Arguments go to nextflow run, e.g.
#   sbatch nextflow/submit.sbatch -profile slurm,pilot -resume
#   sbatch nextflow/submit.sbatch -profile slurm -resume
set -euo pipefail
module load mambaforge/23.11
bash "${PD_REPO:-$HOME/gibas/plasmid_detection}/nextflow/run.sh" "$@"
