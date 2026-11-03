#!/bin/bash
# submit_group_binning.sh -- queue the group-binning re-run (stages 4c -> 9e) in one go, each
# job held by Slurm until the jobs it reads from finish OK (--dependency=afterok). This
# script only submits and exits (seconds); nothing keeps running on the login node.
# If a job fails, everything downstream of it is cancelled (--kill-on-invalid-dep=yes).
# Runs: paths.run_samples (pilot or all 92). Job ids are printed and appended to RUNNING.md.
#   bash scripts/submit_group_binning.sh
#
# Order (what each job waits for):
#   04c map ─┬─ 05d coverage/tag ─┬─ 06c BAM cleanup (also waits for 06a)
#            │                    ├─ 09d rollup ── 09e stats
#            └─ 06a binning ── 06b CheckM2 ── 07a dRep ─┬─ 07b GTDB-Tk ─┬─ 07c abundance + tables
#                                                       │               └─ 08b hosts (also 05d)
#                                                       └─ 08a MAG AMR

set -euo pipefail
PD_REPO=${PD_REPO:-$HOME/gibas/plasmid_detection}
cd "$PD_REPO"
RUNS=$(PYTHONPATH="$PD_REPO/src" "$PD_REPO/conda/envs/pd-core/bin/python" -m plasmid_detection.config get paths.run_samples)
N=$(grep -c . "$RUNS")
S="$PD_REPO/scripts"

# sub <afterok job ids, comma-separated or ""> <sbatch args...>  -> prints the job id
sub(){
    local dep=$1; shift
    if [[ -n "$dep" ]]; then
        sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:${dep//,/:}" "$@"
    else
        sbatch --parsable "$@"
    fi
}

map=$(sub ""            --array=1-"$N" "$S/04c_map.sbatch")
cov=$(sub "$map"        --array=1-"$N" "$S/05d_coverage_tag.sbatch")
bin=$(sub "$map"        --array=1-"$N" "$S/06a_binning.sbatch")
clean=$(sub "$cov,$bin" "$S/06c_cleanup_bams.sbatch")
chk=$(sub "$bin"        --array=1-"$N" "$S/06b_checkm2_mimag.sbatch")
drep=$(sub "$chk"       "$S/07a_drep.sbatch")
gtdb=$(sub "$drep"      "$S/07b_gtdbtk.sbatch")
abund=$(sub "$gtdb"     "$S/07c_mag_abundance.sbatch")
magamr=$(sub "$drep"    "$S/08a_mag_amr.sbatch")
hosts=$(sub "$gtdb,$cov" "$S/08b_plasmid_hosts.sbatch")
roll=$(sub "$cov"       "$S/09d_rollup.sbatch")
stats=$(sub "$roll"     "$S/09e_stats.sbatch")

summary="04c $map | 05d $cov | 06a $bin | 06c $clean | 06b $chk | 07a $drep | 07b $gtdb | 07c $abund | 08a $magamr | 08b $hosts | 09d $roll | 09e $stats"
echo "runs: $N ($RUNS)"
echo "$summary"
echo "- $(date +%F) group-binning re-run submitted ($N runs, afterok chain): $summary" >> "$PD_REPO/RUNNING.md"
