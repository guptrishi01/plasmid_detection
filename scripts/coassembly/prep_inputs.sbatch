#!/bin/bash
#SBATCH --job-name=pd_coasm_prep
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# prep_inputs.sbatch -- inputs of the co-assembly workflow (nextflow/coassembly.nf), taken
# from the finished 92-run Nextflow run without re-running it:
#   <coasm_inputs>/hostfree/<run>_{1,2}.fastq.gz   host-free reads (HOST_DEPLETION outputs)
#   <coasm_inputs>/site_mags/<group>/*.fa          per-site MAGs (DREP_GROUP, 11 site groups)
# Both exist only in the Nextflow work dir (results_nf keeps logs and tables), so they are
# hard-linked out: no extra space, and they survive a cleanup of work_nf. Task work dirs come
# from `nextflow log <run name>` (default: the full 92-run run, 2026-11-08).
# Checks: every run's R1 and R2 read counts equal retained_pairs in tables/03_host_depletion.tsv,
# and every group's MAG files are exactly that group's group MAGs in tables/07_group_mags.tsv.
#     sbatch scripts/coassembly/prep_inputs.sbatch [nextflow run name]
set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
RUN_NAME=${1:-suspicious_murdock}
N=$PATHS_NF_RESULTS
P=$PATHS_COASM_INPUTS
mkdir -p "$P/hostfree" "$P/site_mags"
export NXF_HOME=${NXF_HOME:-$PATHS_DATA_ROOT/nxf_home}
cd "$PD_REPO"   # nextflow log reads the run history in ./.nextflow
"$PATHS_CONDA_ENVS/pd-nextflow/bin/nextflow" log "$RUN_NAME" -f name,status,workdir > "$P/nextflow_log.tsv"

# ---- host-free reads ------------------------------------------------------------------------
awk -F'\t' '$1 ~ /^HOST_DEPLETION / && ($2 == "COMPLETED" || $2 == "CACHED") {
    r = $1; gsub(/.*\(|\)/, "", r); print r "\t" $3 }' "$P/nextflow_log.tsv" | sort -u > "$P/hostfree_workdirs.tsv"
dups=$(cut -f1 "$P/hostfree_workdirs.tsv" | sort | uniq -d | wc -l)
(( dups == 0 )) || { echo "runs with more than one finished HOST_DEPLETION task: $dups" >&2; exit 1; }
n=$(wc -l < "$P/hostfree_workdirs.tsv")
expected=$(awk 'NR>1' "$N/tables/03_host_depletion.tsv" | wc -l)
(( n == expected )) || { echo "HOST_DEPLETION tasks: $n, runs in 03_host_depletion.tsv: $expected" >&2; exit 1; }
while IFS=$'\t' read -r run dir; do
    for m in 1 2; do
        src=$dir/hostfree/${run}_$m.fastq.gz dst=$P/hostfree/${run}_$m.fastq.gz
        [[ -s $src ]] || { echo "missing $src" >&2; exit 1; }
        ln -f "$src" "$dst"
        [[ $(stat -c %i "$src") == $(stat -c %i "$dst") ]] || { echo "not a hard link: $dst" >&2; exit 1; }
    done
done < "$P/hostfree_workdirs.tsv"

activate "$ENVS_ASSEMBLY"   # seqkit
seqkit stats -T -j "${SLURM_CPUS_PER_TASK:-4}" "$P"/hostfree/*.fastq.gz > "$P/hostfree_seqkit.tsv"
deactivate
python3 - "$P/hostfree_seqkit.tsv" "$N/tables/03_host_depletion.tsv" <<'EOF'
import csv, os, sys
seq = {os.path.basename(r["file"]): int(r["num_seqs"]) for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
bad = 0
for r in csv.DictReader(open(sys.argv[2]), delimiter="\t"):
    want = int(r["retained_pairs"])
    got = (seq.get(f"{r['run']}_1.fastq.gz"), seq.get(f"{r['run']}_2.fastq.gz"))
    if got != (want, want):
        bad += 1
        print(f"MISMATCH {r['run']}: R1/R2 {got}, retained_pairs {want}")
print(f"host-free reads: {len(seq) // 2} runs, {sum(seq.values()) // 2:,} pairs, {bad} mismatches")
sys.exit(1 if bad else 0)
EOF

# ---- per-site MAGs (level-1 dRep, one dir per site group) -------------------------------------
awk -F'\t' '$1 ~ /^DREP_GROUP / && ($2 == "COMPLETED" || $2 == "CACHED") {
    g = $1; gsub(/.*\(|\)/, "", g); print g "\t" $3 }' "$P/nextflow_log.tsv" | sort -u > "$P/site_mag_workdirs.tsv"
while IFS=$'\t' read -r g dir; do
    rm -rf "${P:?}/site_mags/$g"; mkdir -p "$P/site_mags/$g"
    ln -f "$dir"/drep/dereplicated_genomes/*.fa "$P/site_mags/$g/"
done < "$P/site_mag_workdirs.tsv"
python3 - "$P/site_mags" "$N/tables/07_group_mags.tsv" <<'EOF'
import csv, os, sys
from collections import defaultdict
want = defaultdict(set)
for r in csv.DictReader(open(sys.argv[2]), delimiter="\t"):
    want[r["group"]].add(r["group_mag"])
bad, total = 0, 0
for g in sorted(want):
    d = os.path.join(sys.argv[1], g)
    got = {f[:-3] for f in os.listdir(d) if f.endswith(".fa")} if os.path.isdir(d) else set()
    total += len(got)
    if got != want[g]:
        bad += 1
        print(f"MISMATCH {g}: {len(got)} files, {len(want[g])} group MAGs in 07_group_mags.tsv")
print(f"per-site MAGs: {len(want)} groups, {total} genomes, {bad} mismatches")
sys.exit(1 if bad else 0)
EOF
echo "[$(date)] done: $P"
