#!/bin/bash
#SBATCH --job-name=pd_data
#SBATCH --partition=Orion
#SBATCH --array=1-92%12
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%A_%a.err

# 01_download_reads.sbatch -- stage 1, one SRP121672 WGS run per array task:
#   prefetch -> vdb-validate -> fasterq-dump -> pigz -> manifest -> read-only
# Pair counts are checked against the NCBI SRA spot count (paths.sra_runinfo). The .sra
# and fasterq-dump temp files are removed once the run passes. A run whose manifest
# already says PASS is skipped, so re-submitting only redoes what's missing or failed.
# Combine the per-run manifests afterwards with 01b_collect_manifest.sbatch.

set -euo pipefail
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
conda activate "$ENVS_DATA"
THREADS=$PARAMS_DATA_THREADS

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PATHS_RUN_ACCESSIONS")
[[ -n "$SRR" ]] || { echo "no accession for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
SPOTS=$(awk -F, -v r="$SRR" 'NR>1 && $1==r {print $4}' "$PATHS_SRA_RUNINFO")
[[ -n "$SPOTS" ]] || { echo "$SRR not in $PATHS_SRA_RUNINFO" >&2; exit 1; }

R1="$PATHS_RAW_READS/${SRR}_1.fastq.gz"
R2="$PATHS_RAW_READS/${SRR}_2.fastq.gz"
MANIFEST="$PATHS_MANIFEST_DIR/${SRR}.tsv"
WORK="$PATHS_SRA_CACHE/$SRR"

echo "[$(date)] Starting ${SRR} (array task ${SLURM_ARRAY_TASK_ID}), SRA spots=${SPOTS}"
if [[ -f "$MANIFEST" && -f "$R1" && -f "$R2" ]] && grep -q $'\tPASS$' "$MANIFEST"; then
    echo "$SRR already PASS in $MANIFEST, skipping"
    exit 0
fi

# Start clean: earlier outputs for this run may be partial (and read-only)
chmod u+w "$R1" "$R2" 2>/dev/null || true
rm -rf "$R1" "$R2" "$MANIFEST" "$WORK"
mkdir -p "$PATHS_RAW_READS" "$PATHS_MANIFEST_DIR" "$WORK/tmp" "$WORK/fastq"
trap 'rm -rf "$WORK"' EXIT

# Per-task sra-tools settings (leaves ~/.ncbi alone; one file per task because concurrent
# vdb-config writes to a shared file race). Request full-quality .sra:
# SRA Lite replaces every quality score with one constant, which would break stage 2.
export NCBI_SETTINGS="$WORK/.ncbi/user-settings.mkfg"
mkdir -p "$(dirname "$NCBI_SETTINGS")"
vdb-config --simplified-quality-scores no > /dev/null

# --- prefetch + vdb-validate ---
prefetch --max-size u -O "$WORK" "$SRR"
SRA="$WORK/$SRR/$SRR.sra"
if [[ ! -f "$SRA" ]]; then
    echo "expected full-quality $SRA; prefetch got: $(ls "$WORK/$SRR")" >&2
    exit 1
fi
vdb-validate "$SRA"

# --- fasterq-dump ---
fasterq-dump --split-3 --skip-technical --threads "$THREADS" \
    --temp "$WORK/tmp" -O "$WORK/fastq" "$SRA"
produced=$(ls "$WORK/fastq")
if [[ "$produced" != "${SRR}_1.fastq"$'\n'"${SRR}_2.fastq" ]]; then
    echo "expected exactly ${SRR}_1.fastq and ${SRR}_2.fastq, got: $produced" >&2
    exit 1
fi

# Quality scores must vary (a constant score means SRA Lite slipped through)
nqual=$(head -n 40000 "$WORK/fastq/${SRR}_1.fastq" | awk 'NR%4==0' | fold -w1 | sort -u | wc -l)
if (( nqual < 3 )); then
    echo "only $nqual distinct quality characters in ${SRR}_1 -- simplified qualities?" >&2
    exit 1
fi

# --- pair counts vs SRA spot count ---
N1=$(awk 'END{print NR/4}' "$WORK/fastq/${SRR}_1.fastq")
N2=$(awk 'END{print NR/4}' "$WORK/fastq/${SRR}_2.fastq")
STATUS=PASS
[[ "$N1" == "$N2" && "$N1" == "$SPOTS" ]] || STATUS=FAIL

# --- pigz, md5, read-only ---
pigz -p "$THREADS" "$WORK/fastq/${SRR}_1.fastq" "$WORK/fastq/${SRR}_2.fastq"
mv "$WORK/fastq/${SRR}_1.fastq.gz" "$R1"
mv "$WORK/fastq/${SRR}_2.fastq.gz" "$R2"
MD5_1=$(md5sum "$R1" | cut -d' ' -f1)
MD5_2=$(md5sum "$R2" | cut -d' ' -f1)
chmod a-w "$R1" "$R2"

printf 'run\tsra_spots\tpairs_r1\tpairs_r2\tdistinct_quals\tr1\tr1_md5\tr2\tr2_md5\tstatus\n' > "$MANIFEST"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$SRR" "$SPOTS" "$N1" "$N2" "$nqual" \
    "$(basename "$R1")" "$MD5_1" "$(basename "$R2")" "$MD5_2" "$STATUS" >> "$MANIFEST"
cat "$MANIFEST"

echo "[$(date)] Finished ${SRR}: ${STATUS}"
[[ "$STATUS" == PASS ]]
