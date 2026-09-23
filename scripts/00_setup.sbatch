#!/bin/bash
#SBATCH --job-name=pd_setup
#SBATCH --partition=Orion
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.out
#SBATCH --error=/users/rgupta25/gibas/plasmid_detection/logs/%x_%j.err

# 00_setup.sbatch -- stage 0: conda envs, human reference + Bowtie2 index, databases.
#     sbatch scripts/00_setup.sbatch
# Every step is idempotent: a step is skipped once its .done marker (or env) exists,
# so a re-run only retries what failed. Ends by rewriting versions.txt and
# envs/lock/*.yml, which record every tool and database version actually installed.

set -uo pipefail

PD_REPO=${PD_REPO:-$HOME/gibas/plasmid_detection}
FAILED=()
sec(){ printf '\n==== [%s] %s\n' "$(date '+%F %T')" "$1"; }
fail(){ echo "  !! $1" >&2; FAILED+=("$1"); }
done_marker(){ [[ -e "$1/.done" ]]; }

module load mambaforge/23.11
source "$(conda info --base)/etc/profile.d/conda.sh"
export CONDA_PKGS_DIRS=$PD_REPO/conda/pkgs

# Build (or skip) one env from envs/<spec>.yml into prefix
make_env(){
    local spec=$1 prefix=$2
    if [[ -d "$prefix/conda-meta" ]]; then
        echo "  $(basename "$prefix") exists, skipping"
    else
        mamba env create -q -p "$prefix" -f "$PD_REPO/envs/$spec.yml" \
            || { fail "env $spec"; rm -rf "$prefix"; }
    fi
}

sec "1. Core env (config loader)"
make_env core "$PD_REPO/conda/envs/pd-core"
source "$PD_REPO/scripts/lib/common.sh" || { echo "cannot load config" >&2; exit 1; }
THREADS=$PARAMS_SETUP_THREADS
mkdir -p "$PATHS_CONDA_ENVS" "$PATHS_DB_ROOT" "$PATHS_WORK" "$PATHS_ENV_LOCKS"
echo "  repo:    $PATHS_REPO"
echo "  sw_root: $PATHS_SW_ROOT"
echo "  data:    $PATHS_DATA_ROOT"

sec "2. Stage envs"
for spec in data host genomad checkv checkm2 rgi amrfinderplus mobsuite gtdbtk; do
    make_env "$spec" "$PATHS_CONDA_ENVS/pd-$spec"
done

sec "3. Human reference $REFERENCE_HUMAN_NAME ($REFERENCE_HUMAN_ACCESSION) + Bowtie2 index"
if done_marker "$REFERENCE_HUMAN_DIR"; then
    echo "  exists, skipping"
else
    (
        set -e
        conda activate "$ENVS_HOST"
        mkdir -p "$REFERENCE_HUMAN_DIR" "$(dirname "$REFERENCE_HUMAN_BOWTIE2_INDEX")"
        cd "$REFERENCE_HUMAN_DIR"
        if [[ ! -s "$REFERENCE_HUMAN_FASTA" ]]; then
            datasets download genome accession "$REFERENCE_HUMAN_ACCESSION" \
                --include genome --filename ncbi_dataset.zip
            unzip -o -q ncbi_dataset.zip
            mv "ncbi_dataset/data/$REFERENCE_HUMAN_ACCESSION/"*_genomic.fna "$REFERENCE_HUMAN_FASTA"
            cp ncbi_dataset/data/assembly_data_report.jsonl .
            rm -rf ncbi_dataset ncbi_dataset.zip README.md md5sum.txt
        fi
        bowtie2-build --threads "$THREADS" "$REFERENCE_HUMAN_FASTA" "$REFERENCE_HUMAN_BOWTIE2_INDEX"
        md5sum "$REFERENCE_HUMAN_FASTA" > fasta.md5
        touch .done
    ) || fail "human reference / bowtie2 index"
fi

sec "4. geNomad database"
if done_marker "$DB_GENOMAD"; then echo "  exists, skipping"; else
    ( set -e; conda activate "$ENVS_GENOMAD"
      mkdir -p "$(dirname "$DB_GENOMAD")"
      genomad download-database "$(dirname "$DB_GENOMAD")"
      touch "$DB_GENOMAD/.done" ) || fail "geNomad db"
fi

sec "5. CheckV database"
if done_marker "$DB_CHECKV"; then echo "  exists, skipping"; else
    ( set -e; conda activate "$ENVS_CHECKV"
      root=$(dirname "$DB_CHECKV"); mkdir -p "$root"
      checkv download_database "$root"
      ln -sfn "$(basename "$(ls -d "$root"/checkv-db-v* | sort -V | tail -1)")" "$DB_CHECKV"
      touch "$DB_CHECKV/.done" ) || fail "CheckV db"
fi

sec "6. CheckM2 database"
if done_marker "$DB_CHECKM2"; then echo "  exists, skipping"; else
    ( set -e; conda activate "$ENVS_CHECKM2"
      mkdir -p "$DB_CHECKM2"
      checkm2 database --download --path "$DB_CHECKM2"
      dmnd=$(find "$DB_CHECKM2" -name '*.dmnd' | head -1)
      [[ -n "$dmnd" ]]
      conda env config vars set -p "$ENVS_CHECKM2" CHECKM2DB="$dmnd"
      touch "$DB_CHECKM2/.done" ) || fail "CheckM2 db"
fi

sec "7. CARD (RGI local database)"
if done_marker "$DB_CARD_DIR"; then echo "  exists, skipping"; else
    ( set -e; conda activate "$ENVS_RGI"
      mkdir -p "$DB_CARD_DIR"; cd "$DB_CARD_DIR"
      wget -q -O card-data.tar.bz2 "$DB_CARD_URL"
      tar -xjf card-data.tar.bz2
      rgi card_annotation -i card.json > card_annotation.log 2>&1
      rgi load --card_json card.json --card_annotation card_database_v*[0-9].fasta --local
      touch .done ) || fail "CARD db"
fi

sec "8. AMRFinderPlus database"
if done_marker "$DB_AMRFINDERPLUS"; then echo "  exists, skipping"; else
    ( set -e; conda activate "$ENVS_AMRFINDERPLUS"
      mkdir -p "$(dirname "$DB_AMRFINDERPLUS")"
      amrfinder_update -d "$(dirname "$DB_AMRFINDERPLUS")"
      touch "$DB_AMRFINDERPLUS/.done" ) || fail "AMRFinderPlus db"
fi

sec "9. MOB-suite database"
if done_marker "$DB_MOB_SUITE"; then echo "  exists, skipping"; else
    ( set -e; conda activate "$ENVS_MOBSUITE"
      mkdir -p "$DB_MOB_SUITE"
      mob_init -d "$DB_MOB_SUITE"
      touch "$DB_MOB_SUITE/.done" ) || fail "MOB-suite db"
fi

sec "10. GTDB $DB_GTDB_RELEASE (GTDB-Tk package, ~141 GB compressed)"
if done_marker "$DB_GTDB_DIR"; then echo "  exists, skipping"; else
    (
        set -e
        # Tarball goes to scratch (resumable, md5-checked) and is removed once extracted
        tarball="$PATHS_WORK/downloads/$(basename "$DB_GTDB_URL")"
        mkdir -p "$(dirname "$tarball")" "$DB_GTDB_DIR"
        wget -q -c -O "$tarball" "$DB_GTDB_URL"
        echo "$DB_GTDB_MD5  $tarball" | md5sum -c -
        tar -xzf "$tarball" -C "$DB_GTDB_DIR" --strip-components=1
        conda env config vars set -p "$ENVS_GTDBTK" GTDBTK_DATA_PATH="$DB_GTDB_DIR"
        touch "$DB_GTDB_DIR/.done"
        rm -f "$tarball"
    ) || fail "GTDB $DB_GTDB_RELEASE"
fi

sec "11. Provenance: versions.txt + envs/lock/"
pkgver(){ conda list -p "$1" "^$2\$" 2>/dev/null | awk -v p="$2" '$1==p{print $2}'; }
for prefix in "$PATHS_CONDA_ENVS"/pd-*; do
    [[ -d "$prefix/conda-meta" ]] || continue
    conda env export -p "$prefix" | grep -v '^prefix:' > "$PATHS_ENV_LOCKS/$(basename "$prefix").yml"
done
{
    echo "# Written by scripts/00_setup.sbatch -- do not edit by hand"
    echo "date: $(date -Iseconds)"
    echo "host: $(hostname)"
    echo "conda: $(conda --version | awk '{print $2}')  mamba: $(mamba --version | head -1 | awk '{print $2}')"
    echo
    echo "## Tools (full solved envs in envs/lock/)"
    for pair in core:python core:pyyaml data:sra-tools data:pigz host:ncbi-datasets-cli host:bowtie2 host:samtools \
                genomad:genomad checkv:checkv checkm2:checkm2 rgi:rgi \
                amrfinderplus:ncbi-amrfinderplus mobsuite:mob_suite gtdbtk:gtdbtk; do
        env=${pair%%:*} pkg=${pair#*:}
        echo "$pkg: $(pkgver "$PATHS_CONDA_ENVS/pd-$env" "$pkg")  (pd-$env)"
    done
    echo
    echo "## Reference / databases"
    echo "human: $REFERENCE_HUMAN_NAME $REFERENCE_HUMAN_ACCESSION  fasta_md5=$(cut -d' ' -f1 "$REFERENCE_HUMAN_DIR/fasta.md5" 2>/dev/null)"
    echo "genomad_db: $(cat "$DB_GENOMAD/version.txt" 2>/dev/null)"
    echo "checkv_db: $(readlink "$DB_CHECKV" 2>/dev/null)"
    echo "checkm2_db: $(basename "$(find "$DB_CHECKM2" -name '*.dmnd' 2>/dev/null | head -1)")"
    echo "card: $("$ENVS_CORE/bin/python" -c 'import json,sys; print(json.load(open(sys.argv[1]))["_version"])' "$DB_CARD_DIR/card.json" 2>/dev/null)"
    echo "amrfinderplus_db: $(cat "$DB_AMRFINDERPLUS/version.txt" 2>/dev/null)"
    # mob_init writes no version file: record its download date + a checksum of the plasmid DB
    echo "mob_suite_db: $(grep -o 'Download date: [0-9-]*' "$DB_MOB_SUITE/status.txt" 2>/dev/null)  ncbi_plasmid_full_seqs_md5=$(md5sum "$DB_MOB_SUITE/ncbi_plasmid_full_seqs.fas" 2>/dev/null | cut -d' ' -f1)"
    echo "gtdb: $DB_GTDB_RELEASE  tarball_md5=$DB_GTDB_MD5  $(grep -h VERSION_DATA "$DB_GTDB_DIR"/metadata/metadata.txt 2>/dev/null)"
} > "$PATHS_VERSIONS"
cat "$PATHS_VERSIONS"

sec "Done"
if (( ${#FAILED[@]} )); then
    printf '  FAILED: %s\n' "${FAILED[@]}"
    exit 1
fi
echo "  all steps complete"
