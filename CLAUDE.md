# CLAUDE.md

Project-specific context for this repo. See [README.md](README.md) for the project summary.

## What this is

Plasmid detection in the 92 shotgun metagenomes (WGS) of SRA study **SRP121672**
(BioProject PRJNA415974; Gibas Lab), to study plasmid–host interactions in
environments under antibiotic stress. The data come from:

> Lambirth K., Tsilimigras M., Lulla A., … Fodor A., Gibas C. (2018). *Microbial Community
> Composition and Antibiotic Resistance Genes within a North Carolina Urban Water System.*
> Water 10(11):1539. doi:10.3390/w10111539

The PDF lives at `/users/rgupta25/gibas/` (outside the repo; `*.pdf` is gitignored anyway).
The paper did **no assembly and no plasmid analysis** — its methods were Trimmomatic + PEAR,
MetaPhlAn2, ShortBRED (CARD 1.1.0 + Lahey β-lactamases), HUMAnN2, and `lm` in R with BH
correction. Plasmid work here is new, not a reproduction, so there's no published number
to diff against for the plasmid results themselves; the paper's ARG findings (ShortBRED
RPKM by site) are the nearest ground truth for sanity-checking ARG calls.

## The 92 samples

`metadata/srp121672_wgs_samples.tsv` is the sample sheet: one row per run with
`location`, `site`, `timepoint` parsed from NCBI SRA `SampleName`
(`Sample_<Location>_<Site>_<TP>`). Verified counts:

- Time points 1–3: 22 runs each (11 Mallard Creek + 11 Sugar Creek); time point 4: 26
  (the same 22 + 4 remote controls MNT_A/B, UW_A/B).
- Mallard Creek sites: INF, PCI, PCE, ATE, FCE, HOSP, RES, UP_A/B, DS_A/B.
  Sugar Creek has UV instead of PCI.
- The paper sequenced triplicates (66–78 libraries per time point); **SRA holds one run
  per site × time point** — consistent with the paper keeping the most deeply sequenced
  replicate as representative. Don't go looking for the missing replicates.
- All 92: Illumina HiSeq 2500, paired-end 125 bp, 15.1–30.9 M read pairs (NCBI SRA
  `spots`), ~163 GB of `.sra` total. The other 496 runs in SRP121672 are 16S amplicon — out of scope.

## Stage 1 follows the spec exactly

The workflow spec is the source of truth: `prefetch` → `vdb-validate` → `fasterq-dump` →
`pigz`, a manifest per run (pair count, md5), pair counts checked against SRA spot counts
(`metadata/srp121672_wgs_sra_runinfo.csv`), raw FASTQ read-only. **Don't substitute ENA
or any other mirror** — an ENA-based download was built once and deleted for not
following the spec. Two things the spec doesn't say but the data forced:

- **SRA Lite.** NCBI's default download for these runs is `.sralite` (one constant
  quality score). The full-quality `.sra` exists (AWS ODP bucket). `01_download_reads`
  sets `--simplified-quality-scores no` in a project-local `NCBI_SETTINGS` file (never
  touch `~/.ncbi`), fails if prefetch returns anything but `<run>.sra`, and records the
  number of distinct quality characters in the manifest (must be ≥3).
- **`.sra` files are deleted** after a run PASSes (the spec doesn't ask to keep them;
  re-running the task re-fetches).

## Never assume — verify and prove it

- Prefer a direct check over reasoning: run the thing against real data and look at the
  output. A claim like "the DB installed" means the version string is in `versions.txt`,
  not that the command exited 0.
- When a check surfaces a difference, root-cause it before making it disappear.
- State what was actually verified ("92/92 runs PASS in raw_manifest.tsv") rather than what
  should follow from it.

## Cluster environment

Two roots, both set in `config/config.yaml`:

- `paths.sw_root` (default = this repo, `/users/rgupta25/gibas/plasmid_detection`, home,
  500G, not purged) — `conda/envs/`, `conda/pkgs/`, `db/`. Built by `00_setup.sbatch`.
- `paths.data_root` (`/scratch/rgupta25/plasmid_detection`, 5T, **assume it can be
  purged**) — raw reads, intermediates, results. `data/raw` in the repo is a symlink to
  `$data_root/raw`. Raw reads are re-fetchable with `01_download_reads.sbatch`.

Modules: `mambaforge/23.11` (conda for every env). The cluster also has `gtdbtk/2.4.0`,
`bowtie2`, `samtools` modules — **don't use them**; every tool comes from its pinned env.
All compute runs via `sbatch` on the `Orion` partition; compute nodes have internet.

## Config — one file, every path and parameter

`config/config.yaml` is the only place paths/parameters live. `${dotted.key}` references
are resolved by `src/plasmid_detection/config.py`. Every batch script starts with:

```bash
source "${PD_REPO:-$HOME/gibas/plasmid_detection}/scripts/lib/common.sh"
```

which loads conda and exports every key as `UPPER_SNAKE` (`db.gtdb.dir` → `DB_GTDB_DIR`,
`envs.genomad` → `ENVS_GENOMAD`). New stages add their parameters under `params.<stage>`,
never as literals in the script. Two exceptions, both unavoidable:

- `#SBATCH` lines are parsed before any code runs, so log paths/partition are literal there.
- `common.sh` must find the `pd-core` env before it can read the config, so that env's
  location is fixed at `$PD_REPO/conda/envs/pd-core`.

Never add a top-level config section that would export `SLURM_*` (Slurm reads those as
input); `tests/test_config.py` enforces this.

## Environments — one per stage, pinned

`envs/<name>.yml` pins the tool version (conda-forge + bioconda, `nodefaults`). "One env
per stage" in practice: an env per stage where the tools co-solve, split per tool where
they don't or where stage 0 already built one (SCAPP needs Python 3.7; stage 5 uses the
geNomad/CheckV/MOB-suite/AMRFinderPlus/RGI envs from stage 0).
`00_setup.sbatch` builds each into `conda/envs/pd-<stage>` and exports the fully solved
env to `envs/lock/pd-<stage>.yml` (tracked). To change a version: edit the spec, delete
the built env, re-run setup, commit spec + lock + `versions.txt` together.

| env | tools | database |
|---|---|---|
| pd-core | python, pyyaml (config loader, `tables.py`) | — |
| pd-data | sra-tools, pigz | — (stage 1) |
| pd-qc | FastQC, fastp, MultiQC (stage 2; MultiQC also for stage 3 logs) | — |
| pd-host | ncbi-datasets-cli, bowtie2, samtools | GRCh38.p14 (GCF_000001405.40) + Bowtie2 index (`.bt2`, not `.bt2l`) |
| pd-assembly | MEGAHIT, SPAdes, Bowtie2, samtools, QUAST, Nonpareil, seqkit (stage 4; seqkit for stage 9 depth) | — |
| pd-scapp | SCAPP (own env: pins Python 3.7) | — |
| pd-coverm | CoverM (stages 5, 7) | — |
| pd-binning | MetaBAT2 + `jgi_summarize_bam_contig_depths` | — |
| pd-drep | dRep | — |
| pd-stats | R, ALDEx2, ANCOMBC, vegan, ggplot2 (stage 9) | — |
| pd-genomad | geNomad | `db/genomad/genomad_db` |
| pd-checkv | CheckV | `db/checkv/current` → `checkv-db-vX.Y` |
| pd-checkm2 | CheckM2 | `db/checkm2` (`CHECKM2DB` set in env) |
| pd-rgi | RGI | CARD, `db/card/localDB` — rgi `--local` reads `./localDB`, so run rgi with cwd `db/card` or copy/link `localDB` into the work dir |
| pd-amrfinderplus | AMRFinderPlus | `db/amrfinderplus/latest` |
| pd-mobsuite | MOB-suite | `db/mob_suite` |
| pd-gtdbtk | GTDB-Tk **2.6.1** | GTDB **r226** (`GTDBTK_DATA_PATH` set in env) |

**GTDB-Tk is pinned by the GTDB release, not by "latest".** R226 supports GTDB-Tk
2.4.1–2.6.1; 2.7.x requires R232. Don't bump gtdbtk without deciding to move GTDB too.

## Provenance

`versions.txt` (tracked, rewritten by `00_setup.sbatch` step 11) records every tool
version and every database version actually installed. Setup steps are idempotent via
`.done` markers in each DB dir; a re-run only retries what failed and always refreshes
`versions.txt` and `envs/lock/`. After any setup run, check `versions.txt` has no blank
fields before trusting it.

## Pipeline execution order

`RUNNING.md` (gitignored, local) is the run log: current stage, next command, job ids.
Read it first and update it after each step. Keep replies to the user to ~20-25 short
lines: they work in a PuTTY terminal.


`scripts/NN*_*.sbatch` run in numeric/letter order; the full list with what each does is
in README.md ("Workflow" + "Running the pilot"). The user submits every stage by hand,
one step at a time, and tunes `params.<stage>` between runs — **don't chain or launch
pipeline stages unless asked.** Stage 0/setup jobs and tool `--help` checks are fine.

Conventions every stage script follows (keep them when adding or editing one):

- Starts with `set -euo pipefail` and `source .../scripts/lib/common.sh`.
- Activates envs with the `activate`/`deactivate` helpers from `common.sh`, never bare
  `conda activate`: some envs' activate.d scripts read unset variables and abort under
  `set -u`.
- Samples come from `paths.run_samples` (`metadata/pilot_samples.txt` for the pilot; point
  it at `srp121672_wgs_run_accessions.txt` for all 92). Per-run scripts have **no**
  `#SBATCH --array`; they're submitted as `sbatch --array=1-N` and read their run with
  `task_sample`. Cross-sample steps are single jobs.
- "System mates" (stage 4 cross-mapping) = runs in `paths.run_samples` sharing the sample
  sheet's `location` (`system_mates` in `common.sh`).
- Inputs from earlier stages are checked with `require` (fails with "run the earlier stage
  first"). Trimmed/host-depleted reads and BAMs are written to a temp name and moved into
  place, so a killed job doesn't leave a complete-looking read file; tool output dirs
  (MEGAHIT, geNomad, ...) are removed at the start of a re-run instead.
- Tables are built by `src/plasmid_detection/tables.py` (stdlib only, run via `pdpy` in
  pd-core); R lives in `scripts/R/`. Tool output formats those parsers rely on were
  checked against the installed tools' source (AMRFinderPlus `Contig id`; RGI `Contig`
  already has the ORF suffix stripped; geNomad `*_summary/*_{plasmid,virus}_summary.tsv`,
  proviruses as `contig|provirus_start_end`; `samtools bedcov -c` = 4 BED columns, depth
  per BAM, then read count per BAM). Re-check when bumping a tool version.

**Never run `find . -name __pycache__ -exec rm -rf {} +` (or any recursive delete) from
the repo root**: `conda/pkgs/` lives inside the repo, and deleting `__pycache__` there
corrupted the package cache once (CondaVerificationError on the next env build; fixed
by removing the extracted package dirs, keeping the archives). Clean `src/ tests/` only.

## Git scope

`results/` in the repo is tracked on purpose (user decision 2026-09-29: shown in lab
meetings). It is filled only by `scripts/sync_results.sh`, which copies reports/logs/tables
from scratch and excludes sequence/alignment files and anything > 20 MB. Never copy
assemblies, BAMs, or reads into it by hand.


Remote: `git@github.com:guptrishi01/plasmid_detection.git` (SSH key `~/.ssh/id_rsa`).
**Don't push unless asked.** `.gitignore` excludes `data/`, `conda/`, `db/`, `logs/`,
FASTQ/FASTA/BAM/index files, and `*.pdf`. Never force-add from those — anything rebuildable
belongs in `00_setup.sbatch`, not committed. Metadata TSVs, env specs + locks, config,
and `versions.txt` are tracked.

## Python

`src/plasmid_detection/` is the Python: `config.py` (config loader; lists flatten to
space-separated strings, and `${...}` inside a YAML `[..]` list must be quoted) and
`tables.py` (per-stage tables, AMR tagging, plasmid-host inference, stdlib only).
`pyyaml` is the sole dependency (`requirements.txt` / `pyproject.toml`, matching
`envs/core.yml`). Tests:
`pytest` from the repo root (`pythonpath = ["src"]` is set in `pyproject.toml`); lint:
`ruff check . && ruff format --check src tests`.
