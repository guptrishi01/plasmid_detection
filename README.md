# plasmid_detection

Gibas Lab project. This project takes in metagenomic sequences and extracts plasmid
sequences in order to better understand the interactions between plasmids and their
bacterial hosts in environments under antibiotic stress.

Data source:

> Lambirth K, Tsilimigras M, Lulla A, Johnson J, Al-Shaer A, Wynblatt O, Sypolt S,
> Brouwer C, Clinton S, Keen O, Redmond M, Fodor A, Gibas C. **"Microbial Community
> Composition and Antibiotic Resistance Genes within a North Carolina Urban Water
> System."** *Water* 10(11):1539 (2018). doi:[10.3390/w10111539](https://doi.org/10.3390/w10111539)

Project-specific working rules (config, envs, verification, git scope) are in
[CLAUDE.md](CLAUDE.md).

## Objective

Recover plasmids (and, alongside them, viruses and chromosomal contigs) from the 92
shotgun metagenomes of an urban water system — two activated-sludge wastewater
treatment plants, their sewer trunklines, and the creeks they discharge into — then:

1. place antibiotic-resistance genes (ARGs) on plasmid, virus, or chromosome contigs;
2. link plasmids to candidate bacterial hosts (MAGs);
3. relate plasmid and mobile-resistome patterns to plant stage, system, time point, and
   the antibiotic concentrations measured on each sample.

**What the source paper did and didn't do.** The paper profiled the same samples with
Trimmomatic + PEAR (read QC/merging), MetaPhlAn2 2.5.0 (taxonomy), ShortBRED against a
custom CARD 1.1.0 + Lahey β-lactamase marker set (ARG quantification, RPKM), HUMAnN2
(pathways), and per-feature linear models in R with Benjamini-Hochberg correction. It did
**no assembly, binning, or plasmid analysis** — everything from stage 4 onward here is new,
not a reproduction. The paper's ShortBRED ARG results are the nearest published reference
point (see stage 9's read-level vs. contig-level comparison).

## Dataset

SRA study **SRP121672** (BioProject **PRJNA415974**). The study holds 588 runs; the **92
WGS runs** are in scope (the other 496 are `AMPLICON`). All 92: Illumina HiSeq 2500,
paired-end 2 × 125 bp, 15.05–30.86 M read pairs per run (NCBI SRA `spots`), 162.7 GB of
`.sra` in total (NCBI `size_MB`).

| System | Sites | Runs |
|---|---|---|
| Mallard Creek (MC) | INF, PCI, PCE, ATE, FCE, HOSP, RES, UP_A, UP_B, DS_A, DS_B (11) | 44 (11 × 4 time points) |
| Sugar Creek (SC) | INF, PCE, ATE, FCE, UV, HOSP, RES, UP_A, UP_B, DS_A, DS_B (11) | 44 (11 × 4 time points) |
| Mountain control | MNT_A, MNT_B | 2 (time point 4 only) |
| Uwharrie control | UW_A, UW_B | 2 (time point 4 only) |

Sites: INF raw influent; PCI/PCE primary clarifier influent/effluent; ATE aeration tank
effluent; FCE final clarifier effluent; UV UV-disinfected effluent (SC only; PCI was not
sampleable at SC, UV not at MC); HOSP/RES hospital-adjacent/residential sewer trunklines;
UP_A/B and DS_A/B two creek sites upstream and downstream of effluent release. Time points
(2016): 1 late winter, 2 early spring, 3 late spring, 4 mid-summer — 22 runs each at 1–3,
26 at 4.

The paper sequenced triplicates; SRA holds **one run per site × time point**, consistent
with the paper keeping the most deeply sequenced replicate as representative.

Each BioSample carries the paper's LC-MS antibiotic measurements as attributes, in ng/L
(`ND` = not detected): ertapenem, amoxicillin, ciprofloxacin, doxycycline, azithromycin,
clindamycin, sulfamethoxazole, cephalexin, trimethoprim, levofloxacin — plus
`collection_date`, `lat_lon`, `temp`, `ph` (checked on SAMN07839266 and SAMN07839313).

Metadata in `metadata/`:

| File | Contents |
|---|---|
| `srp121672_wgs_samples.tsv` | Sample sheet: run, BioSample, `SampleName`, library, SRA spots, parsed `location` / `site` / `timepoint` |
| `srp121672_wgs_sra_runinfo.csv` | NCBI SRA run info for the 92 runs (spot counts used by stage 1) |
| `raw_manifest.tsv` | Stage 1 output: per run pair counts vs SRA spots, MD5s, PASS/FAIL |
| `srp121672_wgs_only.tsv`, `srp121672_wgs_run_accessions.txt` | The 92 WGS runs |
| `srp121672_all_metadata.tsv` | All 588 runs in the study |

## Repo structure

```
plasmid_detection/
├── config/config.yaml        # The one config: every path and parameter
├── envs/                     # One pinned conda env spec per stage (<stage>.yml)
│   └── lock/                 #   fully solved envs exported by 00_setup (tracked)
├── scripts/
│   ├── lib/common.sh         # Sourced by every stage: conda + config -> env vars
│   ├── 00_setup.sbatch       # Stage 0: envs, human reference + index, databases
│   ├── 01_download_reads.sbatch  # Stage 1, 92-task array: prefetch → vdb-validate →
│   │                             #   fasterq-dump → pigz → per-run manifest → read-only
│   └── 01b_collect_manifest.sbatch  # Combine + check the 92 per-run manifests
├── src/plasmid_detection/    # Config loader (Python, pyyaml)
├── tests/                    # pytest for the config loader
├── metadata/                 # Sample sheet + NCBI SRA metadata + raw manifest (above)
├── versions.txt              # Every tool + database version installed (written by 00_setup)
├── requirements.txt, pyproject.toml
└── data/raw -> /scratch/rgupta25/plasmid_detection/raw   (gitignored)
```

Not tracked in git: `conda/` (envs + package cache), `db/` (reference + databases),
`data/`, `logs/`, sequence/index files, `*.pdf`.

**Storage split.** `paths.sw_root` = this repo in home (500G, not purged) holds `conda/`
and `db/`. `paths.data_root` = `/scratch/rgupta25/plasmid_detection` (5T, purgeable) holds
raw reads, intermediates, and results. **Config.** Every stage sources
`scripts/lib/common.sh`, which exports `config/config.yaml` as env vars
(`db.gtdb.dir` → `DB_GTDB_DIR`); stage parameters go under `params.<stage>`, never as
literals. **Envs.** One conda env per stage, `envs/<stage>.yml` pins the tool version;
the solved env is exported to `envs/lock/`. Details in [CLAUDE.md](CLAUDE.md).

## Workflow

Status legend: **implemented** = script exists and has run; **planned** = env/database
installed, stage script not written; **not yet installed** = tool not yet in `envs/` or
`00_setup`.

### 0. Setup — implemented (`scripts/00_setup.sbatch`)

Git repository for all code. One config file for paths and parameters
(`config/config.yaml`). One conda environment per stage, versions pinned in a yml.
Human reference GRCh38.p14 (GCF_000001405.40) via NCBI Datasets CLI, with a Bowtie2
index. Databases: geNomad, CheckV, CheckM2, CARD, AMRFinderPlus, MOB-suite, GTDB r226.
Every version used is recorded in `versions.txt`.

<!-- STAGE0_VERSIONS -->
Installed (from `versions.txt`, setup jobs 27008574 + 27010950, 2026-09-22):

| Component | Version | Env / location |
|---|---|---|
| python / pyyaml | 3.12.3 / 6.0.3 | pd-core |
| sra-tools / pigz | 3.4.1 / 2.8 | pd-data |
| ncbi-datasets-cli / bowtie2 / samtools | 18.37.0 / 2.5.5 / 1.24 | pd-host |
| geNomad (db) | 1.12.0 (db 1.9) | pd-genomad |
| CheckV (db) | 1.1.1 (checkv-db-v1.5) | pd-checkv |
| CheckM2 (db) | 1.1.0 (uniref100.KO.1) | pd-checkm2 |
| RGI (CARD) | 6.0.8 (CARD 4.0.2) | pd-rgi |
| AMRFinderPlus (db) | 4.2.7 (2026-08-07.1) | pd-amrfinderplus |
| MOB-suite (db) | 3.1.9 (downloaded 2026-09-22) | pd-mobsuite |
| GTDB-Tk (GTDB) | 2.6.1 (r226, tarball MD5 verified) | pd-gtdbtk |
| Human reference | GRCh38.p14, GCF_000001405.40 + Bowtie2 index | `db/human/GRCh38.p14` |
<!-- /STAGE0_VERSIONS -->

GTDB-Tk is pinned to 2.6.1 because GTDB r226 supports GTDB-Tk 2.4.1–2.6.1 only; the
2.7.x series requires r232.

### 1. Data — implemented (`scripts/01_download_reads.sbatch`, `01b_collect_manifest.sbatch`)

sra-tools: `prefetch`, `vdb-validate`, `fasterq-dump`. `pigz` to compress. Manifest per
run: pair count, md5. Check pair counts against SRA spot counts. Raw FASTQ read-only.

Per run (array task): `prefetch` the **full-quality** `.sra` (sra-tools is set, through a
project-local `NCBI_SETTINGS` file, to prefer full base qualities over SRA Lite, whose
constant quality scores would invalidate stage 2) → `vdb-validate` →
`fasterq-dump --split-3 --skip-technical` (must yield exactly `_1` and `_2`) → count pairs
in R1 and R2 and compare with the NCBI SRA spot count → `pigz` → MD5 → `chmod a-w`. The
`.sra` and temp files are deleted once the run is done. Each run writes
`$data_root/manifest/<run>.tsv` (SRA spots, R1/R2 pairs, distinct quality characters,
file MD5s, PASS/FAIL); `01b_collect_manifest.sbatch` combines them into
`metadata/raw_manifest.tsv` and checks all 92 are PASS and read-only.

| Tools | Env | Status |
|---|---|---|
| sra-tools 3.4.1, pigz 2.8 | pd-data | implemented |

### 2. Read QC — planned

FastQC on raw reads. fastp. FastQC on trimmed reads. MultiQC over everything.

Table: raw pairs, trimmed pairs, Q30, duplication, insert size, per sample.

| Tools | Env | Status |
|---|---|---|
| FastQC, fastp, MultiQC | — | not yet installed |

### 3. Host depletion — planned

Bowtie2 against the human index. Keep pairs where neither mate maps.

Table: raw, trimmed, human, retained pairs per sample. MultiQC of the Bowtie2 logs.

| Tools | Env | Status |
|---|---|---|
| Bowtie2, samtools; GRCh38.p14 index | `pd-host` | installed by stage 0 |
| MultiQC | — | not yet installed |

### 4. Assembly — planned

MEGAHIT, per sample. Contig headers prefixed with the sample id. Keep the assembly graph.

metaSPAdes: one attempt on the smallest sample (SRR6231137, Mallard Creek UP_A time
point 3, 15.05 M pairs). Record what happens.

Map each sample's reads back to its own assembly, and map every other sample of the same
system to it too. One BAM gives one coverage number per contig, so two organisms at
similar abundance and GC merge into one bin. With every system mate mapped, each contig
gets a coverage profile across samples, and contigs from the same genome rise and fall
together across timepoints and plant stages while contigs from different genomes do not.

QUAST. Nonpareil. Table: contigs, total bp, N50, assembled fraction per sample.

SCAPP on the graph for circular plasmids.

| Tools | Env | Status |
|---|---|---|
| MEGAHIT, metaSPAdes, read mapper for back-/cross-mapping, QUAST, Nonpareil, SCAPP | — | not yet installed |

### 5. Contig-level plasmids and AMR — planned

geNomad, end-to-end. Plasmid and virus contigs. CheckV on the viral contigs. MOB-suite
`mob_typer` on the plasmid contigs. AMRFinderPlus with `--plus` on the whole assembly.
RGI with CARD on the whole assembly. Tag each AMR hit as plasmid, virus, or chromosome by
its contig. CoverM contig for coverage and TPM.

PlasLLM for the plasmid work — source/version TBD, to confirm with user (see open items).

| Tools | Env | Status |
|---|---|---|
| geNomad + DB | `pd-genomad` | installed by stage 0 |
| CheckV + DB | `pd-checkv` | installed by stage 0 |
| MOB-suite + DB | `pd-mobsuite` | installed by stage 0 |
| AMRFinderPlus + DB | `pd-amrfinderplus` | installed by stage 0 |
| RGI + CARD | `pd-rgi` | installed by stage 0 |
| CoverM | — | not yet installed |
| PlasLLM | — | not yet installed; source unconfirmed |

### 6. Binning and MAG QC — planned

`jgi_summarize_bam_contig_depths` over all BAMs on each assembly. MetaBAT2. CheckM2.
MIMAG tiers: high, medium, low. Medium and up go on.

| Tools | Env | Status |
|---|---|---|
| MetaBAT2 (incl. `jgi_summarize_bam_contig_depths`) | — | not yet installed |
| CheckM2 + DB | `pd-checkm2` | installed by stage 0 |

### 7. MAG taxonomy — planned

dRep at 95 percent ANI. GTDB-Tk `classify_wf`, release 226. CoverM genome over all
samples for MAG abundance.

Table: quality, taxonomy, abundance per MAG.

| Tools | Env | Status |
|---|---|---|
| dRep, CoverM | — | not yet installed |
| GTDB-Tk 2.6.1 + GTDB r226 | `pd-gtdbtk` | installed by stage 0 |

### 8. MAG-level AMR and plasmids — planned

AMRFinderPlus per MAG. Which plasmid contigs sit in which MAG. Plasmid to host:
co-binning, MOB-suite host range, coverage correlation across samples. All three are
inferences.

| Tools | Env | Status |
|---|---|---|
| AMRFinderPlus, MOB-suite | `pd-amrfinderplus`, `pd-mobsuite` | installed by stage 0 |

### 9. Rollup — planned

Plasmid summary per sample. Resistome matrix with the mobile subset separated.

Comparison with the paper: RGI at read level (`bwt` mode) vs at contig level.

Depth experiment: subsample two samples with seqkit at 5, 10, 15, 20 M pairs and full,
assemble, geNomad, plot plasmid yield against depth.

Statistics: CLR, ALDEx2 or ANCOM-BC2, PERMANOVA as exploratory, Benjamini-Hochberg.
Comparisons within system and timepoint. Controls are a reference, never subtracted.

Antibiotic concentrations from the BioSample records as a covariate.

| Tools | Env | Status |
|---|---|---|
| RGI `bwt` + CARD | `pd-rgi` | installed by stage 0 (see open items re: `bwt` data) |
| seqkit; R with ALDEx2 / ANCOM-BC2 / vegan | — | not yet installed |

## Order of work

1. Setup, download all 92.
2. Pilot three samples through stages 2 to 8, with the metaSPAdes attempt at the same
   time:

   | Pilot id | Run | BioSample | Sample | Read pairs |
   |---|---|---|---|---|
   | SC_UPB_T4 | SRR6231181 | SAMN07839307 | Sugar Creek, Upstream B, time point 4 | 15,375,933 |
   | MC_INF_T4 | SRR6231191 | SAMN07839316 | Mallard Creek, Influent, time point 4 | 28,073,637 |
   | MNT_MNTA_T4 | SRR6231140 | SAMN07839331 | Mountain Control A, time point 4 | 29,495,353 |

3. Full run as job arrays.
4. Depth experiment.
5. Rollup.

## Status

- **Stage 0:** complete. Every step finished; `versions.txt` has no blank fields.
- **Stage 1:** complete. 92/92 runs PASS in `metadata/raw_manifest.tsv`: R1 pairs = R2
  pairs = SRA spots for every run, full-quality scores (6 distinct quality characters, HiSeq
  2500 binning), 184 FASTQ.gz read-only (252 GB), `.sra` temp files removed.
- Stages 2–9: not started.

## Known discrepancies / open items

1. **Tools not yet in `envs/` or `00_setup`:** FastQC, fastp, MultiQC,
   MEGAHIT, metaSPAdes, QUAST, Nonpareil, SCAPP, CoverM, MetaBAT2
   (`jgi_summarize_bam_contig_depths`), dRep, seqkit, R with ALDEx2 / ANCOM-BC2 / vegan,
   PlasLLM. Stage 4's back-/cross-mapping also needs a read mapper choice (not named in
   the spec).
2. **PlasLLM — source/version TBD, to confirm with user.** No tool by that name was found.
   Searched: web (`PlasLLM plasmid language model`, `"PlasLLM" plasmid`,
   `PlasLLM github metagenome plasmid detection large language model`); GitHub repository
   search API (`q=plasllm`: 0 results); conda-forge + bioconda (`mamba search *plasllm*`:
   no match); PyPI (`plasllm`: not found). Similarly named language-model plasmid tools
   that did turn up: PLASMe (transformer-based plasmid contig identification from
   short-read assemblies), PlasmidGPT (generative, trained on Addgene plasmids), PlasRAG
   (plasmid characterization/retrieval), PlasGO (plasmid protein GO prediction). Need the
   repo/paper link for the intended tool.
3. **"Same system" for cross-mapping (stage 4)** is clear for Mallard Creek (44 runs) and
   Sugar Creek (44 runs); the four controls (MNT_A/B, UW_A/B, time point 4 only) need a
   rule — e.g. map within each control pair, or against each other.
4. **RGI `bwt` mode (stage 9)** runs against CARD's `card_annotation` reference, which
   `00_setup` loads; CARD's WildCARD/variants data is not loaded. Decide whether the
   read-level comparison should include it.
5. **Antibiotic covariate:** the concentrations are present in the BioSample attributes
   (above), with `ND` for not-detected values; they are not yet pulled into
   `metadata/`. How to encode `ND` (zero, LOD/2, censored) is a modeling choice to make
   before stage 9.
