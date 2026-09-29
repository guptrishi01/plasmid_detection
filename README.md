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
├── config/config.yaml        # The one config: every path and parameter (params.<stage>)
├── envs/                     # Pinned conda env specs (<name>.yml); envs/lock/ = solved envs
├── scripts/
│   ├── lib/common.sh         # Sourced by every stage: conda, config -> env vars, helpers
│   ├── 00_setup.sbatch       # Stage 0: envs, human reference + index, databases
│   ├── 01*_*.sbatch          # Stage 1: SRA download + manifest
│   ├── 02*_ … 09*_*.sbatch   # Stages 2–9, one script per step (tables below)
│   └── R/                    # stats.R (stage 9 statistics), depth_plot.R
├── src/plasmid_detection/    # config.py (config loader), tables.py (per-stage tables)
├── tests/                    # pytest: config loader + table/tagging logic
├── metadata/                 # Sample sheet, NCBI SRA metadata, raw manifest, pilot list
├── versions.txt              # Every tool + database version installed (written by 00_setup)
├── requirements.txt, pyproject.toml
└── data/raw -> /scratch/rgupta25/plasmid_detection/raw   (gitignored)
```

Not tracked in git: `conda/` (envs + package cache), `db/` (reference + databases),
`data/`, `logs/`, sequence/index files, `*.pdf`. Stage outputs go to scratch under
`/scratch/rgupta25/plasmid_detection/{work,results}`; summary tables to `results/tables/`.

**`results/` in the repo** is a copy of the shareable part of the scratch results (reports,
logs, QC summaries, tables) for lab meetings, made with `bash scripts/sync_results.sh
[stage ...]`. It never copies reads, contigs, graphs, or BAMs, and skips files over 20 MB;
those stay on scratch.

**Storage split.** `paths.sw_root` = this repo in home (500G, not purged) holds `conda/`
and `db/`. `paths.data_root` = `/scratch/rgupta25/plasmid_detection` (5T, purgeable) holds
raw reads, intermediates, and results. **Config.** Every stage sources
`scripts/lib/common.sh`, which exports `config/config.yaml` as env vars
(`db.gtdb.dir` → `DB_GTDB_DIR`); stage parameters go under `params.<stage>`, never as
literals. **Envs.** Pinned specs in `envs/`; the solved env is exported to `envs/lock/`.
Details in [CLAUDE.md](CLAUDE.md).

## Workflow

Status legend: **implemented** = script exists and has produced its output;
**script written** = script exists, tools installed and pinned, not yet run on data.

### 0. Setup — implemented (`scripts/00_setup.sbatch`)

Git repository for all code. One config file for paths and parameters
(`config/config.yaml`). One conda environment per stage, versions pinned in a yml.
Human reference GRCh38.p14 (GCF_000001405.40) via NCBI Datasets CLI, with a Bowtie2
index. Databases: geNomad, CheckV, CheckM2, CARD, AMRFinderPlus, MOB-suite, GTDB r226.
Every version used is recorded in `versions.txt`.

<!-- STAGE0_VERSIONS -->
Installed for stages 0–1 (from `versions.txt`; stage 2–9 tools are listed per stage below):

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

### 2. Read QC — script written

FastQC on raw reads. fastp. FastQC on trimmed reads. MultiQC over everything.
Table: raw pairs, trimmed pairs, Q30, duplication, insert size, per sample.

| Script | What it does | Env |
|---|---|---|
| `02_read_qc` | per run: FastQC raw → fastp → FastQC trimmed | pd-qc |
| `02b_qc_summary` | MultiQC over all three + `tables/02_read_qc.tsv` | pd-qc |

fastp options (`params.qc.fastp_args`) mirror the paper's Trimmomatic intent: adapter
detection for PE, sliding window Q20, minimum length 50. No `--dedup`: duplication is
reported, not removed. Tools: FastQC 0.12.1, fastp 1.3.7, MultiQC 1.35.

### 3. Host depletion — script written

Bowtie2 against the human index. Keep pairs where neither mate maps.
Table: raw, trimmed, human, retained pairs per sample. MultiQC of the Bowtie2 logs.

| Script | What it does | Env |
|---|---|---|
| `03_host_depletion` | per run: Bowtie2 (`--very-sensitive`) vs GRCh38.p14; keep pairs with SAM flags 4 and 8 both set | pd-host |
| `03b_host_summary` | MultiQC of the Bowtie2 logs + `tables/03_host_depletion.tsv` | pd-qc |

"Neither mate maps" is `samtools view -f 12`, not Bowtie2's `--un-conc` (which also keeps
pairs where only one mate maps). The table step fails if Bowtie2's input pair count
differs from fastp's output, so stage 3 can't silently run on stale reads.

### 4. Assembly — script written

MEGAHIT, per sample. Contig headers prefixed with the sample id. Keep the assembly graph.
metaSPAdes: one attempt on the smallest sample. Record what happens.

Map each sample's reads back to its own assembly, and map every other sample of the same
system to it too. One BAM gives one coverage number per contig, so two organisms at
similar abundance and GC merge into one bin. With every system mate mapped, each contig
gets a coverage profile across samples, and contigs from the same genome rise and fall
together across timepoints and plant stages while contigs from different genomes do not.

QUAST. Nonpareil. Table: contigs, total bp, N50, assembled fraction per sample.
SCAPP on the graph for circular plasmids.

| Script | What it does | Env |
|---|---|---|
| `04a_megahit` | per run: MEGAHIT (`--presets meta-large`, contigs ≥ 1000 bp); headers `>SRR…_k127_n`; FASTG graph from the largest-k contig set | pd-assembly |
| `04b_metaspades` | one metaSPAdes attempt; outcome (exit code, hours, peak RAM, contigs, last error) → `metaspades/attempt.tsv` whether it succeeds or not | pd-assembly |
| `04c_map` | per assembly: Bowtie2 of its own reads + every system mate's reads → sorted, indexed BAMs | pd-assembly |
| `04d_assembly_qc` | per run: QUAST + Nonpareil (host-depleted R1) | pd-assembly |
| `04e_scapp` | per run: SCAPP on the k-max FASTG with the run's reads | pd-scapp |
| `04f_assembly_summary` | `tables/04_assembly.tsv`; assembled fraction = own reads' Bowtie2 overall alignment rate to own assembly | pd-core |

Mapper: **Bowtie2** (chosen for stage 4; the spec doesn't name one). "Same system" =
same `location` in the sample sheet among `paths.run_samples` (Sugar Creek, Mallard Creek,
Mountain_Control, Uwharrie_Control, so the controls cross-map only within their pair).
The graph is `intermediate_contigs/k<K>.contigs.fa` for the largest K MEGAHIT reached,
converted with `megahit_toolkit contig2fastg` (`final.contigs.fa` merges contigs from all
k, so it isn't a single graph). Tools: MEGAHIT 1.2.9, SPAdes 4.3.0, Bowtie2 2.5.5,
samtools 1.24, QUAST 5.3.0, Nonpareil 3.5.5, SCAPP 0.1.4 (own env: it pins Python 3.7).

### 5. Contig-level plasmids and AMR — script written (PlasLLM pending)

geNomad, end-to-end. Plasmid and virus contigs. CheckV on the viral contigs. MOB-suite
`mob_typer` on the plasmid contigs. AMRFinderPlus with `--plus` on the whole assembly.
RGI with CARD on the whole assembly. Tag each AMR hit as plasmid, virus, or chromosome by
its contig. CoverM contig for coverage and TPM.

| Script | What it does | Env |
|---|---|---|
| `05a_genomad` | per run: geNomad end-to-end on all contigs | pd-genomad |
| `05b_checkv_mobtyper` | per run: CheckV on virus contigs; `mob_typer --multi` on plasmid contigs | pd-checkv, pd-mobsuite |
| `05c_amr` | per run: AMRFinderPlus `--plus` and RGI (contig mode, DIAMOND) on the whole assembly | pd-amrfinderplus, pd-rgi |
| `05d_coverage_tag` | per run: CoverM contig (mean, TPM, read count) over all stage-4 BAMs; tag AMR hits; reads per hit per BAM (`samtools bedcov -c`) | pd-coverm, pd-core, pd-assembly |

Tagging: a hit is **virus** if its contig is a geNomad virus or the hit overlaps a geNomad
provirus region, **plasmid** if its contig is a geNomad plasmid, otherwise
**chromosome** (= not called plasmid or virus). Reads over each hit's own coordinates are
the integer counts behind the stage 9 resistome and statistics. **PlasLLM** belongs next
to geNomad here; it isn't installed until its source is confirmed (open items). Tools:
geNomad 1.12.0, CheckV 1.1.1, MOB-suite 3.1.9, AMRFinderPlus 4.2.7, RGI 6.0.8, CoverM 0.8.0.

### 6. Binning and MAG QC — script written

`jgi_summarize_bam_contig_depths` over all BAMs on each assembly. MetaBAT2. CheckM2.
MIMAG tiers: high, medium, low. Medium and up go on.

| Script | What it does | Env |
|---|---|---|
| `06a_binning` | per assembly: depths over all its BAMs → MetaBAT2 (min contig 1500, fixed seed) | pd-binning |
| `06b_checkm2_mimag` | per run: CheckM2 → MIMAG tier per bin → medium+ copied to `binning/pass/` | pd-checkm2, pd-core |

Tiers from CheckM2 completeness/contamination (`params.binning.mimag_*`): high ≥ 90 / < 5,
medium ≥ 50 / < 10, else low. MIMAG's high-quality tier also requires rRNAs and ≥ 18 tRNAs,
which CheckM2 doesn't report (open items). Tools: MetaBAT2 2.18, CheckM2 1.1.0.

### 7. MAG taxonomy — script written

dRep at 95 percent ANI. GTDB-Tk `classify_wf`, release 226. CoverM genome over all
samples for MAG abundance. Table: quality, taxonomy, abundance per MAG.

| Script | What it does | Env |
|---|---|---|
| `07a_drep` | all medium+ MAGs → dRep `-sa 0.95`, CheckM2 values via `--genomeInfo` | pd-drep |
| `07b_gtdbtk` | GTDB-Tk classify_wf on the dereplicated MAGs (r226) | pd-gtdbtk |
| `07c_mag_abundance` | CoverM genome (all runs) + `tables/07_mags.tsv` | pd-coverm, pd-core |

dRep's own quality filter is set to the MIMAG medium thresholds (its defaults, 75/25,
would silently drop MAGs stage 6 passed). CoverM genome maps reads itself (minimap2-sr);
Bowtie2 is the stage-4 mapper only. Tools: dRep 3.7.1, GTDB-Tk 2.6.1, CoverM 0.8.0.

### 8. MAG-level AMR and plasmids — script written

AMRFinderPlus per MAG. Which plasmid contigs sit in which MAG. Plasmid to host:
co-binning, MOB-suite host range, coverage correlation across samples. All three are
inferences.

| Script | What it does | Env |
|---|---|---|
| `08a_mag_amr` | AMRFinderPlus `--plus` per dereplicated MAG → `tables/08_mag_amr.tsv` | pd-amrfinderplus |
| `08b_plasmid_hosts` | `tables/08_plasmids_in_mags.tsv` + `08_plasmid_host_inference.tsv` (the three inferences side by side) | pd-core |

Coverage correlation = Spearman between a plasmid contig's coverage profile and each
medium+ bin's mean profile across the BAMs on that assembly; left blank below
`params.mag_links.min_samples_for_correlation` (5) samples, so with the pilot (2 BAMs per
assembly) it is always blank, by design.

### 9. Rollup — script written

Plasmid summary per sample. Resistome matrix with the mobile subset separated.
Comparison with the paper: RGI at read level (`bwt` mode) vs at contig level.
Depth experiment: subsample two samples with seqkit at 5, 10, 15, 20 M pairs and full,
assemble, geNomad, plot plasmid yield against depth.
Statistics: CLR, ALDEx2 or ANCOM-BC2, PERMANOVA as exploratory, Benjamini-Hochberg.
Comparisons within system and timepoint. Controls are a reference, never subtracted.
Antibiotic concentrations from the BioSample records as a covariate.

| Script | What it does | Env |
|---|---|---|
| `09a_rgi_bwt` | per run: RGI bwt (Bowtie2) on host-depleted reads | pd-rgi |
| `09b_depth_experiment` | per (run, level): seqkit subsample → MEGAHIT → geNomad → yield | pd-assembly, pd-genomad |
| `09c_depth_summary` | adds the "full" point (stage 4/5 outputs), table + plot | pd-core, pd-stats |
| `09d_rollup` | plasmid summary; resistome matrices (all / mobile / chromosome, for AMRFinderPlus and RGI); RGI read vs contig | pd-core |
| `09e_stats` | BioSample antibiotics → `metadata/stats_samples.tsv`; CLR, ALDEx2/ANCOM-BC2, PERMANOVA, BH | pd-core, pd-stats |

Resistome values are read counts over each AMR hit (own reads on own assembly). Depth
levels above a run's host-depleted pair count are skipped and logged (SC_UPB_T4 has
15.4 M raw pairs, so its 15 M and 20 M points won't exist). Statistics skip any test with
fewer than `params.stats.min_group_size` (3) runs per group, so the 4-run pilot skips all
of them: it exercises the code path without producing results. Both the ALDEx2 and the
ANCOM-BC2 paths were checked end to end on synthetic counts. Tools: seqkit 2.14.0,
R 4.5.3, ALDEx2 1.42.0, ANCOMBC 2.14.0, vegan 2.7-5, ggplot2 4.0.3.

## Running the pilot, one step at a time

Runs: `metadata/pilot_samples.txt` = SC_UPB_T4 (SRR6231181), SC_INF_T4 (SRR6231184),
MC_UPB_T4 (SRR6231217), MC_INF_T4 (SRR6231191): two per system at time point 4, so
within-system cross-mapping is exercised in both systems. Every script
reads `paths.run_samples`; for the full run, point it at
`metadata/srp121672_wgs_run_accessions.txt` and use `--array=1-92`. Per-run scripts need
`--array=1-N` (N = runs in the list); the rest are single jobs. Check each step's log in
`logs/` and its table in `results/tables/` before the next.

```bash
cd ~/gibas/plasmid_detection
N=$(grep -c . metadata/pilot_samples.txt)
sbatch --array=1-$N scripts/02_read_qc.sbatch;        sbatch scripts/02b_qc_summary.sbatch
sbatch --array=1-$N scripts/03_host_depletion.sbatch; sbatch scripts/03b_host_summary.sbatch
sbatch --array=1-$N scripts/04a_megahit.sbatch
sbatch scripts/04b_metaspades.sbatch                   # alongside the rest of stage 4
sbatch --array=1-$N scripts/04c_map.sbatch             # needs 04a done for every run
sbatch --array=1-$N scripts/04d_assembly_qc.sbatch
sbatch --array=1-$N scripts/04e_scapp.sbatch
sbatch scripts/04f_assembly_summary.sbatch
sbatch --array=1-$N scripts/05a_genomad.sbatch
sbatch --array=1-$N scripts/05b_checkv_mobtyper.sbatch
sbatch --array=1-$N scripts/05c_amr.sbatch
sbatch --array=1-$N scripts/05d_coverage_tag.sbatch
sbatch --array=1-$N scripts/06a_binning.sbatch
sbatch --array=1-$N scripts/06b_checkm2_mimag.sbatch
sbatch scripts/07a_drep.sbatch
sbatch scripts/07b_gtdbtk.sbatch
sbatch scripts/07c_mag_abundance.sbatch
sbatch scripts/08a_mag_amr.sbatch
sbatch scripts/08b_plasmid_hosts.sbatch
sbatch --array=1-$N scripts/09a_rgi_bwt.sbatch
sbatch --array=1-8 scripts/09b_depth_experiment.sbatch   # 2 runs x 4 levels
sbatch scripts/09c_depth_summary.sbatch
sbatch scripts/09d_rollup.sbatch
sbatch scripts/09e_stats.sbatch
```

Submit each line after the previous step finishes (or chain with
`--dependency=afterok:<jobid>`). A script whose input from an earlier step is missing
stops with "missing input: … (run the earlier stage first)". Everything tunable is in
`config/config.yaml` under `params.<stage>`; tuning never needs a script edit.

## Order of work

1. Setup, download all 92.
2. Pilot three samples through stages 2 to 8, with the metaSPAdes attempt at the same
   time. **Current pilot (decided 2026-09-28): upstream B + influent at time point 4 in
   both systems, SC_UPB_T4 (SRR6231181), SC_INF_T4 (SRR6231184), MC_UPB_T4 (SRR6231217),
   MC_INF_T4 (SRR6231191), through stages 2–9.** The spec's original three:

   | Pilot id | Run | BioSample | Sample | Read pairs |
   |---|---|---|---|---|
   | SC_UPB_T4 | SRR6231181 | SAMN07839307 | Sugar Creek, Upstream B, time point 4 | 15,375,933 |
   | MC_INF_T4 | SRR6231191 | SAMN07839316 | Mallard Creek, Influent, time point 4 | 28,073,637 |
   | MNT_MNTA_T4 | SRR6231140 | SAMN07839331 | Mountain Control A, time point 4 | 29,495,353 |

3. Full run as job arrays.
4. Depth experiment.
5. Rollup.

## Status

- **Stage 0:** complete. `versions.txt` has no blank fields (every env, including stages 2–9).
- **Stage 1:** complete. 92/92 runs PASS in `metadata/raw_manifest.tsv`: R1 pairs = R2
  pairs = SRA spots for every run, full-quality scores (6 distinct quality characters, HiSeq
  2500 binning), 184 FASTQ.gz read-only (252 GB), `.sra` temp files removed.
- **Stages 2–9:** scripts written; envs built and pinned; tool flags checked against each
  installed tool's `--help` or source; table and tagging logic unit-tested; R statistics
  and plot checked on synthetic data. **Not yet run on real data**: the pilot is run step
  by step by hand.

## Known discrepancies / open items

1. **PlasLLM: source/version TBD.** No tool by that name was found (web, GitHub API,
   conda-forge/bioconda, PyPI). Similarly named: PLASMe, PlasmidGPT, PlasRAG, PlasGO. Needs
   the repo/paper link; it goes in stage 5 next to geNomad.
2. **metaSPAdes sample.** The spec says the smallest sample: of all 92 that is SRR6231137
   (MC UP_A T3, 15.05 M pairs), which isn't in the pilot. `04b` defaults to the smallest
   pilot run (SRR6231181); set `params.assembly.metaspades_run: SRR6231137` (after running
   stages 2–3 on it) to follow the spec literally.
3. **MIMAG "high"** also requires 23S/16S/5S rRNA and ≥ 18 tRNAs; CheckM2 reports only
   completeness/contamination, so "high" here is that part only. barrnap + tRNAscan-SE
   would make it the full criterion.
4. **Statistics design.** `params.stats.comparisons` (upstream vs downstream creek;
   influent vs plant) and `params.stats.site_groups` are a proposal for "comparisons
   within system and timepoint". Confirm before the full run.
5. **Antibiotic covariate encoding.** BioSample values are numbers, `ND` (not detected)
   or `<LOQ` (below quantification). `params.stats.nd_encoding` / `loq_encoding` choose
   na | zero | half_min (default na); the raw strings are kept alongside.
6. **RGI `bwt` (stage 9)** runs against CARD canonical references; WildCARD/variants
   data is not loaded. Decide whether the read-level comparison should include it.
7. **Control cross-mapping (stage 4)** defaults to within-pair (MNT_A↔MNT_B,
   UW_A↔UW_B) via the sample sheet's `location`. Only matters for the full run.
