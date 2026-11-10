# Nextflow version

The same pipeline as `scripts/01`–`09` (stages 1–9), run by Nextflow instead of by hand.
Nextflow submits one Slurm job per sample per step, starts each step as soon as its inputs
exist, runs samples in parallel (at most 20 jobs queued at once), and resumes after a
failure without redoing finished work. The sbatch scripts stay as the reference version.

## Input: an NCBI identifier

The first process asks NCBI E-utilities for the run table of `--id` and keeps one library
strategy (`--library_strategy`, default `WGS`):

| `--id` | Runs |
|---|---|
| `PRJNA415974` (default) | the 92 WGS runs of Lambirth et al. 2018 (Gibas lab); its 496 16S runs are dropped |
| `SRP121672` | the same study, by SRA study accession |
| `SRR6231181,SRR6231191` | just those runs |
| any other BioProject / SRA study | its WGS runs |

`--runs` (comma list or file) narrows the set; the `pilot` profile sets it to
`metadata/pilot_samples.txt`. Location, site and timepoint come from the SRA SampleName
(`--sample_name_regex`, default `Sample_<Location>_<Site>_<TP>`). Each run's binning group
is the `params.binning.groups` entry listing its site (`config/config.yaml`, the same 11
groups as the sbatch version); a parsed site in no group stops the run. For a study where
the regex doesn't match, every run goes into one group, `all`.

Reads are downloaded with sra-tools only (prefetch → vdb-validate → fasterq-dump → pigz,
full-quality `.sra`, pair counts checked against SRA spots), into `--raw_dir`
(`/scratch/.../raw`) as a permanent cache: runs already there are never fetched again.

## Running

Submit the head job (2 CPUs, up to 7 days); it submits everything else:

```bash
cd ~/gibas/plasmid_detection
sbatch nextflow/submit.sbatch -profile slurm,pilot -resume     # 4-run pilot
sbatch nextflow/submit.sbatch -profile slurm -resume           # all 92 WGS runs (default id)
sbatch nextflow/submit.sbatch -profile slurm -resume --id PRJNAxxxxxx
```

`nextflow/run.sh` (called by `submit.sbatch`) turns `config/config.yaml` into the
Nextflow params-file (`python -m plasmid_detection.config nextflow`), so every tool
parameter, conda env and database path still comes from that one file. Run-level options
live in `nextflow/nextflow.config`:

| Option | Default |
|---|---|
| `--outdir` | `/scratch/rgupta25/plasmid_detection/results_nf` (`results_nf_pilot` with `pilot`) |
| `--metaspades` | `true`: one metaSPAdes attempt on the smallest run (or `params.assembly.metaspades_run`); `--metaspades false` skips it (the 92-run run does) |
| `--depth_experiment` | `false`: opt-in stage 9 depth experiment (`params.rollup.depth_*`) |

After a failure, fix the cause and resubmit the same command with `-resume`: finished
tasks are reused from the work dir (`/scratch/.../work_nf`).

## Outputs

`--outdir` mirrors the sbatch `results/` layout (`02_qc/`, `03_host/`, `04_assembly/<run>/`,
`05_contigs/<run>/`, `06_binning/<run>/`, `07_mags/`, `08_mag_links/`, `09_rollup/`,
`tables/`, plus `metadata/` with the fetched run table and sample sheet). Large files
(reads, contigs, BAMs) stay in the work dir. `pipeline_info/` has Nextflow's HTML report
(time, CPU, memory per task), timeline and trace.

## Binning groups and dRep (2026-11-03 design)

Same as the sbatch version: each assembly is mapped with the reads of every run in its
binning group, dRep runs per group and then across groups (the catalog), and GTDB-Tk,
CoverM genome and stage 8 use the catalog. Two Nextflow-specific choices:

- `MAP_GROUP` does everything that reads the cross-sample BAMs in one task (Bowtie2 of
  the group's reads, `jgi_summarize_bam_contig_depths`, CoverM contig, `samtools bedcov`
  over the AMR hits) and keeps only the run's own BAM. A later task can't delete another
  task's outputs without breaking `-resume`, so this replaces the sbatch `06c` cleanup.
  It waits for `TAG_AMR` (the hits' BED), which is the tagging half of the sbatch `05d`.
- `DREP_GROUP` (one task per group with medium+ MAGs) → `DREP_CATALOG`; `07_mags/groups/
  <group>/drep/data_tables/` and `07_mags/drep/data_tables/` as in the sbatch layout.

| sbatch | Nextflow |
|---|---|
| `04c_map` + `05d` (CoverM, bedcov) + `06a` (depth) + `06c` | `MAP_GROUP` |
| `05d` (tagging) | `TAG_AMR` |
| `06a` (MetaBAT2) | `BINNING` |
| `07a` (two levels) | `GENOME_INFO` → `DREP_GROUP` → `DREP_CATALOG` |

## Notes

- Shared filesystem assumed: summary steps rebuild the per-run layout `tables.py` expects
  with symlinks to the per-run task outputs (true on the cluster; all work dirs on scratch).
- Envs and databases are the ones built by `scripts/00_setup.sbatch` (`conda/envs/pd-*`,
  `db/`); Nextflow itself runs from `pd-nextflow` (Nextflow 26.10.2). Setup stays outside
  the pipeline.
- Not included: the PLASMe comparison (`05e`/`05f`) and the under-sequencing side analysis
  (`scripts/side/`). Both were supplementary analyses for this study; geNomad is the
  plasmid caller.
- Checked with `nextflow lint` (no errors or warnings) and a `-stub` run of the full DAG
  on the pilot (89/89 tasks, including the depth experiment; group binning: 4 `MAP_GROUP`
  tasks with the right group runs, 2 `DREP_GROUP`, 1 `DREP_CATALOG`).
- Publishing: a `publishDir` pattern that points inside an emitted directory (e.g.
  `gtdbtk/*.summary.tsv` when `gtdbtk/` is also emitted) publishes nothing, so the run of
  2026-10-21 never wrote `07_mags/gtdbtk/` or the dRep tables to `--outdir`. GTDB-Tk and
  dRep now copy those files out as their own outputs.

## Validation of group binning (2026-11-08)

Pilot run (head job 27707762, `-resume 8001ba9b-…` = the 2026-10-20 session: 53 tasks
cached, 27 run, 0 failed, 70 min) against the sbatch group-binning results:

| Check | Result |
|---|---|
| Mapping sets | CoverM and `depth.txt` columns = the run's group mates (upstream 181+217, influent 184+191); Bowtie2 log per mate; counted = tagged AMR hits |
| Stage 7 layout | each group's `Cdb.csv` = its runs' medium+ MAGs; catalog = the group MAGs' representatives; GTDB-Tk summary published (57/57 bacteria); one AMRFinderPlus file per catalog MAG |
| Medium+ MAGs → per group → catalog | 73 → upstream 30 → 26, influent 43 → 31 → 57 (sbatch: identical counts) |
| Per run high/medium | 2/11, 1/16, 2/15, 4/22 (sbatch 2/11, 1/16, 1/15, 3/24) |
| Reads in MAGs | 17.1 / 19.3 / 33.3 / 33.2% (sbatch 17.8 / 20.0 / 32.1 / 33.5%) |
| Plasmid contigs in a medium+ MAG | 549 (sbatch 487) |
| AMR-type hits in MAGs | 13 in 8 MAGs (sbatch 10 in 5); all 10 sbatch genes present, same hosts (Aeromonas caviae, Phocaeicola vulgatus, Bacteroides uniformis, Acinetobacter johnsonii, Moraxella_A); extra: erm(B) in Acinetobacter sp., tet(Q) in Streptococcus suis, tet(32) in P. vulgatus |

The differences come from the assemblies, as in the run below: these MEGAHIT assemblies are
the 2026-10-20 ones, not the sbatch ones, so bins and MAGs differ slightly.

Launch notes: group inputs are sorted by run before `MAP_GROUP` and `DREP_GROUP` (groupTuple keeps arrival order, which is part of the cache key: a resume of the 92-run run re-ran every `MAP_GROUP` until this was fixed, 2026-11-10). Keep the binning group out of `meta` (it is part of every task's cache key;
adding it re-runs everything), and resume a specific session (`-resume <session id>` from
`nextflow log`) when other runs, e.g. a `-stub` test, came after the one to resume.

## Validation against the sbatch pilot (2026-10-21, creek-pair design)

Pilot run (job 27394761, 0 failed) compared with the sbatch results on scratch:

| Stage | Result |
|---|---|
| 2–3 | `02_read_qc.tsv`, `03_host_depletion.tsv` identical |
| 4 | Same contig counts and N50 (±1 bp); 99.97–99.98% of contig sequences identical, but MEGAHIT numbers contigs differently between runs (`k127_<n>`) and 6–27 contigs per run differ |
| 4 metaSPAdes | Identical (exit 0, 38,657 contigs ≥ 1 kb) |
| 5 | geNomad plasmid/virus counts and AMR tags (source × plasmid/virus/chromosome) identical per run; SCAPP 3→2 (181), 47→46 (184) |
| 6–8 | Bins and MAGs differ slightly (45 → 43 dereplicated MAGs; 7 high-quality both runs) |
| 9 | Resistome matrices: same genes (95 AMRFinderPlus, 35 mobile, 127 RGI); statistics skipped as designed |

The stage 4–8 differences come from MEGAHIT, which is not fully reproducible with
multiple threads (same command, memory and 32 threads in both runs): contig order and
names change, and MetaBAT2 bins depend on contig order. Not a pipeline difference.
