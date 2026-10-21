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
(`--sample_name_regex`, default `Sample_<Location>_<Site>_<TP>`); for a study where it
doesn't match, every run goes into one group, which only affects within-system
cross-mapping and the site statistics.

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
| `--metaspades` | `true`: one metaSPAdes attempt on the smallest run (or `params.assembly.metaspades_run`) |
| `--depth_experiment` | `false`: opt-in stage 9 depth experiment (`params.rollup.depth_*`) |

After a failure, fix the cause and resubmit the same command with `-resume`: finished
tasks are reused from the work dir (`/scratch/.../work_nf`).

## Outputs

`--outdir` mirrors the sbatch `results/` layout (`02_qc/`, `03_host/`, `04_assembly/<run>/`,
`05_contigs/<run>/`, `06_binning/<run>/`, `07_mags/`, `08_mag_links/`, `09_rollup/`,
`tables/`, plus `metadata/` with the fetched run table and sample sheet). Large files
(reads, contigs, BAMs) stay in the work dir. `pipeline_info/` has Nextflow's HTML report
(time, CPU, memory per task), timeline and trace.

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
  on the pilot (86/86 tasks, including the depth experiment).

## Validation against the sbatch pilot (2026-10-21)

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
