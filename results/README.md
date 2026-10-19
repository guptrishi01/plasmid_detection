# Pilot results (stages 2–7) and side analyses

Four runs at time point 4: Sugar Creek (SC) and Mallard Creek (MC), upstream B (UP) and
influent (INF). Every number below comes from the files in this folder (tables in `tables/`, per-run
outputs in `<stage>/<run>/`); reads, contigs, graphs and BAMs stay on scratch.

## Stages 2–3: read QC and host depletion

(`tables/02_read_qc.tsv`, `tables/03_host_depletion.tsv`)

| Run | Sample | Raw pairs | Trimmed pairs | Q30 raw → trimmed | Human pairs | Retained pairs (% raw) |
|---|---|---|---|---|---|---|
| SRR6231181 | SC UP_B | 15,375,933 | 14,673,925 | 96.0 → 98.1% | 5,093 | 14,668,832 (95.4%) |
| SRR6231184 | SC INF | 21,031,655 | 20,237,685 | 95.7 → 97.8% | 3,076 | 20,234,609 (96.2%) |
| SRR6231217 | MC UP_B | 26,183,846 | 24,936,856 | 94.9 → 97.6% | 7,827 | 24,929,029 (95.2%) |
| SRR6231191 | MC INF | 28,073,637 | 26,815,410 | 96.0 → 98.1% | 5,397 | 26,810,013 (95.5%) |

Duplication 0.4–0.7%. FastQC on the host-depleted reads (`03_host/multiqc/hostfree_qc.html`)
matches the trimmed reads: same GC, length and duplication; only GC-content and
length-distribution flags (expected for a mixed community after trimming).

## Stage 4: assembly

(`tables/04_assembly.tsv`, `04_assembly/<run>/quast/`, `nonpareil/`)

| Run | Contigs ≥ 1 kb | Contigs ≥ 10 kb | Largest | Total | N50 | Assembled | Nonpareil coverage | SCAPP circular plasmids |
|---|---|---|---|---|---|---|---|---|
| SRR6231181 SC UP_B | 47,604 | 496 | 164 kb | 95.2 Mb | 1,982 | 25.4% | 34% | 3 |
| SRR6231184 SC INF | 71,052 | 740 | 80 kb | 146.5 Mb | 2,134 | 46.6% | 59% | 47 |
| SRR6231217 MC UP_B | 81,609 | 918 | 246 kb | 166.2 Mb | 2,065 | 35.2% | 46% | 12 |
| SRR6231191 MC INF | 92,591 | 1,210 | 117 kb | 198.1 Mb | 2,260 | 49.6% | 61% | 47 |

Assembled = own reads mapped back to own assembly; cross-sample mapping 0.4–0.8%.
Nonpareil coverage is the modelled value (R `Nonpareil`, `C`); reaching 95% would need
144–307 Gbp per run, so rare taxa are under-sampled, more so upstream (higher diversity).
metaSPAdes on SRR6231181: exit 0, 1.25 h, 48.7 GB peak RAM, 38,657 contigs ≥ 1 kb
(MEGAHIT: 47,604), so MEGAHIT is used for all runs.

## Stage 5: plasmids and AMR on contigs

(`05_contigs/<run>/`)

| Run | geNomad plasmid contigs | Mobilizable / conjugative (MOB-typer) | geNomad virus contigs | AMRFinderPlus AMR-type hits | … on plasmid contigs | RGI hits (Perfect + Strict) |
|---|---|---|---|---|---|---|
| SRR6231181 SC UP_B | 124 | 0 / 0 | 648 | 0 | 0 | 24 |
| SRR6231184 SC INF | 1,863 | 111 / 1 | 691 | 80 | 22 | 181 |
| SRR6231217 MC UP_B | 159 | 0 / 0 | 2,441 | 0 | 0 | 51 |
| SRR6231191 MC INF | 2,067 | 94 / 4 | 1,542 | 82 | 30 | 169 |

Influent carries 13–15× more plasmid contigs than upstream, and every AMRFinderPlus
AMR-type hit on a plasmid is in influent: tetracycline (tet(A), tet(C), tet(O), tet(X)),
sulfonamide (sul1, sul2), quinolone (qnrS2, qnrB19), macrolide (mph(A), msr(E)),
β-lactam (blaMCA). `--plus` adds stress (metal, biocide) and virulence hits (about half
of all AMRFinderPlus hits), so resistome analyses filter `Type == AMR`. geNomad counts
plasmid contigs (fragments), SCAPP counts complete circular plasmids. Most plasmid
contigs type as non-mobilizable because fragments often lack the transfer genes.

### Plasmid callers: geNomad vs PLASMe

(`tables/05_plasmid_callers.tsv`, `05_contigs/<run>/plasme*/`)

**geNomad is the main caller** (decided 2026-10-18; the goal is comparing plasmids across
locations). PLASMe was run in its default (balance) and high-precision modes, each with and
without dropping calls that are ≥ 50% chromosome-shared (PLASMe `amb_region`).
Cells: contigs called · also called by geNomad · contigs with AMR-type genes · geNomad's
median plasmid score of the calls.

| Caller | 181 SC UP | 217 MC UP | 184 SC INF | 191 MC INF |
|---|---|---|---|---|
| geNomad | 124 · 124 · 0 · 0.82 | 159 · 159 · 0 · 0.84 | 1,863 · 1,863 · 15 · 0.97 | 2,067 · 2,067 · 23 · 0.97 |
| PLASMe balance | 154 · 7 · 0 · 0.20 | 317 · 9 · 0 · 0.24 | 5,395 · 904 · 46 · 0.58 | 4,503 · 823 · 46 · 0.61 |
| chromosome filter only | 126 · 7 · 0 · 0.20 | 249 · 7 · 0 · 0.25 | 3,413 · 748 · 29 · 0.71 | 3,160 · 705 · 31 · 0.71 |
| high precision only | 74 · 6 · 0 · 0.25 | 171 · 7 · 0 · 0.27 | 4,886 · 857 · 46 · 0.60 | 3,882 · 767 · 41 · 0.65 |
| high precision + filter | 56 · 6 · 0 · 0.25 | 131 · 5 · 0 · 0.26 | 3,035 · 708 · 29 · 0.75 | 2,674 · 657 · 27 · 0.77 |

In upstream creek water, under 11% of PLASMe's calls are geNomad calls in every mode, and
geNomad scores them ~0.25 (chromosome-like); in influent ~20% overlap. The chromosome filter
is the fix that matters; high precision adds a little. A caller whose false-positive rate
differs by habitat would bias the cross-location comparison, hence geNomad.

## Stage 6: binning and MAG quality

(`06_binning/<run>/mimag.tsv`, `checkm2/quality_report.tsv`)

| Run | Bins (MetaBAT2) | High | Medium | Low | Kept (medium+) |
|---|---|---|---|---|---|
| SRR6231181 SC UP_B | 28 | 2 | 6 | 20 | 8 |
| SRR6231184 SC INF | 41 | 2 | 15 | 24 | 17 |
| SRR6231217 MC UP_B | 38 | 1 | 9 | 28 | 10 |
| SRR6231191 MC INF | 50 | 4 | 20 | 26 | 24 |

High = ≥ 90% complete and < 5% contaminated, medium = ≥ 50% / < 10% (CheckM2; MIMAG's rRNA
and tRNA requirement for "high" is not checked). 59 medium-or-better MAGs; dRep (95% ANI)
reduces them to 45, as 14 are the same species in both runs of a creek pair. With only two
BAMs per assembly in the pilot, coverage-based binning is weaker than it will be in the full run.

## Stage 7: MAG taxonomy and abundance

(`tables/07_mags.tsv`, `07_mags/drep/`, `07_mags/gtdbtk/`, `07_mags/coverm_genome.tsv`)

45 dereplicated MAGs (7 high, 38 medium), all classified by GTDB-Tk 2.6.1 (GTDB r226),
all bacterial: Gammaproteobacteria 16, Actinomycetes 13, Bacteroidia 6, Clostridia 5,
Campylobacteria 2, Bacilli 2, Fusobacteriota 1. Influent MAGs include gut taxa (Prevotella
copri, Agathobacter rectalis, Bifidobacterium adolescentis, Alistipes putredinis).

| Run | Reads assigned to the 45 MAGs (CoverM relative abundance) |
|---|---|
| SRR6231181 SC UP_B | 15.5% |
| SRR6231184 SC INF | 29.5% |
| SRR6231217 MC UP_B | 16.2% |
| SRR6231191 MC INF | 31.9% |

The MAGs capture about twice as much of the influent community as of the creek water,
consistent with the Nonpareil coverage and assembled fractions.

## Audits

 Stage 2: all FastQC/fastp reports present, trimmed = fastp output. Stage 3:
retained = Bowtie2 input − human, R1 = R2, read ids paired. Stage 4: 18/18 jobs exit 0,
contig FASTA = QUAST counts and lengths, mapping input = retained pairs, 8/8 BAMs pass
`samtools quickcheck`. Stage 5: 16/16 jobs exit 0, every RGI contig is in the assembly,
tagged hits = counted hits. RGI logs "Requested rname SRR…_k127 does not exist" a few
times per influent run: its ORF-extension step cuts names at underscores; the hits are
still reported with full contig names. Stage 6: 8/8 jobs exit 0, CheckM2 rows = bins,
`pass/` = 59 = high + medium; dRep kept all 59 through its filters (`Cdb.csv` = 59 rows).
Stage 7: dRep and GTDB-Tk exit 0; 45/45 MAGs classified; `07_mags.tsv` has 45 rows. 07c
exited 13 after writing its table (SIGPIPE from `column | head` in the display line; fixed).

## Side analysis: odds of under-sequencing (all 92 runs)

(`side_undersequencing/`; scripts `scripts/side/us01–us02`, `scripts/R/undersequencing.R`)

Nonpareil coverage C (share of the community the reads sampled) for all 92 runs, on raw R1
(only the pilot has been trimmed; on the 4 pilot runs raw and host-depleted C differ by
≤ 0.013). Under-sequenced = C < 0.5 (`params.undersequencing.coverage_threshold`; refit at
0.4–0.7). Predictors: log2 sequencing effort, site group, creek system.

| Site group | Runs | Under-sequenced | Mean C |
|---|---|---|---|
| creek upstream | 16 | 10 | 0.43 |
| remote controls | 4 | 2 | 0.39 |
| plant | 32 | 3 | 0.61 |
| creek downstream | 16 | 2 | 0.59 |
| influent | 8 | 0 | 0.61 |
| sewer | 16 | 0 | 0.60 |

| Model | Upstream vs plant | Per doubling of effort | LOOCV AUC (C < 0.5) |
|---|---|---|---|
| Logistic | not estimable (separation: no influent/sewer run is under-sequenced) | — | 0.725 |
| Firth logistic | odds of under-sequencing ×12.3 (3.1–60) | ×0.90 (0.08–10) | 0.712 |
| Beta regression on C | odds of coverage ×0.50 (0.39–0.65) | ×1.45 (1.03–2.06) | **0.755** |

Upstream creek water is the habitat most likely to be under-sequenced at these depths;
influent and sewer never are. Beta regression uses C as a continuous value, ranks runs best
and is the only model that detects the effect of depth (which varies only ~2-fold, 1.9–3.9
Gbp R1). 17 events limit the logistic models (≈ 2 per term). On the pilot, C ranks the runs
in the same order as the stage 4 assembled fraction.

## Layout

`results/` is a copy of the shareable part of the scratch results, made by
`bash scripts/sync_results.sh [stage ...]`: reports, logs, QC summaries and tables.
Reads, contigs, graphs, BAMs and tool working dirs stay on scratch.

| Folder | Contents |
|---|---|
| `02_qc/` | FastQC raw + trimmed, fastp reports, MultiQC `read_qc.html` |
| `03_host/` | Bowtie2 logs, host-free FastQC, MultiQC `host_depletion.html` / `hostfree_qc.html` |
| `04_assembly/<run>/` | QUAST, Nonpareil curves, mapping logs, SCAPP logs; `metaspades/attempt.tsv` |
| `05_contigs/<run>/` | geNomad summaries, CheckV, MOB-typer, AMRFinderPlus, RGI, CoverM, `amr_tagged.tsv`, `amr_counts.tsv` |
| `05_contigs/<run>/plasme*/` | PLASMe reports (balance, high-precision) |
| `06_binning/<run>/` | MetaBAT2 bin lists, CheckM2 quality report, `mimag.tsv` |
| `07_mags/` | dRep tables, GTDB-Tk summaries, CoverM genome abundance |
| `side_undersequencing/` | Nonpareil curves (92 runs), model comparison, odds ratios, per-run predictions, plots |
| `tables/` | one summary table per stage |
