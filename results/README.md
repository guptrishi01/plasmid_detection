# Pilot results (stages 2–9) and side analyses

**Full run (all 92 runs, Nextflow, 2026-11-09): see [full_92/README.md](full_92/README.md).**
This page is the 4-run pilot used to tune the pipeline.

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

Assembled = own reads mapped back to own assembly. Cross-sample mapping within a binning
group: upstream ↔ upstream 18.5–20.6%, influent ↔ influent 40.7–44.6% (the earlier creek
pairs, upstream ↔ influent of the same creek, mapped only 0.4–0.8%): the same site type
shares far more of its community across creeks than different site types within a creek.
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

**Group binning (2026-11-03).** Each assembly is binned with the coverage of every run in
its binning group (`params.binning.groups`: site types pooled across both creeks and all
time points). In the pilot the groups are upstream (SRR6231181 + SRR6231217) and influent
(SRR6231184 + SRR6231191); before, each assembly was paired with the other run of its creek.
The earlier outputs are kept on scratch in `results_pre_group/`.

| Run | Bins (MetaBAT2) | High | Medium | Low | Kept (medium+) | Before (creek pairs) |
|---|---|---|---|---|---|---|
| SRR6231181 SC UP_B | 36 | 2 | 11 | 23 | 13 | 28 bins → 8 |
| SRR6231184 SC INF | 44 | 1 | 15 | 28 | 16 | 41 → 17 |
| SRR6231217 MC UP_B | 53 | 1 | 16 | 36 | 17 | 38 → 10 |
| SRR6231191 MC INF | 60 | 3 | 24 | 33 | 27 | 50 → 24 |

High = ≥ 90% complete and < 5% contaminated, medium = ≥ 50% / < 10% (CheckM2; MIMAG's rRNA
and tRNA requirement for "high" is not checked). 73 medium-or-better MAGs (was 59). The
assemblies are unchanged, so MAGs match across the two designs by shared contigs: 51 of the
73 are the earlier MAGs (contig Jaccard ≥ 0.5), 22 are new (mostly upstream), and 8 earlier
medium MAGs dropped to low. Mean contamination of kept MAGs 3.9% (was 4.2%); mean
completeness 74.9% (was 78.4%; the new MAGs are more partial). Each group still has only
two runs in the pilot; the full run gives 4–16 per group.

## Stage 7: MAG taxonomy and abundance

(`tables/07_mags.tsv`, `tables/07_group_mags.tsv`, `07_mags/groups/<group>/drep/`,
`07_mags/drep/`, `07_mags/gtdbtk/`, `07_mags/coverm_genome.tsv`)

dRep (95% ANI) in two levels: per group, upstream 30 → 26 and influent 43 → 31; then the
catalog over the 57 group MAGs → 57 (no species shared between upstream and influent;
no cross-group pair passes dRep's 90% Mash pre-screen; the most similar is ~76%). The 57 catalog MAGs (7 high, 50 medium; was
45) are all classified by GTDB-Tk 2.6.1 (GTDB r226), all bacterial: Gammaproteobacteria 21,
Actinomycetes 19, Bacteroidia 7, Clostridia 5, Bacilli 2, Campylobacteria 2, Fusobacteriia
1; 17 have no GTDB species name. Influent MAGs include gut taxa (Prevotella copri_B,
Agathobacter rectalis, Bifidobacterium adolescentis, Alistipes putredinis, Phocaeicola
vulgatus, Bacteroides uniformis). Two influent MAGs (SRR6231184.bin.38, SRR6231191.bin.45)
are both Acinetobacter johnsonii by GTDB (96.5/97.0% ANI to the reference) but 94.0% ANI
to each other over 36% of their length, so dRep keeps both; sum them for species-level
abundance.

| Run | Reads assigned to the 57 MAGs (CoverM relative abundance) | Before (45 MAGs) |
|---|---|---|
| SRR6231181 SC UP_B | 17.8% | 15.5% |
| SRR6231184 SC INF | 32.1% | 29.5% |
| SRR6231217 MC UP_B | 20.0% | 16.2% |
| SRR6231191 MC INF | 33.5% | 31.9% |

The MAGs capture about twice as much of the influent community as of the creek water,
consistent with the Nonpareil coverage and assembled fractions.

## Stage 8: MAG-level AMR and plasmid hosts

(`tables/08_mag_amr.tsv`, `tables/08_plasmids_in_mags.tsv`, `tables/08_plasmid_host_inference.tsv`)

AMRFinderPlus on the 57 catalog MAGs: 22 hits, 10 AMR-type in 5 influent MAGs: Aeromonas
caviae (SRR6231191.bin.51: blaMOX, blaOXA ×2), Phocaeicola vulgatus (SRR6231191.bin.13:
lnu(AN2), lnu(C), mef(En2)), Bacteroides uniformis (SRR6231191.bin.56: blaMUN-1, cblA),
Acinetobacter johnsonii (SRR6231184.bin.38: blaOXA), Moraxella_A (SRR6231184.bin.23:
icr-Mo); 12 stress (metal/biocide). No upstream MAG carries an AMR gene.

| geNomad plasmid contigs (4 runs) | In any MetaBAT2 bin | In a medium+ MAG (co-binning) | MOB-suite host range | Coverage correlation |
|---|---|---|---|---|
| 4,213 | 1,462 (was 1,731) | 487 (was 386) | 439 | blank (needs ≥ 5 samples; pilot groups have 2) |

Co-binned plasmid contigs per run: 11, 194, 8, 274 (SRR6231181, 184, 217, 191). A plasmid's
host MAG is the catalog MAG its bin was dereplicated into (bin → group MAG → catalog MAG).
Most frequent co-binned host genera: Acinetobacter 98, Blautia_A 48, Tolumonas 38,
Phocaeicola 35, Paenacidovorax 30, Moraxella_A 30. Of the 38 influent plasmid contigs with
AMR-type genes, 3 are co-binned: lnu(AN2)/mef(En2) with Phocaeicola vulgatus in both
influent runs (as before) and blaMUN-1 with Bacteroides uniformis (new); the earlier mph(A)
→ Bifidobacterium and sul2 → Rivihabitans links are gone. All three inferences can disagree
and none is proof.

## Stage 9: rollup, read vs contig, depth

(`tables/09_plasmid_summary.tsv`, `09_rollup/resistome_*.tsv`, `tables/09_rgi_read_vs_contig.tsv`,
`tables/09_depth_experiment.tsv`, `09_rollup/depth/plasmid_yield_vs_depth.png`, `09_rollup/stats/`)

Resistome matrices count AMR-type AMRFinderPlus hits only (`params.rollup.amrfinder_types`):
95 genes with reads (35 on mobile contigs, 74 on chromosomal ones).

RGI read level (bwt) vs contig level, genes (ARO terms) per run:

| Run | Both | Read only | Contig only |
|---|---|---|---|
| SRR6231181 SC UP_B | 0 | 130 | 7 |
| SRR6231184 SC INF | 82 | 1,030 | 22 |
| SRR6231217 MC UP_B | 0 | 140 | 10 |
| SRR6231191 MC INF | 67 | 1,046 | 20 |

Read-level counts include ARO terms with a single mapped read. The creek runs' contig-only
genes are weak (Strict) hits to models absent from RGI's read reference (folP mutation,
vanY/vanT/vanW).

Depth experiment (subsample → MEGAHIT → geNomad, same parameters as stages 4–5):

| Run | 5 M pairs | 10 M | 15 M | 20 M | Full |
|---|---|---|---|---|---|
| SRR6231181 plasmid contigs | 16 | 81 | — | — | 124 (14.7 M) |
| SRR6231184 plasmid contigs | 466 | 904 | 1,358 | 1,835 | 1,863 (20.2 M) |

No plateau: SRR6231184 gains ~450 plasmid contigs per 5 M pairs at every step, and its
20 M subsample reproduces the full-depth count within 2%. SRR6231181 has 14.7 M pairs, so
15 and 20 M were skipped.

Statistics (`09_rollup/stats/`): every ALDEx2 comparison is skipped in the pilot (1–2 runs
per group, minimum 3); the PERMANOVA on 4 runs is exploratory only (23 permutations).

## Audits

 Stage 2: all FastQC/fastp reports present, trimmed = fastp output. Stage 3:
retained = Bowtie2 input − human, R1 = R2, read ids paired. Stage 4: 18/18 jobs exit 0,
contig FASTA = QUAST counts and lengths, mapping input = retained pairs, 8/8 BAMs pass
`samtools quickcheck`. Stage 5: 16/16 jobs exit 0, every RGI contig is in the assembly,
tagged hits = counted hits. RGI logs "Requested rname SRR…_k127 does not exist" a few
times per influent run: its ORF-extension step cuts names at underscores; the hits are
still reported with full contig names. Stage 6 (group binning, chain 27651454–65, all exit 0): CheckM2 rows = bins, `pass/` =
73 = high + medium; CoverM and `depth.txt` columns = the run's group mates; `map/` keeps the
own BAM only (06c; deleted BAMs listed in `map/DELETED_BAMS.txt`). Stage 7: each group's
`Cdb.csv` = its runs' medium+ MAGs; catalog `Cdb.csv` = all group MAGs (57); 57/57
classified; `07_group_mags.tsv` = all 73 medium+ MAGs. Stage 8: one AMRFinderPlus file per
catalog MAG and no others (08a first kept 31 files from the earlier MAG set; fixed to clear
its folder, rerun 27707596), table rows = per-MAG hits (22); 08b inference rows = geNomad
plasmid contigs (4,213). Stages 4–5 and 9 outputs that don't depend on binning (assembly
table, own-read mapping logs, AMR tagging, plasmid summary, resistome, RGI read vs contig)
are byte-identical to the earlier design. Full check: `docs/slides/final_audit.py`, ALL PASS. Stage 9: 4/4 RGI bwt
outputs; no stress genes in the AMRFinderPlus resistome; plasmid summary = stage 5 counts;
09b first failed on a SIGPIPE in subsampling (fixed, rerun); all depth points present.

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
| Logistic | ×16.1 (3.3–78); influent and sewer not estimable (separation: none under-sequenced) | ×0.88 (0.07–11) | 0.725 |
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
| `08_mag_links/` | AMRFinderPlus per MAG, plasmid-in-MAG and host-inference tables |
| `09_rollup/` | resistome matrices, RGI bwt summaries, depth experiment yields and plot, statistics |
| `side_undersequencing/` | Nonpareil curves (92 runs), model comparison, odds ratios, per-run predictions, plots |
| `tables/` | one summary table per stage |
