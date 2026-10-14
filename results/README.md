# Pilot results (stages 2–5)

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

## Audits

 Stage 2: all FastQC/fastp reports present, trimmed = fastp output. Stage 3:
retained = Bowtie2 input − human, R1 = R2, read ids paired. Stage 4: 18/18 jobs exit 0,
contig FASTA = QUAST counts and lengths, mapping input = retained pairs, 8/8 BAMs pass
`samtools quickcheck`. Stage 5: 16/16 jobs exit 0, every RGI contig is in the assembly,
tagged hits = counted hits. RGI logs "Requested rname SRR…_k127 does not exist" a few
times per influent run: its ORF-extension step cuts names at underscores; the hits are
still reported with full contig names.

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
| `tables/` | one summary table per stage |
