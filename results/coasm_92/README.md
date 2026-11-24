# Co-assembly: MAGs per group (pre-treatment, plant, stream)

PI direction (2026-11-19): pool the reads of each of three groups into one co-assembly, find
the group's MAGs, compare them with the per-site MAGs of the full run (`../full_92`), then
map every run on its own to its group's MAGs and see what it hits. Workflow
`nextflow/coassembly.nf` (profile `coasm`, see `nextflow/README.md`), `--outdir`
`/scratch/rgupta25/plasmid_detection/results_nf_coasm`; inputs are the host-free reads and
the 870 per-site MAGs of the full run (hard-linked, nothing re-run upstream). Every number
below comes from the files in this folder (synced with `scripts/sync_results_nf.sh`);
contigs, bins, BAMs and depth tables stay on scratch. Audit: 67/67 checks pass.

| Group | Sites | Runs | Host-free read pairs |
|---|---|---|---|
| pre_treatment | HOSP, RES, INF, PCI | 28 | 587 M |
| plant | PCE, ATE, FCE, UV | 28 | 582 M |
| stream | UP_A/B, DS_A/B, MNT_A/B, UW_A/B | 36 | 725 M |

## Co-assemblies

MEGAHIT 1.2.9, `--presets meta-large` (k 27–127), contigs ≥ 1 kb
(`04_assembly/<group>/assembly_stats.tsv`, `runs.txt`, per-run `map/*.bowtie2.log`).

| Group | Contigs | Total | N50 | Largest | MEGAHIT time / peak RAM | Runs' reads aligned (Bowtie2) |
|---|---|---|---|---|---|---|
| pre_treatment | 1,024,058 | 2.62 Gbp | 3,018 | 681 kb | 35 h / 160 GB | 74.4–84.2% |
| plant | 1,492,466 | 3.99 Gbp | 3,229 | 1.31 Mb | 46 h / 158 GB | 53.7–88.8% |
| stream | 1,261,376 | 3.05 Gbp | 2,740 | 596 kb | 72 h / 231 GB | 5.0–77.9% |

The two lowest stream runs are the time-point-4 controls MNT_B (5.0%) and UW_A (9.2%): their
communities are largely absent from a co-assembly dominated by the 34 Sugar/Mallard Creek
runs.

## MAGs per group

MetaBAT2 on each co-assembly with every run of the group as a coverage sample, CheckM2
MIMAG tiers, dRep at 95% ANI within the group (one MAG per species), GTDB-Tk
(`06_binning/<group>/mimag.tsv`, `07_mags/`, `tables/coasm_mags.tsv`).

| Group | Bins | Medium+ | MAGs (high) | No GTDB species | Top phyla |
|---|---|---|---|---|---|
| pre_treatment | 529 | 190 | **190** (44) | 31 | Bacillota 71, Pseudomonadota 50, Bacteroidota 37 |
| plant | 1,096 | 507 | **507** (155) | 291 | Pseudomonadota 205, Bacteroidota 80, Patescibacteriota 65 |
| stream | 1,056 | 369 | **369** (65) | 225 | Pseudomonadota 130, Actinomycetota 121, Bacteroidota 71 |

No two medium+ bins of a group were the same species, so dRep kept all of them.
Stream MAGs are freshwater taxa (Rhodoluna 53, Planktophila 27, Limnohabitans 37,
Polynucleobacter 20); pre-treatment MAGs are gut and sewer taxa (Prevotella, Acinetobacter,
Bacteroides, Bifidobacterium).

## Group MAGs vs per-site MAGs

dRep compare at 95% ANI over the 1,066 group MAGs and the 870 per-site MAGs (level-1 dRep
of the 11 site groups of the full run) (`tables/coasm_compare.tsv`,
`coasm_compare_mags.tsv`).

| Group | Group species | Per-site species (MAGs) | Shared | Group only | Per-site only | Completeness, shared pairs (median group − site) | Resistant per-site MAGs recovered |
|---|---|---|---|---|---|---|---|
| pre_treatment | 190 | 101 (223) | 49 | 141 | 52 | −0.4 | 5 / 14 |
| plant | 507 | 234 (321) | 185 | 322 | 49 | +1.3 | 4 / 12 |
| stream | 369 | 248 (326) | 151 | 218 | 97 | −1.6 | 0 / 0 |

Co-assembly adds 681 species that no per-site MAG had (36 of them were binned per-site in
another group's sites); per-site binning keeps 198 species the co-assemblies did not
recover. Where both have a species, genome quality is about the same. Most resistant per-site MAGs were not recovered by the groups.

## What each run hits

CoverM maps each run separately to its own group's MAGs; a read counts at ≥ 95% identity
over ≥ 75% of its length (`params.coassembly.coverm_read_filters`), a MAG is hit at mean
coverage ≥ 1× (`tables/coasm_run_hits.tsv`). Medians per site:

| Site | Runs | % reads on group MAGs | MAGs hit | Resistant MAGs hit |
|---|---|---|---|---|
| HOSP | 8 | 17.1 | 84 | 7.5 |
| RES | 8 | 17.8 | 77 | 6.5 |
| INF | 8 | 16.4 | 75 | 8 |
| PCI | 4 | 17.7 | 65 | 6.5 |
| PCE | 8 | 15.7 | 71 | 8.5 |
| ATE | 8 | 34.2 | 104 | 1 |
| FCE | 8 | 19.6 | 57 | 1.5 |
| UV | 4 | 23.9 | 63 | 1 |
| UP_A / UP_B | 8 / 8 | 11.7 / 12.8 | 62 / 69 | 0 |
| DS_A / DS_B | 8 / 8 | 12.2 / 14.5 | 52 / 66 | 0 |
| MNT_A / MNT_B | 1 / 1 | 17.6 / 1.7 | 43 / 2 | 0 |
| UW_A / UW_B | 1 / 1 | 1.3 / 19.7 | 12 / 57 | 0 |

PCE (primary clarifier effluent, now a plant step) still looks like pre-treatment: it hits
as many resistant MAGs as influent. From aeration on, resistant MAGs drop to ~1 per run.
Most abundant MAGs: Cloacibacterium, *Phocaeicola vulgatus*, *Streptococcus suis*
(pre-treatment); Thiolinea, Gordonia, Aquella (plant); Zoogloea, Limnohabitans, Rhodoluna
(stream). Every group MAG is hit by at least one run.

Without the read filters (CoverM default: any alignment counts) the medians were 44.8 /
58.4 / 44.4% of reads and 130 / 102 / 126 MAGs hit, and the control runs with 5–9%
Bowtie2 alignment showed 18–25% of reads on MAGs, so the filters are used.

## Resistance genes in group MAGs

AMRFinderPlus (`--plus`, AMR type only) on the group MAGs (`tables/08_mag_amr.tsv`): 42
hits in 31 MAGs (12 pre-treatment, 18 plant, 1 stream).

| Group | MAG | Genes | Runs hit |
|---|---|---|---|
| pre_treatment | *Bacteroides uniformis* | cblA | 28 / 28 |
| pre_treatment | *Aeromonas taiwanensis* | blaMOX, blaOXA | 28 / 28 |
| pre_treatment | *Laribacter hongkongensis* | blaLHK | 27 / 28 |
| pre_treatment | *Enterobacter cloacae* | oqxA, oqxB | 23 / 28 |
| pre_treatment | *Acinetobacter junii* | ant(3'')-IIa | 18 / 28 |
| pre_treatment | *Sebaldella termitidis* | erm(52) | 17 / 28 |
| plant | *Phocaeicola vulgatus* | erm(G), lnu(AN2), mef(En2), msr(D) | 12 / 28 |
| plant | *Aeromonas media* | blaMOX, blaOXA | 10 / 28 |
| plant | *Kluyvera intermedia* | blaCTX-M | 4 / 28 |
| plant | JALHNF01 sp. | blaVIM | 3 / 28 |
| stream | *Rhodococcus qingshengii* | iri | 7 / 36 |

The same species carry the same genes in pre-treatment and plant (*B. uniformis*,
*Aeromonas*, *L. hongkongensis*, *A. bouvetii*, *Mangrovibacter*, *Dysgonomonas*) but are
hit in fewer plant runs. The plant co-assembly has genomes the per-site binning missed,
among them a blaCTX-M *Kluyvera* and a blaVIM carrier. The stream has a single resistant
MAG (rifampin monooxygenase iri, intrinsic to *Rhodococcus*).

## Plasmids on the co-assemblies

geNomad, MOB-typer, plasmid-host inference (`05_contigs/<group>/`,
`tables/08_plasmid_host_inference.tsv`).

| Group | Plasmid contigs | Mobilizable / conjugative | In a MAG (co-binned) | Coverage-correlated MAG (ρ ≥ 0.8) |
|---|---|---|---|---|
| pre_treatment | 24,081 | 420 / 18 | 1,839 | 10,792 |
| plant | 22,508 | 316 / 26 | 5,765 | 17,090 |
| stream | 7,861 | 44 / 5 | 2,027 | 2,601 |

Correlation hosts use all runs of the group (28–36 samples) and are weak evidence on their
own (many contigs correlate with some MAG); co-binning is the stronger link.

## Audit

`docs/slides/audit_coasm.py` (local), 67/67 checks pass on all 92 runs: group membership =
config; pooled runs and mapping logs = the group's runs; every run mapped with all its
host-free pairs; bins, tiers and dRep follow the thresholds; group MAG lists = dRep = GTDB-Tk =
CoverM genomes; run-hit and MAG tables recomputed from CoverM without
`coasm_tables.py`; comparison recomputed from dRep `Cdb`/`Ndb`; one host row per geNomad
plasmid, co-binned hosts are MAGs of the plasmid's own group, correlation over every run
of the group; the full-run tables (`results_nf/tables`) unchanged. The co-assembly AMR outputs
were produced twice (a resume re-ran them) and are byte-identical.

## Layout

`tables/` (summary tables), `pipeline_info/`, `04_assembly/<group>/` (assembly stats, run
list, MEGAHIT log, per-run Bowtie2 logs), `05_contigs/<group>/` (AMR tags, MOB-typer,
CheckV), `06_binning/<group>/` (MIMAG, CheckM2), `07_mags/` (group MAG lists, dRep per
group, comparison, GTDB-Tk, CoverM per group), `08_mag_links/amrfinder/`.
