# Full run: all 92 WGS runs (Nextflow, group binning)

All 92 shotgun metagenomes of SRP121672 through stages 1–9 with the Nextflow pipeline
(`nextflow/`, head job 27710034, 2026-11-08 → 11-09, 22 h; 1,444 tasks run, 52 reused from
the pilot, 0 failed). Binning groups as confirmed by the PI on 2026-11-03
(`params.binning.groups`). metaSPAdes was skipped (`--metaspades false`; tested on the
pilot). Every number below comes from the files in this folder; reads, contigs, bins,
BAMs and depth tables stay on scratch (`/scratch/rgupta25/plasmid_detection/results_nf`).
Synced with `scripts/sync_results_nf.sh`. Full audit: 54/54 checks pass (see Audit).

## Stages 2–4: reads and assemblies

92 runs, 15.1–30.9 M raw read pairs; 76.8–96.9% of pairs kept after trimming and human
removal (1,576–188,091 human pairs per run). MEGAHIT: 9,653–95,110 contigs ≥ 1 kb per run,
16–250 Mb, N50 1,603–5,048 bp; 3.6–75.9% of each run's reads map back to its own assembly.
For the 4 pilot runs, the stage 2–3 tables equal the sbatch pilot exactly.

## Plasmids, resistance and MAGs by binning group

| Group | Runs | Plasmid contigs (median per run) | Plasmid contigs with AMR genes (all runs) | Conjugative / mobilizable (all runs) | SCAPP circular plasmids | Assembled (median) | Reads in catalog MAGs (median) | Medium+ MAGs | Group MAGs (dRep) |
|---|---|---|---|---|---|---|---|---|---|
| upstream | 16 | 130 | 7 | 1 / 12 | 156 | 26.1% | 52.6% | 299 | 161 |
| downstream | 16 | 146 | 14 | 2 / 14 | 239 | 51.7% | 54.2% | 258 | 117 |
| HOSP | 8 | 1,668 | 191 | 21 / 738 | 310 | 48.5% | 54.9% | 204 | 58 |
| RES | 8 | 1,597 | 177 | 20 / 629 | 254 | 46.8% | 55.7% | 207 | 67 |
| influent | 8 | 1,820 | 190 | 27 / 729 | 315 | 48.1% | 57.1% | 201 | 61 |
| PCI | 4 | 1,294 | 82 | 10 / 279 | 117 | 49.6% | 62.6% | 85 | 37 |
| PCE | 8 | 1,580 | 182 | 27 / 674 | 276 | 51.8% | 60.2% | 193 | 65 |
| ATE | 8 | 722 | 31 | 2 / 65 | 140 | 37.4% | 58.5% | 144 | 91 |
| FCE | 8 | 638 | 56 | 17 / 169 | 175 | 59.0% | 69.1% | 128 | 96 |
| UV | 4 | 382 | 6 | 2 / 38 | 102 | 46.9% | 46.2% | 81 | 69 |
| control | 4 | 150 | 0 | 0 / 0 | 30 | 27.9% | 46.1% | 48 | 48 |

(`tables/09_plasmid_summary.tsv`, `04_assembly.tsv`, `07_mags.tsv`, `07_group_mags.tsv`.)
Pre-treatment sites (hospital and residential sewer, influent, PCI, PCE) carry ~10× the
plasmid contigs of creek water and ~35× the plasmid contigs with AMR-type genes per run
(22.8 vs 0.66); both fall
through the plant steps (ATE → FCE → UV). Plasmid contigs with AMR genes = contigs carrying
an AMRFinderPlus AMR-type or RGI hit (`params.rollup.amrfinder_types`).

## Stage 7: MAG catalog

3,927 bins → 1,848 medium-or-better (CheckM2 + MIMAG thresholds) → 870 group MAGs (dRep
95% ANI within each group) → **515 catalog MAGs** (154 high, 361 medium; 355 group MAGs
merged across groups). GTDB-Tk r226: 514 bacterial + 1 unclassified bacterium; 242 have
no GTDB species name. Classes: Gammaproteobacteria 170, Actinomycetes 99, Bacteroidia 82,
Alphaproteobacteria 43, Minisyncoccia 19 (Patescibacteria), Polyangia 11. 329 catalog MAGs
were recovered in one group only, 186 in two or more (the most widespread, e.g.
Aliarcobacter, Prevotella, Bacteroides, Rivihabitans, in 6–7 pre-treatment/plant groups).
The catalog captures 30.5–76.0% of each run's reads.

## Stage 8: resistance genes in MAGs and plasmid hosts

AMRFinderPlus on the 515 catalog MAGs: 53 hits, 46 AMR-type (beta-lactam 23, macrolide
12, …) in 26 MAGs, 7 stress. **No creek or control MAG carries an AMR gene.** Hosts:
Aeromonas (blaMOX, blaOXA; influent, HOSP, PCE, FCE), Cloacibacterium (ere(D), mef(C),
mph(G), blaMCA; RES, PCI, PCE), Pseudomonas (blaPDC, aph(3')-IIb, catB7, mexA; FCE),
Phocaeicola (lnu(AN2), mef(En2), erm(G), aadS; RES), Enterobacter (blaACT, oqxB; HOSP),
Moraxella_A (blaCARB-16, icr-Mo), Bacteroides (blaMUN-1, cblA), Acinetobacter (blaOXA),
Propionivibrio, Laribacter, Flavobacterium, Streptococcus, Brachymonas, Casimicrobium.

| geNomad plasmid contigs | In a medium+ MAG (co-binning) | MOB-suite host range | Coverage correlation (best MAG) | ρ ≥ 0.8 |
|---|---|---|---|---|
| 78,716 | 15,901 | 7,408 | 70,862 | 53,586 |

Coverage correlation now runs (groups of 8–16 runs); it is blank for the 12 runs in the
4-run groups (PCI, UV, control), by design. Where a plasmid has both a co-binned host and
a correlation host (14,168), they are the same MAG for 9,219 (65%). Top co-binned host
genera: Acinetobacter, Acidovorax, Blautia_A, Mycobacterium, Phocaeicola, Rivihabitans;
top hosts at ρ ≥ 0.8: Acinetobacter, Moraxella_A, Streptococcus, Acidovorax,
Cloacibacterium, Agathobacter. All three are inferences, not proof.

## Stage 9: resistome statistics

Resistome: 175 AMR-type genes with reads across the 92 runs (`09_rollup/resistome_*`).
22 runs have no reads on any AMR-type gene of their own assembly: 14 upstream, 4
downstream and all 4 controls.

- PERMANOVA (exploratory, Aitchison, 88 non-control runs): site group explains 40.5% of
  resistome variation (p = 0.0001), time point 4.0% (p = 0.037), creek system 2.0% (p = 0.029).
- ALDEx2 (BH): **no gene significant** in any comparison. Upstream vs downstream within
  each creek (8 vs 8; 3–6 genes with reads) and within time points 2–4 (4 vs 4); time point 1
  skipped (1 gene with reads). Influent vs plant within each creek (4 vs 16; 122–135 genes);
  within time points skipped (2 influent runs each, minimum 3).
- **Stats groups are provisional**: `params.stats.site_groups` still puts PCI and PCE in
  "plant", though they are pre-treatment; the comparisons are to be confirmed with the PI.
  Changing them needs only the statistics step re-run.

The Nextflow STATS step of the 92-run run passed no comparisons to `stats.R` (fixed in
33e1af4); the statistics here were produced by the same commands as a standalone job
(`/scratch/.../stats_rerun/run_stats.sbatch`, job 27725774), after a `stats.R` fix for tests
with fewer than 2 genes (325cc5a).

## Audit

`docs/slides/audit_full.py` (local), 54/54 checks pass: 92 WGS runs, all PASS in the
download manifest; stage 2–3 identities (trimmed = fastp, trimmed − human = retained,
host-free FastQC R1 = R2 = retained) for all 92 and equal to the sbatch pilot for its 4
runs; assembly table = QUAST and own Bowtie2 logs; mapping logs, CoverM, read-count and
depth columns = exactly the run's binning group; MOB-typer rows = geNomad plasmids;
tagged = AMRFinderPlus + RGI = counted; MIMAG tiers follow the thresholds; each group's
dRep input = its runs' medium+ MAGs; catalog = all group MAGs = 07_mags = GTDB-Tk; one
AMRFinderPlus file per catalog MAG; host inference covers all 78,716 plasmids, every host
is a catalog MAG, correlation blank exactly for the 4-run groups; plasmid summary = stage
5 counts; no stress genes in the resistome; RGI bwt 92/92; every comparison in the log.

## Layout

`tables/` (all summary tables), `metadata/` (NCBI run table, sample sheet, BioSample
antibiotics), `pipeline_info/` (Nextflow report, timeline, trace), `02_qc/multiqc`,
`03_host/multiqc`, per run `04_assembly/<run>/{quast/report.tsv,map/*.bowtie2.log}`,
`05_contigs/<run>/{amr_tagged.tsv,mob_typer.tsv,checkv/quality_summary.tsv}`,
`06_binning/<run>/{mimag.tsv,checkm2/}`, `07_mags/` (dRep tables per group and catalog,
GTDB-Tk, CoverM), `08_mag_links/amrfinder/`, `09_rollup/` (resistome, stats, RGI bwt).
