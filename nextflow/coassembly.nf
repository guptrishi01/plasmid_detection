#!/usr/bin/env nextflow
// coassembly.nf -- co-assembly exploration next to the 11-group pipeline (PI, 2026-11-19).
// Each group of sites (params.cfg.coassembly.groups: pre_treatment, plant, stream) is pooled
// into one MEGAHIT co-assembly and binned with every run of the group as a coverage sample;
// dRep within the group gives one MAG set per group (no cross-group catalog). The group MAGs
// are compared with the per-site MAGs of the 11-group run, and each run is mapped on its own
// to its group's MAGs (what each run "hits"). Plasmid and AMR calls on the co-assemblies feed
// the plasmid-host links.
// Inputs (scripts/coassembly/prep_inputs.sbatch, from the 92-run run): host-free reads and
// per-site MAGs under params.paths.coasm_inputs; sample sheet and 11-group tables under
// params.paths.nf_results. Launch with the coasm (or coasm_test) profile:
//   PD_MAIN=coassembly.nf bash nextflow/run.sh -profile slurm,coasm_test -resume

include { BINNING; CHECKM2_MIMAG; GENOME_INFO; DREP_GROUP; GTDBTK } from './modules/binning_mags'
include { GENOMAD; CHECKV_MOBTYPER; AMR; TAG_AMR } from './modules/contigs'
include { MAG_AMR; PLASMID_HOSTS } from './modules/links_rollup'
include { COASM_MEGAHIT; COASM_INDEX; COASM_MAP; COASM_DEPTH; COASM_COLLECT; COASM_ABUNDANCE;
          COASM_COMPARE; COASM_TABLES } from './modules/coassembly'

// [key, x] channel -> one sorted list [[key, x], ...] for summary processes
def pairs(ch) {
    ch.toSortedList { a, b -> a[0] <=> b[0] }
}

// Co-assembly group of a site; every site of the sheet must be in exactly one group
def coasmGroup(site) {
    def hits = [params.cfg.coassembly.groups].flatten().findAll { spec -> site in spec.tokenize(':')[1].tokenize(',') }
    if (hits.size() != 1) {
        error "site ${site} is in ${hits.size()} params.cfg.coassembly.groups entries (need exactly 1)"
    }
    hits[0].tokenize(':')[0]
}

workflow {
    if (!params.containsKey('cfg') || !params.containsKey('paths')) {
        error "tool parameters missing: launch with nextflow/run.sh (it passes config/config.yaml)"
    }
    if (file(params.outdir).toAbsolutePath().normalize() == file(params.paths.nf_results).toAbsolutePath().normalize()) {
        error "--outdir is the 11-group results (${params.paths.nf_results}): use the coasm or coasm_test profile"
    }
    def nf = params.paths.nf_results
    def inp = params.paths.coasm_inputs
    sheetFile = file("${nf}/metadata/samples.tsv", checkIfExists: true)

    // Optional subset (coasm_test): comma-separated runs or a file, one per line
    def subset = null
    if (params.runs) {
        def f = file(params.runs.toString())
        subset = (f.exists() ? f.readLines() : params.runs.toString().tokenize(',')).collect { x -> x.trim() }.findAll { x -> x }
    }

    metas = channel.fromPath(sheetFile)
        .splitCsv(sep: '\t', header: true)
        .filter { r -> subset == null || r.run in subset }
        .map { r -> [id: r.run, location: r.location, site: r.site, timepoint: r.timepoint] }
    reads = metas.map { m ->
        [coasmGroup(m.site), m, file("${inp}/hostfree/${m.id}_1.fastq.gz", checkIfExists: true),
         file("${inp}/hostfree/${m.id}_2.fastq.gz", checkIfExists: true)] }

    // ---- pool -> co-assembly -> one independent mapping per run -----------------------------
    pooled = reads.groupTuple(by: 0)
        .map { g, ms, r1s, r2s ->
            def o = (0..<ms.size()).toList().sort { i -> ms[i].id }       // sorted by run: stable cache key
            [g, o.collect { i -> ms[i].id }, o.collect { i -> r1s[i] }, o.collect { i -> r2s[i] }] }
    COASM_MEGAHIT(pooled)
    contigs = COASM_MEGAHIT.out.contigs
    COASM_INDEX(contigs)
    COASM_MAP(reads.combine(COASM_INDEX.out.idx, by: 0))
    groupBams = COASM_MAP.out.bam.groupTuple(by: 0)
        .map { g, ids, bams, bais -> [g, ids, bams, bais] }
    COASM_DEPTH(contigs.join(groupBams))

    // ---- binning, quality, one MAG set per group --------------------------------------------
    BINNING(contigs.join(COASM_DEPTH.out.depth).map { g, c, d -> [[id: g], c, d] })
    CHECKM2_MIMAG(BINNING.out.bins)
    mimag = pairs(CHECKM2_MIMAG.out.mimag.map { meta, m -> [meta.id, m] })
    GENOME_INFO(mimag)
    DREP_GROUP(CHECKM2_MIMAG.out.pass.map { meta, p -> [meta.id, [p]] }, GENOME_INFO.out)
    groupDreps = pairs(DREP_GROUP.out.drep)
    COASM_COLLECT(groupDreps)
    mags = COASM_COLLECT.out.mags
    GTDBTK(mags)
    MAG_AMR(mags)

    // ---- each run on its own against its group's MAGs ---------------------------------------
    COASM_ABUNDANCE(DREP_GROUP.out.drep.join(pooled.map { g, ids, r1s, r2s -> [g, ids, r1s + r2s] }))

    // ---- plasmids and AMR on the co-assemblies, plasmid hosts -------------------------------
    cm = contigs.map { g, c -> [[id: g], c] }
    GENOMAD(cm)
    CHECKV_MOBTYPER(GENOMAD.out.genomad)
    AMR(cm)
    TAG_AMR(GENOMAD.out.genomad.join(AMR.out.amr))
    PLASMID_HOSTS(pairs(BINNING.out.bins.map { meta, b -> [meta.id, b] }), mimag,
                  pairs(GENOMAD.out.genomad.map { meta, x -> [meta.id, x] }),
                  pairs(CHECKV_MOBTYPER.out.mob.map { meta, x -> [meta.id, x] }),
                  pairs(COASM_DEPTH.out.coverm), mags, GTDBTK.out.gtdbtk, groupDreps)

    // ---- comparison with the per-site MAGs, summary tables ----------------------------------
    COASM_COMPARE(mags, file("${inp}/site_mags", checkIfExists: true))
    COASM_TABLES(sheetFile, mimag, pairs(COASM_ABUNDANCE.out.coverm), COASM_COLLECT.out[1], GTDBTK.out.gtdbtk,
                 MAG_AMR.out[1], COASM_COMPARE.out.tables, COASM_COMPARE.out[1],
                 file("${nf}/tables/07_group_mags.tsv", checkIfExists: true),
                 file("${nf}/tables/08_mag_amr.tsv", checkIfExists: true))
}
