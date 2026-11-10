#!/usr/bin/env nextflow
// plasmid_detection -- Nextflow version of scripts/01-09 (see nextflow/README.md).
// Launch with nextflow/run.sh (it passes config/config.yaml as -params-file).

include { FETCH_RUNINFO; SAMPLE_SHEET; DOWNLOAD_READS } from './modules/fetch'
include { FASTQC_RAW; FASTP; FASTQC_TRIMMED; QC_SUMMARY;
          HOST_DEPLETION; FASTQC_HOSTFREE; HOST_SUMMARY } from './modules/qc_host'
include { MEGAHIT; METASPADES; MAP_GROUP; ASSEMBLY_QC; SCAPP; ASSEMBLY_TABLE } from './modules/assembly'
include { GENOMAD; CHECKV_MOBTYPER; AMR; TAG_AMR } from './modules/contigs'
include { BINNING; CHECKM2_MIMAG; GENOME_INFO; DREP_GROUP; DREP_CATALOG; GTDBTK; MAG_ABUNDANCE } from './modules/binning_mags'
include { MAG_AMR; PLASMID_HOSTS; RGI_BWT; ROLLUP; STATS;
          DEPTH_SUBSAMPLE; DEPTH_SUMMARY } from './modules/links_rollup'

// [meta, x, ...] channel -> one sorted list [[run, x], ...] for summary processes
def pairs(ch) {
    ch.map { t -> [t[0].id, t[1]] }.toSortedList { a, b -> a[0] <=> b[0] }
}

// On/off options: a command-line value arrives as a string ("--metaspades false" is the
// non-empty, so truthy, String "false"); only true / "true" switch an option on
def flag(v) {
    v in [true, 'true']
}

// Binning group of a site (params.binning.groups, "group:SITE,SITE"), looked up where needed
// rather than stored in meta: meta is part of every task's cache key, so adding a field would
// re-run every finished task on -resume. A SampleName the regex
// doesn't parse (site NA, another study) puts every run in one group, "all"; a parsed site
// missing from the groups is an error, as in the sbatch version.
def binningGroup(site) {
    def g = [params.cfg.binning.groups].flatten().find { spec -> site in spec.tokenize(':')[1].tokenize(',') }
    if (g) {
        return g.tokenize(':')[0]
    }
    if (site == 'NA') {
        return 'all'
    }
    error "site ${site} is in no params.binning.groups entry (config/config.yaml)"
}

// [key, [run, ...]] channel -> [key, [[run, ...], ...]] once every run of that key has
// arrived (nByKey: value channel {key: expected count}), sorted by run: groupTuple keeps
// arrival order, which changes between runs, and list order is part of the task cache key
def groupByKey(ch, nByKey) {
    ch.combine(nByKey)
        .map { k, x, n -> [groupKey(k, n[k]), x] }
        .groupTuple()
        .map { k, xs -> [k.toString(), xs.sort(false) { x -> x[0] }] }
}

workflow {
    if (!params.containsKey('cfg')) {
        error "tool parameters missing: launch with nextflow/run.sh (it passes config/config.yaml)"
    }

    // ---- stage 1: run list, sample sheet, reads -----------------------------------------
    runinfo = FETCH_RUNINFO(params.id, params.library_strategy, params.runs ?: '')
    sheet = SAMPLE_SHEET(runinfo, params.sample_name_regex)
    metas = sheet
        .splitCsv(sep: '\t', header: true)
        .map { r -> [id: r.run, biosample: r.biosample, location: r.location, site: r.site,
                     timepoint: r.timepoint, spots: r.spots as long] }
    raw = DOWNLOAD_READS(metas)

    // ---- stages 2-3 ------------------------------------------------------------------------
    FASTQC_RAW(raw)
    FASTP(raw)
    FASTQC_TRIMMED(FASTP.out.reads)
    runs = metas.map { m -> m.id }.toSortedList()
    QC_SUMMARY(runs, FASTQC_RAW.out.collect(), FASTP.out.json.collect(), FASTQC_TRIMMED.out.collect())

    HOST_DEPLETION(FASTP.out.reads)
    hostfree = HOST_DEPLETION.out.reads
    FASTQC_HOSTFREE(hostfree)
    HOST_SUMMARY(runs, HOST_DEPLETION.out.log.collect(), HOST_DEPLETION.out.counts.collect(),
                 FASTP.out.json.collect(), FASTQC_HOSTFREE.out.collect())
    counts = HOST_DEPLETION.out.counts.map { f -> [f.name.replace('.counts.tsv', ''), f] }

    // ---- stage 4 ---------------------------------------------------------------------------
    MEGAHIT(hostfree)
    contigs = MEGAHIT.out.contigs

    if (flag(params.metaspades)) {
        pick = params.cfg.assembly.metaspades_run
        one = pick ? hostfree.filter { t -> t[0].id == pick }
                   : hostfree.toSortedList { a, b -> a[0].spots <=> b[0].spots }.map { l -> l[0] }
        METASPADES(one)
    }

    ASSEMBLY_QC(contigs.join(hostfree))
    SCAPP(MEGAHIT.out.graph.join(hostfree))

    // ---- stage 5 ---------------------------------------------------------------------------
    GENOMAD(contigs)
    CHECKV_MOBTYPER(GENOMAD.out.genomad)
    AMR(contigs)
    TAG_AMR(GENOMAD.out.genomad.join(AMR.out.amr))

    // ---- stage 4c + 5 coverage: map each assembly with every run of its binning group ------
    // Reads grouped as soon as all of a group's runs are host-depleted (groupKey carries the
    // expected count); MAP_GROUP also needs the run's AMR hits (TAG_AMR) for the read counts
    nPerGroup = metas.map { m -> binningGroup(m.site) }.toList().map { l -> l.countBy { x -> x } }
    groupReads = groupByKey(hostfree.map { meta, r1, r2 -> [binningGroup(meta.site), [meta.id, r1, r2]] }, nPerGroup)
        .map { g, xs -> [g, xs.collect { x -> x[0] }, xs.collect { x -> x[1] }, xs.collect { x -> x[2] }] }
    MAP_GROUP(contigs.join(TAG_AMR.out.tagged.map { meta, _t, bed -> [meta, bed] })
                     .map { meta, c, bed -> [binningGroup(meta.site), meta, c, bed] }
                     .combine(groupReads, by: 0)
                     .map { _g, meta, c, bed, ids, r1s, r2s -> [meta, c, bed, ids, r1s, r2s] })

    ASSEMBLY_TABLE(pairs(ASSEMBLY_QC.out.report), pairs(MAP_GROUP.out.own_log),
                   pairs(MEGAHIT.out.graph.map { meta, _g, k -> [meta, k] }))

    // ---- stage 6 ---------------------------------------------------------------------------
    BINNING(contigs.join(MAP_GROUP.out.depth))
    CHECKM2_MIMAG(BINNING.out.bins)
    mimag = pairs(CHECKM2_MIMAG.out.mimag)

    // ---- stage 7: dRep per binning group, then across groups (catalog) ---------------------
    GENOME_INFO(mimag)
    DREP_GROUP(groupByKey(CHECKM2_MIMAG.out.pass.map { meta, p -> [binningGroup(meta.site), [meta.id, p]] }, nPerGroup)
                   .map { g, xs -> [g, xs.collect { x -> x[1] }] },
               GENOME_INFO.out)
    groupDreps = DREP_GROUP.out.drep.toSortedList { a, b -> a[0] <=> b[0] }
    DREP_CATALOG(groupDreps, GENOME_INFO.out)
    catalog = DREP_CATALOG.out.drep
    GTDBTK(catalog)
    MAG_ABUNDANCE(catalog, GTDBTK.out.gtdbtk, mimag, groupDreps,
                  hostfree.flatMap { _meta, r1, r2 -> [r1, r2] }.collect())

    // ---- stage 8 ---------------------------------------------------------------------------
    MAG_AMR(catalog)
    PLASMID_HOSTS(pairs(BINNING.out.bins), mimag, pairs(GENOMAD.out.genomad),
                  pairs(CHECKV_MOBTYPER.out.mob), pairs(MAP_GROUP.out.coverm), catalog, GTDBTK.out.gtdbtk,
                  groupDreps)

    // ---- stage 9 ---------------------------------------------------------------------------
    RGI_BWT(hostfree)
    ROLLUP(pairs(GENOMAD.out.genomad), pairs(CHECKV_MOBTYPER.out.mob),
           pairs(TAG_AMR.out.tagged.map { meta, t, _bed -> [meta, t] }),
           pairs(MAP_GROUP.out.counts),
           pairs(SCAPP.out.scapp), pairs(RGI_BWT.out.genes))
    STATS(sheet, runs, ROLLUP.out.resistome)

    if (flag(params.depth_experiment)) {
        dRuns = [params.cfg.rollup.depth_samples].flatten()
        levels = channel.fromList([params.cfg.rollup.depth_levels_million_pairs].flatten())
        subs = hostfree
            .filter { t -> t[0].id in dRuns }
            .map { meta, r1, r2 -> [meta.id, meta, r1, r2] }
            .join(counts)
            .map { _id, meta, r1, r2, c -> [meta, r1, r2, c] }
            .combine(levels)
        DEPTH_SUBSAMPLE(subs)
        full = contigs.filter { t -> t[0].id in dRuns }
            .join(GENOMAD.out.genomad)
            .map { meta, c, g -> [meta.id, c, g] }
            .join(counts)
            .toSortedList { a, b -> a[0] <=> b[0] }
        DEPTH_SUMMARY(DEPTH_SUBSAMPLE.out.yield.collect(), full)
    }
}
