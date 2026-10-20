#!/usr/bin/env nextflow
// plasmid_detection -- Nextflow version of scripts/01-09 (see nextflow/README.md).
// Launch with nextflow/run.sh (it passes config/config.yaml as -params-file).

include { FETCH_RUNINFO; SAMPLE_SHEET; DOWNLOAD_READS } from './modules/fetch'
include { FASTQC_RAW; FASTP; FASTQC_TRIMMED; QC_SUMMARY;
          HOST_DEPLETION; FASTQC_HOSTFREE; HOST_SUMMARY } from './modules/qc_host'
include { MEGAHIT; METASPADES; MAP_SYSTEM; ASSEMBLY_QC; SCAPP; ASSEMBLY_TABLE } from './modules/assembly'
include { GENOMAD; CHECKV_MOBTYPER; AMR; COVERAGE_TAG } from './modules/contigs'
include { BINNING; CHECKM2_MIMAG; DREP; GTDBTK; MAG_ABUNDANCE } from './modules/binning_mags'
include { MAG_AMR; PLASMID_HOSTS; RGI_BWT; ROLLUP; STATS;
          DEPTH_SUBSAMPLE; DEPTH_SUMMARY } from './modules/links_rollup'

// [meta, x, ...] channel -> one sorted list [[run, x], ...] for summary processes
def pairs(ch) {
    ch.map { t -> [t[0].id, t[1]] }.toSortedList { a, b -> a[0] <=> b[0] }
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

    if (params.metaspades) {
        pick = params.cfg.assembly.metaspades_run
        one = pick ? hostfree.filter { t -> t[0].id == pick }
                   : hostfree.toSortedList { a, b -> a[0].spots <=> b[0].spots }.map { l -> l[0] }
        METASPADES(one)
    }

    // Reads of every run in the same system (location), grouped as soon as that system's
    // runs are all host-depleted (groupKey carries the expected count)
    perLoc = metas.map { m -> m.location }.toList().map { l -> l.countBy { x -> x } }
    byLoc = hostfree
        .combine(perLoc)
        .map { meta, r1, r2, n -> [groupKey(meta.location, n[meta.location]), meta.id, r1, r2] }
        .groupTuple()
        .map { loc, ids, r1s, r2s -> [loc.toString(), ids, r1s, r2s] }
    MAP_SYSTEM(contigs.map { meta, c -> [meta.location, meta, c] }
                      .combine(byLoc, by: 0)
                      .map { _loc, meta, c, ids, r1s, r2s -> [meta, c, ids, r1s, r2s] })
    bams = MAP_SYSTEM.out.bams

    ASSEMBLY_QC(contigs.join(hostfree))
    SCAPP(MEGAHIT.out.graph.join(hostfree))
    ASSEMBLY_TABLE(pairs(ASSEMBLY_QC.out.report), pairs(MAP_SYSTEM.out.own_log),
                   pairs(MEGAHIT.out.graph.map { meta, _g, k -> [meta, k] }))

    // ---- stage 5 ---------------------------------------------------------------------------
    GENOMAD(contigs)
    CHECKV_MOBTYPER(GENOMAD.out.genomad)
    AMR(contigs)
    COVERAGE_TAG(bams.join(GENOMAD.out.genomad).join(AMR.out.amr))

    // ---- stage 6 ---------------------------------------------------------------------------
    BINNING(contigs.join(bams))
    CHECKM2_MIMAG(BINNING.out.bins)
    mimag = pairs(CHECKM2_MIMAG.out.mimag)

    // ---- stage 7 ---------------------------------------------------------------------------
    DREP(mimag, CHECKM2_MIMAG.out.pass.flatten().collect())
    GTDBTK(DREP.out.drep)
    MAG_ABUNDANCE(DREP.out.drep, GTDBTK.out.gtdbtk, mimag,
                  hostfree.flatMap { _meta, r1, r2 -> [r1, r2] }.collect())

    // ---- stage 8 ---------------------------------------------------------------------------
    MAG_AMR(DREP.out.drep)
    coverm = COVERAGE_TAG.out.tagged.map { meta, _t, _c, cov -> [meta, cov] }
    PLASMID_HOSTS(pairs(BINNING.out.bins), mimag, pairs(GENOMAD.out.genomad),
                  pairs(CHECKV_MOBTYPER.out.mob), pairs(coverm), DREP.out.drep, GTDBTK.out.gtdbtk)

    // ---- stage 9 ---------------------------------------------------------------------------
    RGI_BWT(hostfree)
    ROLLUP(pairs(GENOMAD.out.genomad), pairs(CHECKV_MOBTYPER.out.mob),
           pairs(COVERAGE_TAG.out.tagged.map { meta, t, _c, _cov -> [meta, t] }),
           pairs(COVERAGE_TAG.out.tagged.map { meta, _t, c, _cov -> [meta, c] }),
           pairs(SCAPP.out.scapp), pairs(RGI_BWT.out.genes))
    STATS(sheet, runs, ROLLUP.out.resistome)

    if (params.depth_experiment) {
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
