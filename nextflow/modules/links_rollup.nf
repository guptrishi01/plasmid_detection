// Stages 8-9: MAG AMR, plasmid-host inference, RGI read level, rollup, statistics,
// opt-in depth experiment

include { linkLayout; runsFile } from './utils'

def PY() {
    "PYTHONPATH=${params.repo}/src ${params.envs.core}/bin/python -m plasmid_detection.tables"
}

process MAG_AMR {
    conda params.envs.amrfinderplus
    publishDir "${params.outdir}/08_mag_links", mode: 'copy', pattern: 'amrfinder'
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: '08_mag_amr.tsv'

    input:
    path drep

    output:
    path 'amrfinder'
    path '08_mag_amr.tsv'

    script:
    """
    mkdir -p amrfinder
    for g in ${drep}/dereplicated_genomes/*.fa; do
        mag=\$(basename \$g .fa)
        amrfinder --nucleotide \$g ${params.cfg.contigs.amrfinder_args} --database ${params.db.amrfinderplus} \\
            --threads ${task.cpus} --name \$mag --output amrfinder/\$mag.tsv
    done
    first=1
    for f in amrfinder/*.tsv; do
        if (( first )); then cat \$f; first=0; else tail -n +2 \$f; fi
    done > 08_mag_amr.tsv
    """

    stub:
    """
    mkdir -p amrfinder; touch 08_mag_amr.tsv
    """
}

process PLASMID_HOSTS {
    conda params.envs.core
    publishDir "${params.outdir}/08_mag_links", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: '08_*.tsv'

    input:
    val bins         // [[run, bins/], ...]
    val mimag        // [[run, mimag.tsv], ...]
    val genomad      // [[run, genomad/], ...]
    val mob          // [[run, mob_typer.tsv], ...] (runs with plasmid contigs only)
    val coverm       // [[run, coverm_contig.tsv], ...]
    path drep
    path gtdbtk

    output:
    path '*.tsv'

    script:
    def runs = bins.collect { x -> x[0] }
    """
    ${linkLayout(bins, 'binning/{run}/bins')}
    ${linkLayout(mimag, 'binning/{run}/mimag.tsv')}
    ${linkLayout(genomad, 'contigs/{run}/genomad')}
    ${linkLayout(mob, 'contigs/{run}/mob_typer.tsv')}
    ${linkLayout(coverm, 'contigs/{run}/coverm_contig.tsv')}
    mkdir -p mags && ln -s ../${drep} mags/drep && ln -s ../${gtdbtk} mags/gtdbtk
    ${runsFile(runs)}
    ${PY()} plasmid-hosts --samples runs.txt --binning-dir binning --mags-dir mags --contigs-dir contigs \\
        --min-samples ${params.cfg.mag_links.min_samples_for_correlation} --out-dir .
    cp plasmids_in_mags.tsv 08_plasmids_in_mags.tsv
    cp plasmid_host_inference.tsv 08_plasmid_host_inference.tsv
    """

    stub:
    """
    touch plasmids_in_mags.tsv plasmid_host_inference.tsv 08_plasmids_in_mags.tsv 08_plasmid_host_inference.tsv
    """
}

process RGI_BWT {
    tag "${meta.id}"
    conda params.envs.rgi
    publishDir path: { "${params.outdir}/09_rollup/rgi_bwt/${meta.id}" }, mode: 'copy', pattern: '*.{txt,json}'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.gene_mapping_data.txt"), emit: genes
    path '*.txt'

    script:
    def id = meta.id
    def idx = "${params.db.card.dir}/localDB/bwt/card_reference"
    """
    # RGI builds its Bowtie2 index of CARD inside the shared localDB on first use; build it
    # once here, under a lock, so parallel tasks don't all write it at once
    (
        flock 9
        if [[ ! -s ${idx}/bowtie2.rev.2.bt2 ]]; then
            mkdir -p ${idx}
            bowtie2-build --quiet --threads ${task.cpus} ${params.db.card.dir}/localDB/card_reference.fasta ${idx}/bowtie2
        fi
    ) 9> ${params.db.card.dir}/localDB/.bwt_index.lock
    ln -sfn ${params.db.card.dir}/localDB localDB
    rgi bwt --read_one ${r1} --read_two ${r2} --aligner bowtie2 --threads ${task.cpus} \\
        --output_file ${id} --local --clean
    rm -f localDB *.bam *.bai *.allele_mapping_data.json
    """

    stub:
    """
    printf 'ARO Term\\tAll Mapped Reads\\n' > ${meta.id}.gene_mapping_data.txt
    """
}

process ROLLUP {
    conda params.envs.core
    publishDir "${params.outdir}/09_rollup", mode: 'copy', pattern: 'resistome_*.tsv'
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: '09_*.tsv'

    input:
    val genomad      // [[run, genomad/], ...]
    val mob          // [[run, mob_typer.tsv], ...]
    val tagged       // [[run, amr_tagged.tsv], ...]
    val counts       // [[run, amr_counts.tsv], ...]
    val scapp        // [[run, scapp/], ...]
    val bwt          // [[run, gene_mapping_data.txt], ...]

    output:
    path 'resistome_*.tsv', emit: resistome
    path '09_*.tsv'

    script:
    def runs = genomad.collect { x -> x[0] }
    def types = [params.cfg.rollup.amrfinder_types].flatten().join(' ')
    """
    ${linkLayout(genomad, 'contigs/{run}/genomad')}
    ${linkLayout(mob, 'contigs/{run}/mob_typer.tsv')}
    ${linkLayout(tagged, 'contigs/{run}/amr_tagged.tsv')}
    ${linkLayout(counts, 'contigs/{run}/amr_counts.tsv')}
    ${linkLayout(scapp, 'asm/{run}/scapp')}
    ${linkLayout(bwt, 'bwt/{run}/{run}.gene_mapping_data.txt')}
    ${runsFile(runs)}
    ${PY()} plasmid-summary --samples runs.txt --contigs-dir contigs --assembly-dir asm \\
        --amrfinder-types ${types} --out 09_plasmid_summary.tsv
    for src in amrfinder rgi; do
        ${PY()} resistome --samples runs.txt --contigs-dir contigs --source \$src \\
            --amrfinder-types ${types} --out-dir .
    done
    ${PY()} rgi-compare --samples runs.txt --contigs-dir contigs --bwt-dir bwt --out 09_rgi_read_vs_contig.tsv
    """

    stub:
    """
    touch resistome_amrfinder_all.tsv resistome_amrfinder_mobile.tsv 09_plasmid_summary.tsv 09_rgi_read_vs_contig.tsv
    """
}

process STATS {
    conda params.envs.stats
    publishDir "${params.outdir}/09_rollup", mode: 'copy', pattern: 'stats'
    publishDir "${params.outdir}/metadata", mode: 'copy', pattern: 'stats_samples.tsv'

    input:
    path sheet
    val runs
    path resistome

    output:
    path 'stats'
    path 'stats_samples.tsv'

    script:
    def s = params.cfg.stats
    def groups = [s.site_groups].flatten().collect { g -> "'${g}'" }.join(' ')
    """
    ${runsFile(runs)}
    # BioSample antibiotic concentrations are fetched from NCBI by stats-samples
    ${PY()} stats-samples --samples runs.txt --sheet ${sheet} --site-groups ${groups} \\
        --nd-encoding ${s.nd_encoding} --loq-encoding ${s.loq_encoding} --out stats_samples.tsv
    for part in all mobile; do
        Rscript ${params.repo}/scripts/R/stats.R stats_samples.tsv resistome_amrfinder_\${part}.tsv stats/\$part
    done
    """

    stub:
    """
    mkdir -p stats; touch stats_samples.tsv
    """
}

// Depth experiment (opt-in): subsample -> MEGAHIT -> geNomad, one (run, level) per task
process DEPTH_SUBSAMPLE {
    tag "${meta.id} ${level}M"
    conda params.envs.assembly
    publishDir path: { "${params.outdir}/09_rollup/depth/${meta.id}_${level}M" }, mode: 'copy', pattern: '{yield.tsv,SKIPPED}'

    input:
    tuple val(meta), path(r1), path(r2), path(counts), val(level)

    output:
    path 'yield.tsv', optional: true, emit: yield
    path 'SKIPPED', optional: true

    script:
    def id = meta.id
    def n = (level as long) * 1000000L
    """
    total=\$(awk -F'\\t' 'NR==2{print \$4}' ${counts})
    if (( ${n} > total )); then
        echo "skipped: ${level}M > \$total available pairs" > SKIPPED
        exit 0
    fi
    p=\$(awk -v n=${n} -v t=\$total 'BEGIN{p=n*1.02/t; if(p>1)p=1; printf "%.6f", p}')
    for r in 1 2; do
        in=\$( [[ \$r == 1 ]] && echo ${r1} || echo ${r2} )
        # seqkit head stops reading after N records, so seqkit sample dies of SIGPIPE (141):
        # expected; counts and ids are checked below
        seqkit sample -p \$p -s ${params.cfg.rollup.depth_seed} \$in 2>/dev/null \\
            | seqkit head -n ${n} -o sub_\${r}.fastq.gz \\
            || [[ \${PIPESTATUS[0]} == 141 && \${PIPESTATUS[1]} == 0 ]]
    done
    n1=\$(seqkit stats -T sub_1.fastq.gz | awk 'NR==2{print \$4}')
    n2=\$(seqkit stats -T sub_2.fastq.gz | awk 'NR==2{print \$4}')
    l1=\$(seqkit seq -n -i sub_1.fastq.gz | tail -1); l2=\$(seqkit seq -n -i sub_2.fastq.gz | tail -1)
    [[ \$n1 == ${n} && \$n2 == ${n} && \$l1 == \$l2 ]] || { echo "subsample out of sync" >&2; exit 1; }

    mem=${(task.memory.toBytes() * 0.9).longValue()}
    megahit -1 sub_1.fastq.gz -2 sub_2.fastq.gz -o megahit -t ${task.cpus} -m \$mem \\
        --min-contig-len ${params.cfg.assembly.min_contig_len} ${params.cfg.assembly.megahit_args}
    awk -v s=${id} '/^>/{sub(/^>/, ">" s "_")} {print}' megahit/final.contigs.fa > contigs.fa
    ${params.envs.genomad}/bin/genomad end-to-end ${params.cfg.contigs.genomad_args} --threads ${task.cpus} \\
        contigs.fa genomad ${params.db.genomad}
    ${PY()} depth-yield --run ${id} --level ${level}M --pairs ${n} --contigs contigs.fa \\
        --genomad-dir genomad --out yield.tsv
    rm -rf sub_*.fastq.gz megahit genomad contigs.fa
    """

    stub:
    """
    printf 'run\\tlevel\\tpairs\\tcontigs\\tassembly_bp\\tplasmid_contigs\\tplasmid_bp\\n${meta.id}\\t${level}M\\t0\\t0\\t0\\t0\\t0\\n' > yield.tsv
    """
}

process DEPTH_SUMMARY {
    conda params.envs.stats
    publishDir "${params.outdir}/09_rollup/depth", mode: 'copy', pattern: 'plasmid_yield_vs_depth.*'
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: '09_depth_experiment.tsv'

    input:
    path yields, stageAs: 'sub/y?/yield.tsv'
    val full         // [[run, contigs.fa, genomad/, counts.tsv], ...] for the full-depth point

    output:
    path '09_depth_experiment.tsv'
    path 'plasmid_yield_vs_depth.*'

    script:
    def fullCmds = full.collect { run, contigs, genomad, counts ->
        "mkdir -p sub/${run}_full && ${PY()} depth-yield --run ${run} --level full " +
        "--pairs \$(awk -F'\\t' 'NR==2{print \$4}' ${counts}) --contigs ${contigs} " +
        "--genomad-dir ${genomad} --out sub/${run}_full/yield.tsv"
    }.join('\n')
    """
    ${fullCmds}
    ${PY()} depth --depth-dir sub --out 09_depth_experiment.tsv
    Rscript ${params.repo}/scripts/R/depth_plot.R 09_depth_experiment.tsv plasmid_yield_vs_depth
    """

    stub:
    """
    touch 09_depth_experiment.tsv plasmid_yield_vs_depth.png
    """
}
