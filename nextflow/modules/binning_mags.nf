// Stages 6-7: MetaBAT2 binning, CheckM2 + MIMAG tiers, dRep, GTDB-Tk, MAG abundance

include { linkLayout; runsFile } from './utils'

process BINNING {
    tag "${meta.id}"
    conda params.envs.binning
    publishDir path: { "${params.outdir}/06_binning/${meta.id}" }, mode: 'copy', pattern: 'bins/*.txt'

    input:
    tuple val(meta), path(contigs), path(bams), path(bais)

    output:
    tuple val(meta), path('bins'), emit: bins
    path 'bins/*.txt', optional: true

    script:
    def id = meta.id
    """
    jgi_summarize_bam_contig_depths --outputDepth depth.txt ${bams}
    mkdir -p bins
    metabat2 --inFile ${contigs} --abdFile depth.txt --outFile bins/${id}.bin \\
        --minContig ${params.cfg.binning.metabat2_min_contig} --seed ${params.cfg.binning.metabat2_seed} \\
        --numThreads ${task.cpus}
    """

    stub:
    """
    mkdir -p bins; printf '>${meta.id}_k127_1\\nACGT\\n' > bins/${meta.id}.bin.1.fa
    """
}

process CHECKM2_MIMAG {
    tag "${meta.id}"
    conda params.envs.checkm2
    publishDir path: { "${params.outdir}/06_binning/${meta.id}" }, mode: 'copy', pattern: '{mimag.tsv,checkm2/quality_report.tsv}'

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path('mimag.tsv'), emit: mimag
    path 'pass/*.fa', optional: true, emit: pass
    path 'checkm2/quality_report.tsv'

    script:
    def c = params.cfg.binning
    """
    ls ${bins}/*.fa > /dev/null
    checkm2 predict --input ${bins} -x fa --output-directory checkm2 --threads ${task.cpus} --force
    rm -rf checkm2/protein_files checkm2/diamond_output
    PYTHONPATH=${params.repo}/src ${params.envs.core}/bin/python -m plasmid_detection.tables mimag \\
        --report checkm2/quality_report.tsv --run ${meta.id} \\
        --high-comp ${c.mimag_high_completeness} --high-cont ${c.mimag_high_contamination} \\
        --medium-comp ${c.mimag_medium_completeness} --medium-cont ${c.mimag_medium_contamination} \\
        --out mimag.tsv
    mkdir -p pass
    awk -F'\\t' 'NR>1 && (\$5=="high" || \$5=="medium") {print \$2}' mimag.tsv | while read -r b; do
        cp ${bins}/\$b.fa pass/
    done
    """

    stub:
    """
    mkdir -p checkm2 pass
    printf 'Name\\tCompleteness\\tContamination\\n${meta.id}.bin.1\\t95\\t1\\n' > checkm2/quality_report.tsv
    printf 'run\\tbin\\tcompleteness\\tcontamination\\ttier\\n${meta.id}\\t${meta.id}.bin.1\\t95\\t1\\thigh\\n' > mimag.tsv
    cp ${bins}/${meta.id}.bin.1.fa pass/
    """
}

process DREP {
    conda params.envs.drep
    publishDir "${params.outdir}/07_mags", mode: 'copy', pattern: '{genome_info.csv,drep/data_tables,drep/log}'

    input:
    val mimag        // [[run, mimag.tsv], ...]
    path genomes, stageAs: 'pass/*'

    output:
    path 'drep', emit: drep
    path 'genome_info.csv'

    script:
    def runs = mimag.collect { x -> x[0] }
    """
    ${linkLayout(mimag, 'binning/{run}/mimag.tsv')}
    ${runsFile(runs)}
    PYTHONPATH=${params.repo}/src ${params.envs.core}/bin/python -m plasmid_detection.tables drep-info \\
        --binning-dir binning --samples runs.txt --out genome_info.csv
    n=\$(ls pass/*.fa | wc -l)
    if (( n == 1 )); then
        # One MAG: nothing to dereplicate; lay out what dRep would have produced
        mkdir -p drep/dereplicated_genomes drep/data_tables
        cp pass/*.fa drep/dereplicated_genomes/
        g=\$(basename pass/*.fa)
        printf 'genome,secondary_cluster\\n%s,1_1\\n' \$g > drep/data_tables/Cdb.csv
        printf 'genome,cluster\\n%s,1_1\\n' \$g > drep/data_tables/Wdb.csv
    else
        dRep dereplicate drep -g pass/*.fa -p ${task.cpus} -sa ${params.cfg.mags.drep_ani} \\
            --genomeInfo genome_info.csv \\
            -comp ${params.cfg.binning.mimag_medium_completeness} -con ${params.cfg.binning.mimag_medium_contamination}
    fi
    """

    stub:
    """
    mkdir -p drep/dereplicated_genomes drep/data_tables
    cp pass/*.fa drep/dereplicated_genomes/
    printf 'genome,secondary_cluster\\n' > drep/data_tables/Cdb.csv
    printf 'genome,cluster\\n' > drep/data_tables/Wdb.csv
    touch genome_info.csv
    """
}

process GTDBTK {
    conda params.envs.gtdbtk
    publishDir "${params.outdir}/07_mags", mode: 'copy', pattern: 'gtdbtk/*.summary.tsv'

    input:
    path drep

    output:
    path 'gtdbtk', emit: gtdbtk
    path 'gtdbtk/*.summary.tsv'

    script:
    """
    gtdbtk classify_wf --genome_dir ${drep}/dereplicated_genomes -x fa --out_dir gtdbtk \\
        --cpus ${task.cpus} ${params.cfg.mags.gtdbtk_args}
    rm -rf gtdbtk/align gtdbtk/identify
    """

    stub:
    """
    mkdir -p gtdbtk; printf 'user_genome\\tclassification\\n' > gtdbtk/gtdbtk.bac120.summary.tsv
    """
}

process MAG_ABUNDANCE {
    conda params.envs.coverm
    publishDir "${params.outdir}/07_mags", mode: 'copy', pattern: 'coverm_genome.tsv'
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: '07_mags.tsv'

    input:
    path drep
    path gtdbtk
    val mimag        // [[run, mimag.tsv], ...]
    path reads, stageAs: 'reads/*'

    output:
    path 'coverm_genome.tsv'
    path '07_mags.tsv'

    script:
    def runs = mimag.collect { x -> x[0] }
    def coupled = runs.collect { r -> "reads/${r}_1.fastq.gz reads/${r}_2.fastq.gz" }.join(' ')
    """
    coverm genome --coupled ${coupled} --genome-fasta-directory ${drep}/dereplicated_genomes -x fa \\
        --methods ${params.cfg.mags.coverm_genome_methods} --threads ${task.cpus} --output-file coverm_genome.tsv
    ${linkLayout(mimag, 'binning/{run}/mimag.tsv')}
    mkdir -p mags && ln -s ../${drep} mags/drep && ln -s ../${gtdbtk} mags/gtdbtk && ln -s ../coverm_genome.tsv mags/
    ${runsFile(runs)}
    PYTHONPATH=${params.repo}/src ${params.envs.core}/bin/python -m plasmid_detection.tables mags \\
        --samples runs.txt --binning-dir binning --mags-dir mags --out 07_mags.tsv
    """

    stub:
    """
    touch coverm_genome.tsv 07_mags.tsv
    """
}
