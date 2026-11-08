// Stages 6-7: MetaBAT2 binning, CheckM2 + MIMAG tiers, dRep per binning group then across
// groups (catalog), GTDB-Tk, MAG abundance

include { linkLayout; runsFile; dereplicate } from './utils'

process BINNING {
    tag "${meta.id}"
    conda params.envs.binning
    publishDir path: { "${params.outdir}/06_binning/${meta.id}" }, mode: 'copy', pattern: 'bins/*.txt'

    input:
    tuple val(meta), path(contigs), path(depth)    // depth.txt over the group BAMs (MAP_GROUP)

    output:
    tuple val(meta), path('bins'), emit: bins
    path 'bins/*.txt', optional: true

    script:
    def id = meta.id
    """
    mkdir -p bins
    metabat2 --inFile ${contigs} --abdFile ${depth} --outFile bins/${id}.bin \\
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
    tuple val(meta), path('pass'), emit: pass       // medium+ bins; may be empty
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

// dRep --genomeInfo from every run's MIMAG table (CheckM2 completeness / contamination)
process GENOME_INFO {
    conda params.envs.core
    publishDir "${params.outdir}/07_mags", mode: 'copy'

    input:
    val mimag        // [[run, mimag.tsv], ...]

    output:
    path 'genome_info.csv'

    script:
    """
    ${linkLayout(mimag, 'binning/{run}/mimag.tsv')}
    ${runsFile(mimag.collect { x -> x[0] })}
    PYTHONPATH=${params.repo}/src python -m plasmid_detection.tables drep-info \\
        --binning-dir binning --samples runs.txt --out genome_info.csv
    """

    stub:
    """
    touch genome_info.csv
    """
}

// Level 1: one MAG per species within a binning group (no output if the group has none)
process DREP_GROUP {
    tag "${group}"
    conda params.envs.drep
    // Patterns inside an emitted dir aren't published, so the tables are copied out first
    publishDir path: { "${params.outdir}/07_mags/groups/${group}" }, mode: 'copy', pattern: 'drep_data_tables',
        saveAs: { _f -> 'drep/data_tables' }

    input:
    tuple val(group), path(pass, stageAs: 'pass*')   // the group's runs' pass/ dirs
    path genome_info

    output:
    tuple val(group), path('drep'), emit: drep, optional: true
    path 'drep_data_tables', optional: true

    script:
    """
    mkdir -p genomes
    for d in pass*; do
        for f in \$d/*.fa; do [[ -e \$f ]] && ln -s "\$(readlink -f \$f)" genomes/; done
    done
    if [[ -z "\$(ls genomes)" ]]; then echo "group ${group}: no medium+ MAGs"; exit 0; fi
    ${dereplicate('drep', 'genomes/*.fa', task.cpus)}
    cp -r drep/data_tables drep_data_tables
    """

    stub:
    """
    mkdir -p drep/dereplicated_genomes drep/data_tables
    for d in pass*; do cp \$d/*.fa drep/dereplicated_genomes/ 2>/dev/null || true; done
    printf 'genome,secondary_cluster\\n' > drep/data_tables/Cdb.csv
    printf 'genome,cluster\\n' > drep/data_tables/Wdb.csv
    cp -r drep/data_tables drep_data_tables
    """
}

// Level 2: the catalog, one MAG per species across all groups (GTDB-Tk, CoverM, stage 8)
process DREP_CATALOG {
    conda params.envs.drep
    publishDir "${params.outdir}/07_mags", mode: 'copy', pattern: 'drep_data_tables', saveAs: { _f -> 'drep/data_tables' }

    input:
    val groups       // [[group, drep/], ...]
    path genome_info

    output:
    path 'drep', emit: drep
    path 'drep_data_tables'

    script:
    """
    ${linkLayout(groups, 'groups/{run}/drep')}
    ${dereplicate('drep', 'groups/*/drep/dereplicated_genomes/*.fa', task.cpus)}
    cp -r drep/data_tables drep_data_tables
    """

    stub:
    """
    mkdir -p drep/dereplicated_genomes drep/data_tables
    cp -L ${groups.collect { g -> "${g[1]}/dereplicated_genomes/*.fa" }.join(' ')} drep/dereplicated_genomes/ 2>/dev/null || true
    printf 'genome,secondary_cluster\\n' > drep/data_tables/Cdb.csv
    printf 'genome,cluster\\n' > drep/data_tables/Wdb.csv
    cp -r drep/data_tables drep_data_tables
    """
}

process GTDBTK {
    conda params.envs.gtdbtk
    publishDir "${params.outdir}/07_mags/gtdbtk", mode: 'copy', pattern: 'gtdbtk.*.summary.tsv'

    input:
    path drep

    output:
    path 'gtdbtk', emit: gtdbtk
    path 'gtdbtk.*.summary.tsv'      // copies: a pattern inside an emitted dir isn't published

    script:
    """
    gtdbtk classify_wf --genome_dir ${drep}/dereplicated_genomes -x fa --out_dir gtdbtk \\
        --cpus ${task.cpus} ${params.cfg.mags.gtdbtk_args}
    rm -rf gtdbtk/align gtdbtk/identify
    cp gtdbtk/gtdbtk.*.summary.tsv .
    """

    stub:
    """
    mkdir -p gtdbtk; printf 'user_genome\\tclassification\\n' > gtdbtk/gtdbtk.bac120.summary.tsv
    cp gtdbtk/gtdbtk.*.summary.tsv .
    """
}

process MAG_ABUNDANCE {
    conda params.envs.coverm
    publishDir "${params.outdir}/07_mags", mode: 'copy', pattern: 'coverm_genome.tsv'
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: '07_*mags.tsv'

    input:
    path drep
    path gtdbtk
    val mimag        // [[run, mimag.tsv], ...]
    val groups       // [[group, drep/], ...] (level 1)
    path reads, stageAs: 'reads/*'

    output:
    path 'coverm_genome.tsv'
    path '07_mags.tsv'
    path '07_group_mags.tsv'

    script:
    def runs = mimag.collect { x -> x[0] }
    def coupled = runs.collect { r -> "reads/${r}_1.fastq.gz reads/${r}_2.fastq.gz" }.join(' ')
    """
    coverm genome --coupled ${coupled} --genome-fasta-directory ${drep}/dereplicated_genomes -x fa \\
        --methods ${params.cfg.mags.coverm_genome_methods} --threads ${task.cpus} --output-file coverm_genome.tsv
    ${linkLayout(mimag, 'binning/{run}/mimag.tsv')}
    mkdir -p mags && ln -s ../${drep} mags/drep && ln -s ../${gtdbtk} mags/gtdbtk && ln -s ../coverm_genome.tsv mags/
    ${linkLayout(groups, 'mags/groups/{run}/drep')}
    ${runsFile(runs)}
    PYTHONPATH=${params.repo}/src ${params.envs.core}/bin/python -m plasmid_detection.tables mags \\
        --samples runs.txt --binning-dir binning --mags-dir mags --out 07_mags.tsv --group-out 07_group_mags.tsv
    """

    stub:
    """
    touch coverm_genome.tsv 07_mags.tsv 07_group_mags.tsv
    """
}
