// Stages 2-3: read QC (FastQC, fastp, FastQC, MultiQC) and human-read removal (Bowtie2)

include { linkLayout; runsFile } from './utils'

process FASTQC_RAW {
    tag "${meta.id}"
    conda params.envs.qc
    publishDir "${params.outdir}/02_qc/fastqc_raw", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    path '*_fastqc.{zip,html}'

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${r1} ${r2}
    """

    stub:
    """
    touch ${meta.id}_1_fastqc.zip ${meta.id}_1_fastqc.html
    """
}

process FASTP {
    tag "${meta.id}"
    conda params.envs.qc
    publishDir "${params.outdir}/02_qc/fastp", mode: 'copy', pattern: '*.fastp.*'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("trimmed/${meta.id}_1.fastq.gz"), path("trimmed/${meta.id}_2.fastq.gz"), emit: reads
    path "${meta.id}.fastp.json", emit: json
    path "${meta.id}.fastp.{html,log}"

    script:
    def id = meta.id
    """
    mkdir -p trimmed
    fastp -i ${r1} -I ${r2} -o trimmed/${id}_1.fastq.gz -O trimmed/${id}_2.fastq.gz \\
        --thread ${task.cpus} ${params.cfg.qc.fastp_args} \\
        --json ${id}.fastp.json --html ${id}.fastp.html --report_title ${id} 2> ${id}.fastp.log
    """

    stub:
    """
    mkdir -p trimmed; echo | gzip > trimmed/${meta.id}_1.fastq.gz; echo | gzip > trimmed/${meta.id}_2.fastq.gz; touch ${meta.id}.fastp.json ${meta.id}.fastp.html ${meta.id}.fastp.log
    """
}

process FASTQC_TRIMMED {
    tag "${meta.id}"
    conda params.envs.qc
    publishDir "${params.outdir}/02_qc/fastqc_trimmed", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    path '*_fastqc.{zip,html}'

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${r1} ${r2}
    """

    stub:
    """
    touch ${meta.id}_1_fastqc.zip
    """
}

process QC_SUMMARY {
    conda params.envs.qc
    publishDir "${params.outdir}/02_qc/multiqc", mode: 'copy', pattern: 'read_qc*'
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: '02_read_qc.tsv'

    input:
    val runs
    path raw, stageAs: 'fastqc_raw/*'
    path fastp, stageAs: 'fastp/*'
    path trimmed, stageAs: 'fastqc_trimmed/*'

    output:
    path 'read_qc*'
    path '02_read_qc.tsv'

    script:
    """
    # --dirs: raw and trimmed FastQC reports share sample names
    multiqc --force --dirs --dirs-depth 1 --filename read_qc fastqc_raw fastp fastqc_trimmed
    ${runsFile(runs)}
    PYTHONPATH=${params.repo}/src ${params.envs.core}/bin/python -m plasmid_detection.tables qc \\
        --samples runs.txt --fastp-dir fastp --out 02_read_qc.tsv
    """

    stub:
    """
    touch read_qc.html 02_read_qc.tsv
    """
}

process HOST_DEPLETION {
    tag "${meta.id}"
    conda params.envs.host
    publishDir "${params.outdir}/03_host", mode: 'copy', pattern: '*.{counts.tsv,bowtie2.log}'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("hostfree/${meta.id}_1.fastq.gz"), path("hostfree/${meta.id}_2.fastq.gz"), emit: reads
    path "${meta.id}.bowtie2.log", emit: log
    path "${meta.id}.counts.tsv", emit: counts

    script:
    def id = meta.id
    """
    mkdir -p hostfree
    # Keep pairs where neither mate maps (-f 12), primary records only (-F 256)
    bowtie2 -p ${task.cpus} ${params.cfg.host.bowtie2_args} -x ${params.reference.human.bowtie2_index} \\
            -1 ${r1} -2 ${r2} 2> ${id}.bowtie2.log \\
        | samtools view -u -f 12 -F 256 - \\
        | samtools fastq -@ 4 -n -c 6 -1 hostfree/${id}_1.fastq.gz -2 hostfree/${id}_2.fastq.gz -0 /dev/null -s /dev/null -

    trimmed=\$(awk '/reads; of these:/{print \$1; exit}' ${id}.bowtie2.log)
    n1=\$(zcat hostfree/${id}_1.fastq.gz | awk 'END{print NR/4}')
    n2=\$(zcat hostfree/${id}_2.fastq.gz | awk 'END{print NR/4}')
    [[ "\$n1" == "\$n2" ]] || { echo "R1/R2 pair counts differ: \$n1 vs \$n2" >&2; exit 1; }
    printf 'run\\ttrimmed_pairs\\thuman_pairs\\tretained_pairs\\n%s\\t%s\\t%s\\t%s\\n' \\
        ${id} \$trimmed \$((trimmed - n1)) \$n1 > ${id}.counts.tsv
    """

    stub:
    """
    mkdir -p hostfree; echo | gzip > hostfree/${meta.id}_1.fastq.gz; echo | gzip > hostfree/${meta.id}_2.fastq.gz; touch ${meta.id}.bowtie2.log; printf "run\\ttrimmed_pairs\\thuman_pairs\\tretained_pairs\\n${meta.id}\\t10\\t0\\t10\\n" > ${meta.id}.counts.tsv
    """
}

process FASTQC_HOSTFREE {
    tag "${meta.id}"
    conda params.envs.qc
    publishDir "${params.outdir}/03_host/fastqc_hostfree", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    path '*_fastqc.{zip,html}'

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${r1} ${r2}
    """

    stub:
    """
    touch ${meta.id}_1_fastqc.zip
    """
}

process HOST_SUMMARY {
    conda params.envs.qc
    publishDir "${params.outdir}/03_host/multiqc", mode: 'copy', pattern: '{host_depletion,hostfree_qc}*'
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: '03_host_depletion.tsv'

    input:
    val runs
    path logs, stageAs: 'bowtie2_logs/*'
    path counts, stageAs: 'counts/*'
    path fastp, stageAs: 'fastp/*'
    path fastqc, stageAs: 'fastqc_hostfree/*'

    output:
    path 'host_depletion*'
    path 'hostfree_qc*'
    path '03_host_depletion.tsv'

    script:
    """
    multiqc --force --filename host_depletion bowtie2_logs
    multiqc --force --dirs --dirs-depth 1 --filename hostfree_qc bowtie2_logs fastqc_hostfree
    ${runsFile(runs)}
    PYTHONPATH=${params.repo}/src ${params.envs.core}/bin/python -m plasmid_detection.tables host \\
        --samples runs.txt --fastp-dir fastp --counts-dir counts --out 03_host_depletion.tsv
    """

    stub:
    """
    touch host_depletion.html hostfree_qc.html 03_host_depletion.tsv
    """
}
