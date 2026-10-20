// Stage 4: MEGAHIT per run, one metaSPAdes attempt, within-system mapping, QUAST +
// Nonpareil, SCAPP, assembly table

include { linkLayout; runsFile } from './utils'

process MEGAHIT {
    tag "${meta.id}"
    conda params.envs.assembly
    publishDir path: { "${params.outdir}/04_assembly/${meta.id}" }, mode: 'copy', pattern: '{megahit/log,megahit/options.json,graph_k.txt}'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.contigs.fa"), emit: contigs
    tuple val(meta), path("${meta.id}.k*.fastg"), path('graph_k.txt'), emit: graph
    path 'megahit/{log,options.json}'
    path 'graph_k.txt'

    script:
    def id = meta.id
    def mem = (task.memory.toBytes() * 0.9).longValue()   // MEGAHIT sizes itself from node RAM
    """
    megahit -1 ${r1} -2 ${r2} -o megahit -t ${task.cpus} -m ${mem} \\
        --min-contig-len ${params.cfg.assembly.min_contig_len} ${params.cfg.assembly.megahit_args}
    # Sample-prefixed headers: contigs from different runs never share a name
    awk -v s=${id} '/^>/{sub(/^>/, ">" s "_")} {print}' megahit/final.contigs.fa > ${id}.contigs.fa
    # Assembly graph from the largest-k contig set
    k=\$(ls megahit/intermediate_contigs | sed -n 's/^k\\([0-9]*\\)\\.contigs\\.fa\$/\\1/p' | sort -n | tail -1)
    [[ -n "\$k" ]]
    megahit_toolkit contig2fastg \$k megahit/intermediate_contigs/k\$k.contigs.fa > ${id}.k\$k.fastg
    echo \$k > graph_k.txt
    rm -rf megahit/intermediate_contigs
    """

    stub:
    """
    mkdir -p megahit; touch megahit/log megahit/options.json; printf ">${meta.id}_k127_1\\nACGT\\n" > ${meta.id}.contigs.fa; touch ${meta.id}.k127.fastg; echo 127 > graph_k.txt
    """
}

// One attempt, outcome recorded whether it succeeds or not (errorStrategy 'ignore')
process METASPADES {
    tag "${meta.id}"
    conda params.envs.assembly
    publishDir "${params.outdir}/04_assembly/metaspades", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    path 'attempt.tsv'
    path "${meta.id}.time.txt"

    script:
    def id = meta.id
    """
    set +e
    start=\$(date +%s)
    /usr/bin/time -v -o ${id}.time.txt spades.py --meta -1 ${r1} -2 ${r2} -o spades -t ${task.cpus} \\
        -m ${params.cfg.assembly.metaspades_memory_gb}
    rc=\$?
    set -e
    elapsed=\$(( \$(date +%s) - start ))
    rss=\$(awk -F': ' '/Maximum resident set size/{printf "%.1f", \$2/1024/1024}' ${id}.time.txt)
    contigs=NA; ge=NA
    if [[ -s spades/contigs.fasta ]]; then
        contigs=\$(grep -c '^>' spades/contigs.fasta)
        ge=\$(awk -v m=${params.cfg.assembly.min_contig_len} '/^>/{if(l>=m)n++; l=0; next}{l+=length(\$0)} END{if(l>=m)n++; print n+0}' spades/contigs.fasta)
    fi
    err=\$(grep -E "== Error|ERROR|Exception|Killed|not enough memory|finished abnormally" spades/spades.log 2>/dev/null | tail -1 || true)   # no match is the normal case
    printf 'run\\tspots\\texit_code\\telapsed_h\\tmax_rss_gb\\tcontigs\\tcontigs_ge_min_len\\tlast_error_line\\n' > attempt.tsv
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' ${id} ${meta.spots} \$rc \$(awk -v e=\$elapsed 'BEGIN{printf "%.2f", e/3600}') \\
        "\$rss" \$contigs \$ge "\${err:-none}" >> attempt.tsv
    rm -rf spades
    """

    stub:
    """
    touch attempt.tsv ${meta.id}.time.txt
    """
}

// One assembly + the host-free reads of every run in the same system (sample sheet location)
process MAP_SYSTEM {
    tag "${meta.id}"
    conda params.envs.assembly
    publishDir path: { "${params.outdir}/04_assembly/${meta.id}/map" }, mode: 'copy', pattern: '*.bowtie2.log'

    input:
    tuple val(meta), path(contigs), val(mate_ids), path(r1s, stageAs: 'r1/*'), path(r2s, stageAs: 'r2/*')

    output:
    tuple val(meta), path('*.bam'), path('*.bam.bai'), emit: bams
    tuple val(meta), path("${meta.id}.bowtie2.log"), emit: own_log
    path '*.bowtie2.log'

    script:
    """
    mkdir -p idx
    bowtie2-build --threads ${task.cpus} ${contigs} idx/asm > build.log 2>&1
    for r in ${mate_ids.join(' ')}; do
        bowtie2 -p ${task.cpus} ${params.cfg.assembly.bowtie2_args} -x idx/asm -1 r1/\${r}_1.fastq.gz -2 r2/\${r}_2.fastq.gz \\
                2> \$r.bowtie2.log \\
            | samtools sort -@ 4 -m 2G -o \$r.bam -
        samtools index \$r.bam
    done
    rm -rf idx
    """

    stub:
    """
    for r in ${mate_ids.join(" ")}; do touch \$r.bam \$r.bam.bai \$r.bowtie2.log; done
    """
}

process ASSEMBLY_QC {
    tag "${meta.id}"
    conda params.envs.assembly
    publishDir path: { "${params.outdir}/04_assembly/${meta.id}" }, mode: 'copy', pattern: '{quast,nonpareil}'

    input:
    tuple val(meta), path(contigs), path(r1), path(r2)

    output:
    tuple val(meta), path('quast/report.tsv'), emit: report
    path 'quast'
    path 'nonpareil'

    script:
    def id = meta.id
    """
    quast.py ${contigs} -o quast -t ${task.cpus} --min-contig ${params.cfg.assembly.min_contig_len} --labels ${id}
    mkdir -p nonpareil
    zcat ${r1} > ${id}_1.fastq        # Nonpareil reads one uncompressed file
    nonpareil -s ${id}_1.fastq -f fastq ${params.cfg.assembly.nonpareil_args} -b nonpareil/${id} -t ${task.cpus}
    rm ${id}_1.fastq
    """

    stub:
    """
    mkdir -p quast nonpareil; touch quast/report.tsv
    """
}

process SCAPP {
    tag "${meta.id}"
    conda params.envs.scapp
    publishDir path: { "${params.outdir}/04_assembly/${meta.id}" }, mode: 'copy', pattern: 'scapp'

    input:
    tuple val(meta), path(fastg), path(graph_k), path(r1), path(r2)

    output:
    tuple val(meta), path('scapp'), emit: scapp

    script:
    """
    k=\$(cat ${graph_k})
    if ! scapp -g ${fastg} -k \$k -r1 ${r1} -r2 ${r2} -o scapp -p ${task.cpus} ${params.cfg.assembly.scapp_args ?: ''}; then
        # SCAPP 0.1.4 crashes in makeblastdb on an empty *.cycs.fasta: no circular plasmids
        cycs=\$(ls scapp/intermediate_files/*.cycs.fasta 2>/dev/null | head -1 || true)
        [[ -n "\$cycs" && ! -s "\$cycs" ]] || { echo "SCAPP failed" >&2; exit 1; }
        echo "SCAPP: no candidate cycles in the graph -> 0 circular plasmids" > scapp/NO_CYCLES
    fi
    rm -rf scapp/intermediate_files
    """

    stub:
    """
    mkdir -p scapp
    """
}

process ASSEMBLY_TABLE {
    conda params.envs.core
    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    val reports      // [[run, quast/report.tsv], ...]
    val logs         // [[run, own bowtie2 log], ...]
    val graph_ks     // [[run, graph_k.txt], ...]

    output:
    path '04_assembly.tsv'

    script:
    def runs = reports.collect { x -> x[0] }
    """
    ${linkLayout(reports, 'asm/{run}/quast/report.tsv')}
    ${linkLayout(logs.collect { r, f -> [r, f] }, 'asm/{run}/map/{run}.bowtie2.log')}
    ${linkLayout(graph_ks, 'asm/{run}/graph_k.txt')}
    ${runsFile(runs)}
    PYTHONPATH=${params.repo}/src python -m plasmid_detection.tables assembly --samples runs.txt \\
        --assembly-dir asm --out 04_assembly.tsv
    """

    stub:
    """
    touch 04_assembly.tsv
    """
}
