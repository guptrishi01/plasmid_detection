// Stage 1: run list from NCBI, sample sheet, and the SRA download (sra-tools only, per spec)

process FETCH_RUNINFO {
    label 'local'
    conda params.envs.core
    publishDir "${params.outdir}/metadata", mode: 'copy'

    input:
    val id
    val strategy
    val runs

    output:
    path 'runinfo.csv'

    script:
    """
    fetch_runinfo.py '${id}' '${strategy}' runinfo.csv ${runs ? "'${runs}'" : ''}
    """
}

process SAMPLE_SHEET {
    label 'local'
    conda params.envs.core
    publishDir "${params.outdir}/metadata", mode: 'copy'

    input:
    path runinfo
    val regex

    output:
    path 'samples.tsv'

    script:
    """
    sample_sheet.py ${runinfo} '${regex}' samples.tsv
    """
}

// storeDir: a permanent cache. A run whose two FASTQs are already in raw_dir is not
// downloaded again (the 92 runs fetched by 01_download_reads.sbatch are reused as-is).
process DOWNLOAD_READS {
    tag "${meta.id}"
    conda params.envs.data
    storeDir params.raw_dir

    input:
    val meta

    output:
    tuple val(meta), path("${meta.id}_1.fastq.gz"), path("${meta.id}_2.fastq.gz")

    script:
    def run = meta.id
    """
    # Per-task sra-tools settings (never ~/.ncbi): full-quality .sra, not SRA Lite
    export NCBI_SETTINGS=\$PWD/.ncbi/user-settings.mkfg
    mkdir -p .ncbi tmp fastq
    vdb-config --simplified-quality-scores no > /dev/null

    prefetch --max-size u -O . ${run}
    [[ -f ${run}/${run}.sra ]] || { echo "expected full-quality ${run}.sra, got: \$(ls ${run})" >&2; exit 1; }
    vdb-validate ${run}/${run}.sra

    fasterq-dump --split-3 --skip-technical --threads ${task.cpus} --temp tmp -O fastq ${run}/${run}.sra
    [[ "\$(ls fastq)" == "${run}_1.fastq"\$'\\n'"${run}_2.fastq" ]] || { echo "unexpected fasterq-dump output: \$(ls fastq)" >&2; exit 1; }

    nqual=\$(head -n 40000 fastq/${run}_1.fastq | awk 'NR%4==0' | fold -w1 | sort -u | wc -l)
    (( nqual >= 3 )) || { echo "only \$nqual distinct quality characters: simplified qualities?" >&2; exit 1; }

    n1=\$(awk 'END{print NR/4}' fastq/${run}_1.fastq)
    n2=\$(awk 'END{print NR/4}' fastq/${run}_2.fastq)
    status=PASS
    [[ "\$n1" == "\$n2" && "\$n1" == "${meta.spots}" ]] || status=FAIL

    pigz -p ${task.cpus} fastq/${run}_1.fastq fastq/${run}_2.fastq
    mv fastq/${run}_1.fastq.gz fastq/${run}_2.fastq.gz .
    md5_1=\$(md5sum ${run}_1.fastq.gz | cut -d' ' -f1)
    md5_2=\$(md5sum ${run}_2.fastq.gz | cut -d' ' -f1)
    rm -rf ${run} fastq tmp

    mkdir -p ${params.manifest_dir}
    printf 'run\\tsra_spots\\tpairs_r1\\tpairs_r2\\tdistinct_quals\\tr1\\tr1_md5\\tr2\\tr2_md5\\tstatus\\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        ${run} ${meta.spots} \$n1 \$n2 \$nqual ${run}_1.fastq.gz \$md5_1 ${run}_2.fastq.gz \$md5_2 \$status \\
        > ${params.manifest_dir}/${run}.tsv
    [[ \$status == PASS ]] || { cat ${params.manifest_dir}/${run}.tsv >&2; exit 1; }
    chmod a-w ${run}_1.fastq.gz ${run}_2.fastq.gz
    """

    stub:
    """
    echo "stub: ${meta.id} not in raw_dir; refusing to fake reads" >&2; exit 1
    """
}
