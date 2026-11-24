// Co-assembly exploration (PI, 2026-11-19; nextflow/coassembly.nf): one MEGAHIT co-assembly
// per group of sites (params.cfg.coassembly.groups), binned with every run of the group as a
// coverage sample; each run is then mapped on its own to its group's MAGs. The binning,
// CheckM2, dRep, GTDB-Tk, AMR and plasmid-host steps are the 11-group pipeline's processes,
// with the group in the place of a run.

include { linkLayout } from './utils'

def PYCORE() {
    "PYTHONPATH=${params.repo}/src ${params.envs.core}/bin/python"
}

// Every run of the group in one assembly; contigs prefixed with the group name
process COASM_MEGAHIT {
    tag "${group}"
    conda params.envs.assembly
    publishDir path: { "${params.outdir}/04_assembly/${group}" }, mode: 'copy', pattern: '{megahit/log,megahit/options.json,assembly_stats.tsv,runs.txt}'

    input:
    tuple val(group), val(ids), path(r1s, stageAs: 'r1/*'), path(r2s, stageAs: 'r2/*')

    output:
    tuple val(group), path("${group}.contigs.fa"), emit: contigs
    path 'megahit/{log,options.json}'
    path 'assembly_stats.tsv'
    path 'runs.txt'

    script:
    def runs = ids.sort(false)                            // fixed order: same command on -resume
    def mem = (task.memory.toBytes() * 0.9).longValue()   // MEGAHIT sizes itself from node RAM
    """
    printf '%s\\n' ${runs.join(' ')} > runs.txt
    megahit -1 ${runs.collect { r -> "r1/${r}_1.fastq.gz" }.join(',')} \\
            -2 ${runs.collect { r -> "r2/${r}_2.fastq.gz" }.join(',')} \\
        -o megahit -t ${task.cpus} -m ${mem} \\
        --min-contig-len ${params.cfg.assembly.min_contig_len} ${params.cfg.coassembly.megahit_args}
    awk -v s=${group} '/^>/{sub(/^>/, ">" s "_")} {print}' megahit/final.contigs.fa > ${group}.contigs.fa
    rm -rf megahit/intermediate_contigs
    seqkit stats -a -T ${group}.contigs.fa > assembly_stats.tsv
    """

    stub:
    """
    mkdir -p megahit; touch megahit/log megahit/options.json
    printf '%s\\n' ${ids.sort(false).join(' ')} > runs.txt
    printf ">${group}_k141_1\\nACGT\\n" > ${group}.contigs.fa; touch assembly_stats.tsv
    """
}

process COASM_INDEX {
    tag "${group}"
    conda params.envs.assembly

    input:
    tuple val(group), path(contigs)

    output:
    tuple val(group), path('idx'), emit: idx

    script:
    """
    mkdir -p idx
    bowtie2-build --threads ${task.cpus} ${contigs} idx/asm > build.log 2>&1
    """

    stub:
    """
    mkdir -p idx; touch idx/asm.1.bt2
    """
}

// One run's own mapping to its group's co-assembly (the binning coverage sample)
process COASM_MAP {
    tag "${group}:${meta.id}"
    conda params.envs.assembly
    publishDir path: { "${params.outdir}/04_assembly/${group}/map" }, mode: 'copy', pattern: '*.bowtie2.log'

    input:
    tuple val(group), val(meta), path(r1), path(r2), path(idx)

    output:
    tuple val(group), val(meta.id), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    path "${meta.id}.bowtie2.log"

    script:
    """
    bowtie2 -p ${task.cpus} ${params.cfg.assembly.bowtie2_args} -x ${idx}/asm -1 ${r1} -2 ${r2} \\
            2> ${meta.id}.bowtie2.log \\
        | samtools sort -@ 4 -m 2G -o ${meta.id}.bam.tmp -
    mv ${meta.id}.bam.tmp ${meta.id}.bam
    samtools index ${meta.id}.bam
    """

    stub:
    """
    touch ${meta.id}.bam ${meta.id}.bam.bai ${meta.id}.bowtie2.log
    """
}

// Depth over the group's BAMs: MetaBAT2 input and CoverM contig coverage (plasmid hosts)
process COASM_DEPTH {
    tag "${group}"
    conda params.envs.binning
    publishDir path: { "${params.outdir}/06_binning/${group}" }, mode: 'copy', pattern: 'depth.txt'
    publishDir path: { "${params.outdir}/05_contigs/${group}" }, mode: 'copy', pattern: 'coverm_contig.tsv'

    input:
    tuple val(group), path(contigs), val(ids), path(bams, stageAs: 'bam/*'), path(bais, stageAs: 'bam/*')

    output:
    tuple val(group), path('depth.txt'), emit: depth
    tuple val(group), path('coverm_contig.tsv'), emit: coverm

    script:
    def b = ids.sort(false).collect { r -> "bam/${r}.bam" }.join(' ')   // fixed column order
    """
    jgi_summarize_bam_contig_depths --outputDepth depth.txt ${b}
    PATH=${params.envs.coverm}/bin:\$PATH coverm contig --bam-files ${b} --methods ${params.cfg.contigs.coverm_methods} \\
        --threads ${task.cpus} --output-file coverm_contig.tsv
    """

    stub:
    """
    touch depth.txt coverm_contig.tsv
    """
}

// All groups' MAGs in one dRep-shaped dir for GTDB-Tk, MAG AMR and plasmid hosts. Not
// dereplicated across groups (each group keeps its own MAG set): every genome is its own
// cluster, so tables.catalog_links maps bin -> group MAG -> the same group MAG.
process COASM_COLLECT {
    conda params.envs.core
    publishDir "${params.outdir}/07_mags", mode: 'copy', pattern: 'group_mags.tsv'

    input:
    val groups       // [[group, drep/], ...]

    output:
    path 'group_mags', emit: mags
    path 'group_mags.tsv'

    script:
    """
    mkdir -p group_mags/dereplicated_genomes group_mags/data_tables
    printf 'group\\tmag\\n' > group_mags.tsv
    ${groups.collect { g, d -> "for f in ${d}/dereplicated_genomes/*.fa; do cp -L \$f group_mags/dereplicated_genomes/; printf '${g}\\t%s\\n' \$(basename \$f .fa) >> group_mags.tsv; done" }.join('\n    ')}
    n=\$(ls group_mags/dereplicated_genomes | wc -l); m=\$(( \$(wc -l < group_mags.tsv) - 1 ))
    (( n == m )) || { echo "genome name collision across groups: \$n files, \$m MAGs" >&2; exit 1; }
    (( n > 0 )) || { echo "no group MAGs (no medium+ bins in any group)" >&2; exit 1; }
    ( printf 'genome,secondary_cluster\\n'; ls group_mags/dereplicated_genomes | awk '{print \$0 "," NR "_0"}' ) > group_mags/data_tables/Cdb.csv
    ( printf 'genome,cluster\\n'; ls group_mags/dereplicated_genomes | awk '{print \$0 "," NR "_0"}' ) > group_mags/data_tables/Wdb.csv
    """

    stub:
    """
    mkdir -p group_mags/dereplicated_genomes group_mags/data_tables
    printf 'group\\tmag\\n' > group_mags.tsv
    ${groups.collect { g, d -> "for f in ${d}/dereplicated_genomes/*.fa; do cp -L \$f group_mags/dereplicated_genomes/; printf '${g}\\t%s\\n' \$(basename \$f .fa) >> group_mags.tsv; done" }.join('\n    ')}
    ( printf 'genome,secondary_cluster\\n'; ls group_mags/dereplicated_genomes | awk '{print \$0 "," NR "_0"}' ) > group_mags/data_tables/Cdb.csv
    ( printf 'genome,cluster\\n'; ls group_mags/dereplicated_genomes | awk '{print \$0 "," NR "_0"}' ) > group_mags/data_tables/Wdb.csv
    """
}

// Each run of the group mapped on its own to the group's MAGs (CoverM maps every read pair
// set separately): relative abundance and mean coverage per (run, MAG)
process COASM_ABUNDANCE {
    tag "${group}"
    conda params.envs.coverm
    publishDir path: { "${params.outdir}/07_mags/groups/${group}" }, mode: 'copy'

    input:
    tuple val(group), path(drep), val(ids), path(reads, stageAs: 'reads/*')

    output:
    tuple val(group), path('coverm_genome.tsv'), emit: coverm

    script:
    def coupled = ids.sort(false).collect { r -> "reads/${r}_1.fastq.gz reads/${r}_2.fastq.gz" }.join(' ')
    """
    coverm genome --coupled ${coupled} --genome-fasta-directory ${drep}/dereplicated_genomes -x fa \\
        --methods ${params.cfg.mags.coverm_genome_methods} ${params.cfg.coassembly.coverm_read_filters} \\
        --threads ${task.cpus} --output-file coverm_genome.tsv
    """

    stub:
    """
    touch coverm_genome.tsv
    """
}

// Group MAGs vs the per-site MAGs (level-1 dRep of the 11 site groups, prep_inputs.sbatch):
// dRep compare puts all of them into species clusters at params.cfg.coassembly.compare_ani
process COASM_COMPARE {
    conda params.envs.drep
    publishDir "${params.outdir}/07_mags/compare", mode: 'copy'

    input:
    path mags        // COASM_COLLECT group_mags/
    path site_mags   // <coasm_inputs>/site_mags/<site group>/*.fa

    output:
    path 'data_tables', emit: tables
    path 'genomes.tsv'

    script:
    """
    mkdir -p genomes
    printf 'genome\\tset\\tgroup\\n' > genomes.tsv
    for f in ${mags}/dereplicated_genomes/*.fa; do
        ln -s "\$(readlink -f \$f)" genomes/; printf '%s\\tcoassembly\\t\\n' \$(basename \$f) >> genomes.tsv
    done
    for d in ${site_mags}/*/; do
        g=\$(basename \$d)
        for f in \$d*.fa; do ln -s "\$(readlink -f \$f)" genomes/; printf '%s\\tsite\\t%s\\n' \$(basename \$f) \$g >> genomes.tsv; done
    done
    n=\$(ls genomes | wc -l); m=\$(( \$(wc -l < genomes.tsv) - 1 ))
    (( n == m )) || { echo "genome name collision: \$n links, \$m genomes" >&2; exit 1; }
    dRep compare cmp -g genomes/*.fa -p ${task.cpus} -sa ${params.cfg.coassembly.compare_ani}
    cp -r cmp/data_tables data_tables
    """

    stub:
    """
    mkdir -p data_tables; printf 'genome,secondary_cluster\\n' > data_tables/Cdb.csv; touch genomes.tsv
    """
}

// Summary tables: group MAGs, what each run hits, and the comparison with per-site MAGs
process COASM_TABLES {
    conda params.envs.core
    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    path sheet
    val mimag        // [[group, mimag.tsv], ...]
    val coverm       // [[group, coverm_genome.tsv], ...]
    path mag_list    // group_mags.tsv
    path gtdbtk
    path mag_amr     // 08_mag_amr.tsv
    path compare     // dRep compare data_tables
    path genomes     // genomes.tsv (COASM_COMPARE)
    path site_mags_table, stageAs: 'site/07_group_mags.tsv'   // 11-group run
    path site_mag_amr, stageAs: 'site/08_mag_amr.tsv'         // 11-group run (catalog MAGs)

    output:
    path 'coasm_*.tsv'

    script:
    """
    ${linkLayout(mimag, 'mimag/{run}.tsv')}
    ${linkLayout(coverm, 'coverm/{run}.tsv')}
    ${PYCORE()} ${params.repo}/scripts/coassembly/coasm_tables.py --sheet ${sheet} \\
        --groups ${[params.cfg.coassembly.groups].flatten().collect { g -> "'${g}'" }.join(' ')} \\
        --site-groups ${[params.cfg.binning.groups].flatten().collect { g -> "'${g}'" }.join(' ')} \\
        --mimag-dir mimag --coverm-dir coverm --mag-list ${mag_list} --gtdbtk ${gtdbtk} --mag-amr ${mag_amr} \\
        --compare ${compare} --genomes ${genomes} --site-mags ${site_mags_table} --site-mag-amr ${site_mag_amr} \\
        --min-coverage ${params.cfg.coassembly.present_min_coverage} --out-dir .
    """

    stub:
    """
    touch coasm_mags.tsv coasm_run_hits.tsv coasm_compare.tsv coasm_compare_mags.tsv
    """
}
