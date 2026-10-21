// Stage 5: geNomad (main plasmid caller), CheckV + MOB-typer, AMRFinderPlus + RGI,
// CoverM contig coverage and AMR tagging (plasmid / virus / chromosome)

process GENOMAD {
    tag "${meta.id}"
    conda params.envs.genomad
    publishDir path: { "${params.outdir}/05_contigs/${meta.id}" }, mode: 'copy', pattern: 'genomad'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path('genomad'), emit: genomad

    script:
    """
    genomad end-to-end ${params.cfg.contigs.genomad_args} --threads ${task.cpus} ${contigs} genomad ${params.db.genomad}
    """

    stub:
    """
    mkdir -p genomad
    """
}

process CHECKV_MOBTYPER {
    tag "${meta.id}"
    conda params.envs.checkv
    publishDir path: { "${params.outdir}/05_contigs/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(genomad)

    output:
    tuple val(meta), path('mob_typer.tsv'), emit: mob, optional: true
    path 'checkv'
    path 'mob_typer.EMPTY', optional: true

    script:
    """
    virus=\$(ls ${genomad}/*_summary/*_virus.fna)
    plasmid=\$(ls ${genomad}/*_summary/*_plasmid.fna)
    if [[ -s \$virus ]]; then
        checkv end_to_end \$virus checkv -d ${params.db.checkv} -t ${task.cpus}
        rm -rf checkv/tmp
    else
        mkdir -p checkv; echo "no virus contigs" > checkv/EMPTY
    fi
    if [[ -s \$plasmid ]]; then
        # mob_typer calls blastn/mash/etc. by name: its env's bin must lead PATH
        PATH=${params.envs.mobsuite}/bin:\$PATH mob_typer --multi --infile \$plasmid --out_file mob_typer.tsv \\
            --database_directory ${params.db.mob_suite} --num_threads ${task.cpus}
    else
        echo "no plasmid contigs" > mob_typer.EMPTY
    fi
    """

    stub:
    """
    mkdir -p checkv; touch mob_typer.tsv
    """
}

process AMR {
    tag "${meta.id}"
    conda params.envs.amrfinderplus
    publishDir path: { "${params.outdir}/05_contigs/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path('amrfinder.tsv'), path("rgi/${meta.id}.rgi.txt"), emit: amr

    script:
    def id = meta.id
    """
    amrfinder --nucleotide ${contigs} ${params.cfg.contigs.amrfinder_args} --database ${params.db.amrfinderplus} \\
        --threads ${task.cpus} --name ${id} --output amrfinder.tsv
    # RGI --local reads ./localDB (CARD), linked in from db/card
    mkdir -p rgi && cd rgi && ln -sfn ${params.db.card.dir}/localDB localDB
    # RGI calls prodigal/diamond by name: its env's bin must lead PATH (else no ORFs ->
    # FileNotFoundError on *.temp.draft)
    PATH=${params.envs.rgi}/bin:\$PATH rgi main --input_sequence ../${contigs} --output_file ${id}.rgi --input_type contig \\
        --alignment_tool DIAMOND --threads ${task.cpus} --local --clean
    rm -f ${id}.rgi.json localDB
    """

    stub:
    """
    mkdir -p rgi; touch amrfinder.tsv rgi/${meta.id}.rgi.txt
    """
}

process COVERAGE_TAG {
    tag "${meta.id}"
    conda params.envs.coverm
    publishDir path: { "${params.outdir}/05_contigs/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(bams), path(bais), path(genomad), path(amrfinder), path(rgi)

    output:
    tuple val(meta), path('amr_tagged.tsv'), path('amr_counts.tsv'), path('coverm_contig.tsv'), emit: tagged
    path 'amr_hits.bed'
    path 'amr_bedcov.txt'

    script:
    def id = meta.id
    def py = "PYTHONPATH=${params.repo}/src ${params.envs.core}/bin/python -m plasmid_detection.tables"
    """
    coverm contig --bam-files ${bams} --methods ${params.cfg.contigs.coverm_methods} \\
        --threads ${task.cpus} --output-file coverm_contig.tsv
    ${py} tag-amr --run ${id} --genomad-dir ${genomad} --amrfinder ${amrfinder} --rgi ${rgi} \\
        --out amr_tagged.tsv --bed amr_hits.bed
    ${params.envs.assembly}/bin/samtools bedcov -c amr_hits.bed ${bams} > amr_bedcov.txt
    ${py} amr-counts --bedcov amr_bedcov.txt --bams ${bams} --out amr_counts.tsv
    """

    stub:
    """
    touch amr_tagged.tsv amr_counts.tsv coverm_contig.tsv amr_hits.bed amr_bedcov.txt
    """
}
