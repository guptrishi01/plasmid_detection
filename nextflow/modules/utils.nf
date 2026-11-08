// Helpers for summary processes. The table builders (src/plasmid_detection/tables.py) read
// the per-run layout of the sbatch version (<dir>/<run>/...). Summary processes rebuild that
// layout with symlinks to the per-run task outputs: this assumes a shared filesystem (true
// on the cluster; every work dir is on /scratch).

// runs.txt: one run per line, the sample list the table builders expect
def runsFile(runs) {
    "printf '%s\\n' ${runs.join(' ')} > runs.txt"
}

// pairs: [[run, path], ...]; dest: template with {run}, e.g. 'asm/{run}/quast'.
// Emits shell lines: mkdir -p <parent>; ln -s <absolute path> <dest>
def linkLayout(pairs, String dest) {
    pairs.collect { run, p ->
        def target = dest.replace('{run}', run.toString())
        def parent = target.contains('/') ? target.substring(0, target.lastIndexOf('/')) : '.'
        "mkdir -p '${parent}' && ln -sfn '${p.toAbsolutePath()}' '${target}'"
    }.join('\n')
}

// dRep at params.cfg.mags.drep_ani over the genomes matching `glob` into `out`, CheckM2
// values from genome_info.csv; one genome skips dRep and gets the same layout by hand
def dereplicate(String out, String glob, cpus) {
    """
    n=\$(ls ${glob} | wc -l)
    if (( n == 1 )); then
        mkdir -p ${out}/dereplicated_genomes ${out}/data_tables
        cp -L ${glob} ${out}/dereplicated_genomes/
        g=\$(basename ${glob})
        printf 'genome,secondary_cluster\\n%s,1_1\\n' \$g > ${out}/data_tables/Cdb.csv
        printf 'genome,cluster\\n%s,1_1\\n' \$g > ${out}/data_tables/Wdb.csv
    else
        dRep dereplicate ${out} -g ${glob} -p ${cpus} -sa ${params.cfg.mags.drep_ani} \\
            --genomeInfo genome_info.csv \\
            -comp ${params.cfg.binning.mimag_medium_completeness} -con ${params.cfg.binning.mimag_medium_contamination}
    fi
    """
}
