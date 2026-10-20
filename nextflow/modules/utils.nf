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
