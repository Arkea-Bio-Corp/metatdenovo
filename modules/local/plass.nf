process PLASS {
    tag "$meta.id"

    container "soedinglab/plass"
    docker.runOptions = "--entrypoint /bin/bash"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fa.gz")       , emit: fasta
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        reads_args = "${reads}"
    } else {
        reads_args = "${reads[0]} ${reads[1]}"
    }

    """
    plass assemble \\
    $reads_args \\
    ${prefix}_assembly.fa.gz \\
    tmp \\
    --split-memory-limit ${task.memory.giga}G
    --threads $task.cpus \\
    --compress 1 \\
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plass: \$(echo \$(plass | head -3 | tail -1 | sed 's/^Plass Version: //'))
    END_VERSIONS
    """
}
