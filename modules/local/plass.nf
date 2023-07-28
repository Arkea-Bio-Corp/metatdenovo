process PLASS {
    tag "$meta.id"

    container "quay.io/biocontainers/plass:4.687d7--pl5321h6a68c12_5"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fa")       , emit: fasta
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
    plass nuclassemble \\
    $reads_args \\
    ${prefix}_assembly.fa \\
    tmp \\
    --threads $task.cpus \\
    --compressed 0 \\
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plass: \$(echo \$(plass | head -3 | tail -1 | sed 's/^Plass Version: //'))
    END_VERSIONS
    """
}
