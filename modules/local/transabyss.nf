process TRANS_ABYSS {
    tag "$meta.id"

    container "quay.io/biocontainers/transabyss:2.0.1--pyh864c0ab_7"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*-final.fa.gz')        , emit: transcripts
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    transabyss \\
    --se SAMPLE1_PE.fastq.gz \\
    --outdir ./ \\
    --threads $task.cpus \\
    --name $prefix \\
    $args

    gzip -n ${prefix}-final.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trans-abyss: \$(echo \$(transabyss --version))
    END_VERSIONS
    """
}
