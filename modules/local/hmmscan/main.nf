process HMMER_HMMSCAN {
    tag "$meta.id"

    conda "bioconda::hmmer=3.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.3.2--h1b792b2_1' :
        'quay.io/biocontainers/hmmer:3.3.2--h1b792b2_1' }"

    input:
    tuple val(meta), path(pepfile)
    path(hmmdir)
    val(hmmfile)


    output:
    tuple val(meta), path('*.txt')    , emit: output
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    output = "${prefix}.txt"

    """
    hmmscan \\
        $args \\
        --cpu $task.cpus \\
        -o $output \\
        $hmmdir/$hmmfile \\
        $pepfile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmscan -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}
