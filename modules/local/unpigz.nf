process UNPIGZ {
    tag "$file"

    // conda (params.enable_conda ? "pigz=2.3.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4':
        'quay.io/biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("$gunzip") , emit: unzipped
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    gunzip = file.toString() - '.gz'

    """
    unpigz \\
        -c \\
        -p $task.cpus \\
        ${file} > $gunzip
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: 2.3.4  \$( pigz --version)
    END_VERSIONS
    """
}
