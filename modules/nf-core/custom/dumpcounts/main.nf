process CUSTOM_DUMPCOUNTS {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
    path collated_reads
    tuple val(meta), val("null")

    output:
    tuple val(meta), path("pipeline_counts.csv"), emit: readcounts

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat collated_counts.csv > pipeline_counts.csv
    """
}
