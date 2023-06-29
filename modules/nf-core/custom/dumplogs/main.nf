process CUSTOM_DUMPLOGS {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
    path logs
    val log_type
    tuple val(meta), val("null")

    output:
    path "${log_type}_logs.txt", emit: saved_log

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat collected_logs.txt > ${log_type}_logs.txt
    """
}
