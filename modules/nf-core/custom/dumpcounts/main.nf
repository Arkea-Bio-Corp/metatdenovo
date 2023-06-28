process CUSTOM_DUMPCOUNTS {

    input:
    path collated_reads
    tuple val(meta), val("null")

    output:
    path "counts.yml", emit: readcounts

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat collated_counts.yml > counts.yml
    """
}
