process CUSTOM_DUMPCOUNTS {

    input:
    path collated_reads

    output:
    path "counts.yml", emit: readcounts

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat collated_counts.yml > counts.yml
    """
}
