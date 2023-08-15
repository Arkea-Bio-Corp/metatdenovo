process BBMAP_DEDUPE {
    tag "$meta.id"

    conda "bioconda::bbmap=39.01"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    tuple val(meta), val("null")       , emit: meta // passing meta tag only to others
    path "versions.yml"                , emit: versions
    path "counts.csv"                  , emit: readcounts


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def all_reads      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def deduped_reads  = meta.single_end ? "out=${prefix}.fastq.gz" : "out=${prefix}_deduped.fastq.gz"
    """
    echo $meta > meta.log
    dedupe.sh \\
        -Xmx${task.memory.toGiga()}g \\
        $all_reads \\
        $deduped_reads \\
        threads=$task.cpus \\
        $args \\
        &> ${prefix}.dedupe.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    cat <<-END_COUNTS > counts.csv
    "${task.process}",\$(zcat ${prefix}.fastq.gz | grep -c "@" )
    END_COUNTS
    """
}
