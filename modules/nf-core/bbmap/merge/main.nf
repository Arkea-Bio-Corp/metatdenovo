process BBMAP_MERGE {
    tag "$meta.id"

    conda "bioconda::bbmap=39.01"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.merged.fq.gz')   , emit: merged
    tuple val(meta), path('*.unmerged.fq.gz') , emit: unmerged
    tuple val(meta), path('*.log')            , emit: log
    tuple val(meta), path('*.ihist.txt')      , emit: ihist
    path "versions.yml"                       , emit: versions
    path "counts.csv"                         , emit: readcounts

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def all_reads = meta.single_end ? 
        "in=${reads[0]}" : 
        "in1=${reads[0]} in2=${reads[1]}"
    """
    echo $meta > meta.log
    bbmerge.sh \\
        -Xmx${task.memory.toGiga()}g \\
        $all_reads \\
        out=${prefix}.merged.fq.gz \\
        outu=${prefix}.unmerged.fq.gz \\
        threads=$task.cpus \\
        ihist=${prefix}.ihist.txt \\
        $args \\
        &> ${prefix}.merge.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    cat <<-END_COUNTS > counts.csv
    "${task.process}", \$(zcat ${prefix}.merged.fq.gz | grep -c "@" )
    END_COUNTS
    """
}
