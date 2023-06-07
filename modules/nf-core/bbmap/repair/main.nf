process BBMAP_REPAIR {
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
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def all_reads      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def repaired_reads  = meta.single_end ? "out=${prefix}_rep.fastq.gz" : "out1=${prefix}_1_rep.fastq.gz out2=${prefix}_2_rep.fastq.gz"
    """
    repair.sh \\
        -Xmx${task.memory.toGiga()}g \\
        $all_reads \\
        $repaired_reads \\
        $args \\
        &> ${prefix}.repair.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
