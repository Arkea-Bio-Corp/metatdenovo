process SALMON_QUANT {
    tag "$meta.id"
    label "process_medium"

    conda "bioconda::salmon=1.10.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.1--h7e5ed60_0' :
        'quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0' }"

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("${prefix}") , emit: results
    tuple val(meta), path("*info.json"), emit: json_info, optional: true
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    def input_reads = meta.single_end ? "-r $reads" : "-1 ${reads[0]} -2 ${reads[1]}"

    """
    salmon quant \\
        --threads $task.cpus \\
        --libType A \\
        --index $index \\
        $input_reads \\
        $args \\
        -o $prefix

    if [ -f $prefix/aux_info/meta_info.json ]; then
        cp $prefix/aux_info/meta_info.json "${prefix}_meta_info.json"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
