process COUNTS_PLOT {
    tag "$meta.id"

    container "rocker/verse"

    input:
    tuple val(meta), path(counts_txt)

    output:
    tuple val(meta), path(""), emit: counts_png
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    Rscript --vanilla parse_all_counts.R -f $counts_txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( echo \$(Rscript --version 2>&1 | sed 's/.* version//')
    END_VERSIONS
    """
}
