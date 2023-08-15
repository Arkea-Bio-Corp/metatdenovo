process COUNTS_PLOT {
    tag "$meta.id"

    container "rocker/verse"

    input:
    tuple val(meta), path(counts_txt)

    output:
    path("*.png")                               , emit: counts_png
    path("*.html")                              , emit: counts_html
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    parse_all_counts.R $counts_txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( echo \$(Rscript --version 2>&1 | sed 's/.* version//')
    END_VERSIONS
    """
}
