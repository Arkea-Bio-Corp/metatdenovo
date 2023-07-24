process PLOT_CONTIGS {
    tag "$meta.id"

    container "public.ecr.aws/r8r1f0u4/contig_counter:latest"

    input:
    tuple val(meta), path(assembly)

    output:
    path("*.png")                , emit: contig_histo
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    contig_counts.R $counts_txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( echo \$(Rscript --version 2>&1 | sed 's/.* version//')
    END_VERSIONS
    """
}
