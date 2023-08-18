process SALMON_MERGE {
    container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"

    input:
    tuple val(metas), path(quants, stageAs: "?/quant.sf")

    output:
    path  "*.sf"                       , emit: quants
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''

    """
    salmon quantmerge \\
        --quants \$(echo $quants | grep -o "[0-9]\\+/" | awk '\$1=\$1' ORS=' ') \\
        --names ${metas.collect{ it.id }.join(" ")} \\
        --column numreads \\
        --output salmon_quant_out.sf \\
        $args 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
