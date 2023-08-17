process SALMON_MERGE {
    container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"

    input:
    path(quants, stageAs: "?/*")

    output:
    path  "*.sf"                       , emit: quants
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''

    """
    quant_dirs="\$(ls -d -- */)"

    salmon quantmerge \\
        --quants \$(echo $quants | grep -o "[0-9]\\+/" | awk '\$1=\$1' ORS=' ') \\
        --names $quants \\
        --column numreads \\
        --output salmon_qaunt_out.sf \\
        $args 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
