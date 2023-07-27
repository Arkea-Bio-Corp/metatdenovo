process TRANSRATE {
    tag "$meta.id"

    container "quay.io/biocontainers/transrate:1.0.3--h031d066_5"

    input:
    tuple val(meta), path(assembly)

    output:
    path("*assemblies.csv")      , emit: assembly_qc
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    // tolerate both gzipped and normal inputs
    def gzipped = assembly.toString().endsWith(".gz") ? true : false
    """
    if $gzipped; then
      gzip -d -c $assembly > assembly_unzipped.fa
    else
      mv $assembly assembly_unzipped.fa
    fi

    transrate \\
      --threads $task.cpus \\
      --assembly assembly_unzipped.fa \\
      --output ${prefix}_transrate \\
      $args

    mv ${prefix}_transrate/assemblies.csv ${prefix}_assemblies.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transrate: 1.0.3
    END_VERSIONS
    """
    // transrate isn't getting updated and I can't figure out how to print the version
}
