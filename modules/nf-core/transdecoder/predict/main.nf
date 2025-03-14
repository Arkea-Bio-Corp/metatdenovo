process TRANSDECODER_PREDICT {
    tag "$meta.id"
    cache false

    conda "bioconda::transdecoder=5.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl5262hdfd78af_4':
        'quay.io/biocontainers/transdecoder:5.7.0--pl5321hdfd78af_0'}"

    input:
    tuple val(meta), path(fasta)
    path(folder)

    output:
    tuple val(meta), path("*.transdecoder.pep")  , emit: pep
    tuple val(meta), path("*.transdecoder.gff3") , emit: gff3
    tuple val(meta), path("*.transdecoder.cds")  , emit: cds
    tuple val(meta), path("*.transdecoder.bed")  , emit: bed
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    TransDecoder.Predict \\
        $args \\
        --no_refine_starts \\
        -O ${prefix} \\
        -t \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(echo \$(TransDecoder.Predict --version) | sed -e "s/TransDecoder.Predict //g")
    END_VERSIONS
    """
}
