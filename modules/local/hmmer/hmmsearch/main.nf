process HMMER_HMMSEARCH {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::hmmer=3.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.3.2--h1b792b2_1' :
        'quay.io/biocontainers/hmmer:3.3.2--h1b792b2_1' }"

    input:
    tuple val(meta), path(seqdb)
    path hmmdir
    val  hmmpattern

    output:
    tuple val(meta), path('*.tblout'), emit: tblout
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    hmmsearch \\
        --tblout $hmmdir/$hmmpattern".tblout" \\
        $hmmdir/$hmmpattern \\
        $seqdb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}
