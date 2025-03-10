process TRINITY {
    tag "$meta.id"

    conda "bioconda::trinity=2.15.1"
    container "quay.io/biocontainers/trinity:2.15.1--h6ab5fc9_2"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fa.gz")       , emit: transcript_fasta
    tuple val(meta), path('*.log')         , emit: log, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        reads_args = "--single ${reads}"
    } else {
        reads_args = "--left ${reads[0]} --right ${reads[1]}"
    }

    // --seqType argument, fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    seqType_args = reads[0] ==~ /(.*fasta(.gz)?$)|(.*fa(.gz)?$)/ ? "fa" : "fq"

    // Define the memory requirements. Trinity needs this as an option.
    def avail_mem = task.memory.giga

    """
    # Note that Trinity needs the word 'trinity' in the outdir

    Trinity \\
    --seqType ${seqType_args} \\
    --max_memory ${avail_mem}G \\
    ${reads_args} \\
    --output ${prefix}_trinity \\
    --CPU $task.cpus \\
    $args 

    gzip -cf ${prefix}_trinity.Trinity.fasta > ${prefix}.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trinity: \$(echo \$(Trinity --version | head -n 1 2>&1) | sed 's/^Trinity version: Trinity-v//' )
    END_VERSIONS

    # Need to only take the first line of --version since it will warn about not being up-to-date and this messes up the version.yaml.
    """
}
