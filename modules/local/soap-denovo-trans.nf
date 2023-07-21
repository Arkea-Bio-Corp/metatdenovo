process SOAP_DENOVO_TRANS {
    tag "$meta.id"

    container "quay.io/biocontainers/soapdenovo-trans:1.03--he4a0461_2"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.transcripts.fa.gz')  , optional:true, emit: transcripts
    tuple val(meta), path('*.assembly.gfa.gz')    , optional:true, emit: gfa
    tuple val(meta), path('*.log')                , emit: log
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()
    def illumina_reads = reads ? ( meta.single_end ? "-s $reads" : "-1 ${reads[0]} -2 ${reads[1]}" ) : ""

    """
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(rnaspades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS
    """
}