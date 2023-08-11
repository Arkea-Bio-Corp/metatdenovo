process SOAP_DENOVO_TRANS {
    tag "$meta.id"

    container "quay.io/biocontainers/soapdenovo-trans:1.03--he4a0461_2"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.contig')             , emit: contigs
    tuple val(meta), path('*.scafStatistics')     , emit: stats
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()

    """
    gzip -d -c $reads > ${prefix}_unzip.fq
    cat <<-ENDCONF > config.txt
    #maximal read length
    max_rd_len=500
    [LIB]
    #maximal read length in this lib
    rd_len_cutof=500
    #average insert size
    avg_ins=150
    #if sequence needs to be reversed
    reverse_seq=0
    #in which part(s) the reads are used
    asm_flags=1
    map_len=32
    #fastq file for single reads
    q=${prefix}_unzip.fq
    ENDCONF

    SOAPdenovo-Trans-31mer all \\
    -s config.txt \\
    -o $prefix \\
    -R \\
    -p $task.cpus \\
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SOAPdenovo-Trans: \$(SOAPdenovo-Trans-31mer 2>&1 | head -2 | tail -1 | sed 's/^The version //;s/: released .*//')
    END_VERSIONS
    """
}