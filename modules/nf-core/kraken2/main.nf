process KRAKEN2_KRAKEN2 {
    tag "$meta.id"

    conda "bioconda::kraken2=2.1.2 conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' :
        'quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' }"

    input:
    tuple val(meta), path(reads)
    path db
    val save_output_fastqs
    val save_reads_assignment

    output:
    tuple val(meta), path('*.classified{.,_}*')  , emit: classified_reads_fastq
    tuple val(meta), path('*.unclassified{.,_}*'), emit: unclassified_reads_fastq
    tuple val(meta), path('*classifiedreads.txt'), optional:true, emit: classified_reads_assignment
    tuple val(meta), path('*report.txt')         , emit: report
    tuple val(meta), val("null")                 , emit: meta // passing meta tag only to others
    path('*report.txt')                          , emit: collect_log // only report to collect
    path "versions.yml"                          , emit: versions
    path "counts.csv"                            , emit: readcounts

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired       = meta.single_end ? "" : "--paired"
    def classified   = meta.single_end ? "${prefix}.classified.fastq"   : "${prefix}.classified#.fastq"
    def unclassified = meta.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    def classified_option = save_output_fastqs ? "--classified-out ${classified}" : ""
    def unclassified_option = save_output_fastqs ? "--unclassified-out ${unclassified}" : ""
    def readclassification_option = save_reads_assignment ? "--output ${prefix}.kraken2.classifiedreads.txt" : "--output /dev/null"
    def compress_reads_command = save_output_fastqs ? "pigz -p $task.cpus *.fastq" : ""
    def gzswitch = reads.toString().endsWith(".gz") ? "--gzip-compressed" : ""

    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}.kraken2.report.txt \\
        $gzswitch \\
        $unclassified_option \\
        $classified_option \\
        $readclassification_option \\
        $paired \\
        $args \\
        $reads

    $compress_reads_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    cat <<-END_COUNTS > counts.csv
    "${task.process}", \$(zcat *unclassified* | grep -c "@" )
    END_COUNTS
    """
}
