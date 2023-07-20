process RNASPADES {
    tag "$meta.id"

    conda "bioconda::spades=3.15.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spades:3.15.5--h95f258a_1' :
        'quay.io/biocontainers/spades:3.15.5--h95f258a_1' }"

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
    rnaspades.py \\
        $args \\
        --threads $task.cpus \\
        --memory $maxmem \\
        $illumina_reads \\
        -o ./
    mv spades.log ${prefix}.spades.log

    if [ -f transcripts.fasta ]; then
        mv transcripts.fasta ${prefix}.transcripts.fa
        gzip -n ${prefix}.transcripts.fa
    fi
    if [ -f assembly_graph_with_scaffolds.gfa ]; then
        mv assembly_graph_with_scaffolds.gfa ${prefix}.assembly.gfa
        gzip -n ${prefix}.assembly.gfa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(rnaspades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS
    """
}