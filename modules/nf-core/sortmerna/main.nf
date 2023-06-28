process SORTMERNA {
    tag "$meta.id"

    conda "bioconda::sortmerna=4.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sortmerna:4.3.4--h9ee0642_0' :
        'quay.io/biocontainers/sortmerna:4.3.4--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)
    path  fastas
    path index_dir

    output:
    tuple val(meta), path("*non_rRNA.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path  "versions.yml"               , emit: versions
    path  "counts.yml"                 , emit: readcounts

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def all_reads = meta.single_end ? 
        "--reads ${reads}" :
        "--reads ${reads[0]} --reads ${reads[1]}"
    if (meta.single_end) {
        """
        sortmerna \\
            ${'--ref '+fastas.join(' --ref ')} \\
            $all_reads \\
            --threads $task.cpus \\
            --index 0 \\
            --idx-dir ${index_dir} \\
            --workdir . \\
            --aligned rRNA_reads \\
            --fastx \\
            --other non_rRNA_reads \\
            $args

        mv non_rRNA_reads.f*q.gz ${prefix}.non_rRNA.fastq.gz
        mv rRNA_reads.log ${prefix}.sortmerna.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
        END_VERSIONS
        cat <<-END_COUNTS > counts.yml
        "${task.process}_${task.index}":
            \$(zcat *non_rRNA.fastq.gz | grep -c "@" )
        END_COUNTS
        """
    } else {
        """
        sortmerna \\
            ${'--ref '+fastas.join(' --ref ')} \\
            --paired_in \\
            $all_reads \\
            --threads $task.cpus \\
            --index 0 \\
            --idx-dir ${index_dir} \\
            --workdir . \\
            --aligned rRNA_reads \\
            --fastx \\
            --other non_rRNA_reads \\
            --out2 \\
            $args

        mv non_rRNA_reads_fwd.f*q.gz ${prefix}_1.non_rRNA.fastq.gz
        mv non_rRNA_reads_rev.f*q.gz ${prefix}_2.non_rRNA.fastq.gz
        mv rRNA_reads.log ${prefix}.sortmerna.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
        END_VERSIONS
        cat <<-END_COUNTS > counts.yml
        "${task.process}_${task.index}":
            \$(zcat *non_rRNA.fastq.gz | grep -c "@" )
        END_COUNTS
        """
    }
}
