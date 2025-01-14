process EGGNOG_MAPPER {
    tag "$meta.id"
    label 'process_high'

    // conda (params.enable_conda ? "bioconda::eggnog-mapper=2.1.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.9--pyhdfd78af_0':
        'quay.io/biocontainers/eggnog-mapper:2.1.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(db)
    each dbchoice

    output:
    tuple val(meta), path("*.emapper.hits")                , emit: hits
    tuple val(meta), path("*.emapper.seed_orthologs")      , emit: seed_orthologs
    tuple val(meta), path("*.emapper.annotations")         , emit: annotations
    tuple val(meta), path("*.emapper.annotations.xlsx")    , emit: xlsx,      optional: true
    tuple val(meta), path("*.emapper.orthologs")           , emit: orthologs, optional: true
    tuple val(meta), path("*.emapper.genepred.fasta")      , emit: genepred,  optional: true
    tuple val(meta), path("*.emapper.gff")                 , emit: gff,       optional: true
    tuple val(meta), path("*.emapper.no_annotations.fasta"), emit: no_anno,   optional: true
    tuple val(meta), path("*.emapper.pfam")                , emit: pfam,      optional: true
    path "versions.yml"                                    , emit: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    input    = fasta =~ /\.gz$/ ? fasta.name.take(fasta.name.lastIndexOf('.')) : fasta
    gunzip   = fasta =~ /\.gz$/ ? "gunzip -c ${fasta} > ${input}" : ""
    db_label = dbchoice
    dbchoice = dbchoice == "hmmer" ? "hmmer -d Archaea" : dbchoice
    dbchoice = dbchoice == "eggnog" ? "" : "-m ${dbchoice}"
    """
    $gunzip

    emapper.py \\
        $args \\
        --itype proteins \\
        $dbchoice \\
        --cpu $task.cpus \\
        --data_dir $db \\
        --output ${prefix}_${db_label} \\
        -i $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//' | sed 's/\\/ Expected eggNOG DB version: 5.0.2 \\/ Installed eggNOG DB version: unknown \\/ Diamond version found: diamond version 2.0.15 \\/ MMseqs2 version found: 13.45111//g' )
    END_VERSIONS
    """
}
