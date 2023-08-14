process CAT_QUANTS {
    label 'process_low'

    container "quay.io/biocontainers/grep:3.4--hf43ccf4_4"

    input:
    path(quants, stageAs: "?/*")

    output:
    path("combined_quants.sf")      , emit: combined

    script:
    def qlist = quants instanceof List ? 
        quants.collect{ it.toString() } : 
        [quants.toString()]
    """
    # head -1 $quants > combined_quants.sf
    grep -vh "^Name" ${qlist.join(' ')} >> combined_quants.sf
    """
}
