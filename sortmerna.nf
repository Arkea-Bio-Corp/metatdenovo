// SortMeRNA rRNA filtering
// Taylor Falk tfalk@arkeabio.com
// Arkea Bio Corp, May 2023

include { SORTMERNA } from './modules/nf-core/sortmerna/'

read_ch  = Channel.fromFilePairs(params.pairedreads, checkIfExists: true)
ref_ch   = Channel.fromPath(params.reference, checkIfExists: true)
index_ch = Channel.fromPath(params.index, checkIfExists: true)

read_ch.view()
ref_ch.view()

workflow RRNA_REMOVE {
    reads  = read_ch.map  { [[id: 'test'], it[1]]}
    fastas = ref_ch
    index  = index_ch
    SORTMERNA(reads, fastas, index_ch)
}

workflow  {
    RRNA_REMOVE()
}
