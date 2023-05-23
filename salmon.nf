// Salmon indexing and quantification
// Taylor Falk tfalk@arkeabio.com
// Arkea Bio Corp, May 2023

include { SALMON_INDEX } from './modules/nf-core/salmon/index'
include { SALMON_QUANT } from './modules/nf-core/salmon/quant'

txome = Channel.fromPath(params.transcriptome, checkIfExists: true)
read_ch = Channel.fromFilePairs(params.pairedreads, checkIfExists: true)
txome.view()
read_ch.view()

workflow SALMONY {
    reads = read_ch.map  { [[id: 'test'], it[1]]}
    salmon_ind = SALMON_INDEX(txome).index
    SALMON_QUANT(reads, salmon_ind)
}

workflow  {
    SALMONY()
}
