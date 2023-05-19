// Trinity assemble de novo transcriptome
// Taylor Falk tfalk@arkeabio.com
// Arkea Bio Corp, May 2023

include { TRINITY } from './modules/nf-core/trinity'

read_ch = Channel.fromFilePairs(params.pairedreads, checkIfExists: true)

read_ch.view()

workflow TRINITY_TRIN {
    reads = read_ch.map  { [[id: 'test'], it[1]]}
    TRINITY(reads)
}

workflow  {
    TRINITY_TRIN()
}
