// Dedupe - deduplication of sequencing reads
// Laura Holland lholland@arkeabio.com
// Arkea Bio Corp, May 2023

include { BBMAP_DEDUPE } from './modules/nf-core/bbmap/dedupe/'

read_ch = Channel.fromFilePairs(params.pairedreads, checkIfExists: true)

read_ch.view()

workflow DEDUPE {
    reads = read_ch.map  { [[id: 'test'], it[1]]}
    BBMAP_DEDUPE(reads)
}

workflow  {
    DEDUPE()
}
