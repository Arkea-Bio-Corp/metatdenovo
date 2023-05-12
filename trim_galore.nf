// Trim Galore! Illumina adapter trimming
// Taylor Falk tfalk@arkeabio.com
// Arkea Bio Corp, May 2023

include { TRIMGALORE } from './modules/nf-core/trimgalore/'

read_ch = Channel.fromFilePairs(params.pairedreads, checkIfExists: true)

read_ch.view()

workflow TRIMMYTRIM {
    reads = read_ch.map  { [[id: 'test'], it[1]]}
    TRIMGALORE(reads)
}

workflow  {
    TRIMMYTRIM()
}
