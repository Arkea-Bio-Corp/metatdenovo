// CD-HIT-EST protein clustering
// Taylor Falk tfalk@arkeabio.com
// Arkea Bio Corp, May 2023

include { CDHIT_CDHIT } from './modules/nf-core/cdhit/'

read_ch = Channel.fromPath(params.fasta, checkIfExists: true)
read_ch.view()

workflow CDHITEST {
    reads = read_ch.map  { [[id: 'test'], it]}
    CDHIT_CDHIT(reads)
}

workflow  {
    CDHITEST()
}
