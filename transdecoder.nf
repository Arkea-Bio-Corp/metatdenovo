// TransDecoder ORF prediction
// Taylor Falk tfalk@arkeabio.com
// Arkea Bio Corp, May 2023

include { TRANSDECODER_LONGORF } from './modules/nf-core/transdecoder/longorf'
include { TRANSDECODER_PREDICT } from './modules/nf-core/transdecoder/predict'

read_ch = Channel.fromPath(params.fasta, checkIfExists: true)
read_ch.view()

workflow LONGORF_PREDICT {
    reads = read_ch.map  { [[id: 'test'], it]}
    tdecorder_folder = TRANSDECODER_LONGORF(reads).folder
    TRANSDECODER_PREDICT(reads, tdecorder_folder)
}

workflow  {
    LONGORF_PREDICT()
}
