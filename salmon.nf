// Salmon indexing and quantification
// Taylor Falk tfalk@arkeabio.com
// Arkea Bio Corp, May 2023

include { SALMON_INDEX } from './modules/nf-core/salmon/index'

read_ch = Channel.fromPath(params.transcriptome, checkIfExists: true)
read_ch.view()

workflow SALMONY {
    SALMON_INDEX(read_ch)
}

workflow  {
    SALMONY()
}
