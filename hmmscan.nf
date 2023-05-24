// HMMER functional annotation
// Taylor Falk tfalk@arkeabio.com
// Arkea Bio Corp, May 2023

include { HMMER_HMMSCAN } from './modules/local/hmmscan'

peptide_ch = Channel.fromPath(params.peptide, checkIfExists: true)
hmmerdir   = Channel.fromPath(params.hmmerdir, checkIfExists: true)
hmmerfile  = Channel.value(params.hmmer_file_name)
peptide_ch.view()
hmmerdir.view()
hmmerfile.view()

workflow HMMERTIME {
    peptidemap = peptide_ch.map { [[id: 'test'], it]}
    HMMER_HMMSCAN(peptidemap, hmmerdir, hmmerfile)
}

workflow  {
    HMMERTIME()
}
