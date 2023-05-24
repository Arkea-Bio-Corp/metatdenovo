// eggnog functional annotation
// Taylor Falk tfalk@arkeabio.com
// Arkea Bio Corp, May 2023

include { EGGNOG_MAPPER } from './modules/local/eggnog/mapper.nf'

read_ch     = Channel.fromPath(params.reads, checkIfExists: true)
database_ch = Channel.fromPath(params.databases, checkIfExists: true)
read_ch.view()
database_ch.view()

workflow MAPPY {
    reads  = read_ch.map  { [[id: 'test'], it]}
    EGGNOG_MAPPER(reads, database_ch)
}

workflow  {
    MAPPY()
}
