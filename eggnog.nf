// eggnog functional annotation
// Taylor Falk tfalk@arkeabio.com
// Arkea Bio Corp, May 2023

include { EGGNOG_MAPPER } from './modules/local/eggnog/mapper.nf'

read_ch     = Channel.fromPath(params.reads, checkIfExists: true)
database_ch = Channel.fromPath(params.databases, checkIfExists: true)

workflow MAPPY {
    dbchoice = ["diamond", "mmseqs", "hmmer", "novel_fams"]
    read_ch = read_ch.map  { [[id: 'test'], it]}
    EGGNOG_MAPPER(read_ch, database_ch, dbchoice)
}


workflow  {
    MAPPY()
}
