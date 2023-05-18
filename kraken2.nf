// Kraken2 taxonomic sequence identifier
// Taylor Falk tfalk@arkeabio.com
// Arkea Bio Corp, May 2023

include { KRAKEN2_KRAKEN2 } from './modules/nf-core/kraken2/main.nf'

read_ch = Channel.fromFilePairs(params.pairedreads, checkIfExists: true)
db_ch   = Channel.fromPath(params.database, checkIfExists: true)

read_ch.view()

workflow KRAKEN_ID {
    reads = read_ch.map  { [[id: 'test'], it[1]]}
    KRAKEN2_KRAKEN2(reads, db_ch, true, true)
}

workflow  {
    KRAKEN_ID()
}
