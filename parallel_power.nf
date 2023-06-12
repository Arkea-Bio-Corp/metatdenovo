// SortMeRNA rRNA filtering
// Taylor Falk tfalk@arkeabio.com
// Arkea Bio Corp, May 2023

include { SORTMERNA                         } from './modules/nf-core/sortmerna/'
include { KRAKEN2_KRAKEN2 as KRKN_NO_ARCH   } from './modules/nf-core/kraken2/'

read_ch  = Channel.fromFilePairs(params.pairedreads, checkIfExists: true)
ref_ch   = Channel.value(file(params.reference,      checkIfExists: true))
index_ch = Channel.value(file(params.index,          checkIfExists: true))
k2db_ch  = Channel.value(file(params.no_archaea_db,  checkIfExists: true))


workflow RRNA_REMOVE {
    read_ch = read_ch.map  { [[id: 'test'], it[1]]}
    read_ch
        .map { [it.get(0), it.get(1)[0], it.get(1)[1]] }
        .splitFastq(by: 1_000, pe: true, file: true, compress: true, decompress: true)
        .map{ [ [id: 'test'], [it.get(1), it.get(2)] ]}
        .set { split_reads }
    SORTMERNA(split_reads, ref_ch, index_ch)
    KRKN_NO_ARCH(SORTMERNA.out.reads, k2db_ch, true, true)
}

workflow  {
    RRNA_REMOVE()
}