// Bowtie 2 Bos taurus alignment and filtering
// Taylor Falk tfalk@arkeabio.com
// Arkea Bio Corp, May 2023

include { BOWTIE2_ALIGN } from './modules/nf-core/bowtie2/align'

// index_ch = Channel.fromPath("$projectDir/bos_taurus/bowtie2_index/")
// read_ch = Channel.fromFilePairs("$projectDir/localdata/*{R1,R2}.fastq.gz")
index_ch = Channel.fromPath(params.indexdir)
read_ch = Channel.fromFilePairs(params.pairedreads, checkIfExists: true)

index_ch.view()
read_ch.view()

workflow BT2_ALIGN {
    reads = read_ch.map  { [[id: 'test'], it[1]]}
    index = index_ch.map { [[id: 'test'], it] }
    BOWTIE2_ALIGN(reads, index, true, false)
}

workflow  {
    BT2_ALIGN()
}
