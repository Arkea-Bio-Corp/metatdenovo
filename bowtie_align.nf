// Bowtie 2 Bos taurus alignment and filtering
// Taylor Falk tfalk@arkeabio.com
// Arkea Bio Corp, May 2023

include { FASTQ_ALIGN_BOWTIE2 } from './subworkflows/nf-core/fastq_align_bowtie2'
include { BOWTIE2_BUILD       } from './modules/nf-core/bowtie2/build'

params.genome =  file('./bos_taurus/ncbi_dataset/data/GCA_002263795.3/GCA_002263795.3_ARS-UCD1.3_genomic.fna')

workflow BT2_BUILD {
    input = [
        [id:'test'],
        params.genome
    ]
    BOWTIE2_BUILD(input)
}

// Implicit workflow
workflow  {
    BT2_BUILD()
}
