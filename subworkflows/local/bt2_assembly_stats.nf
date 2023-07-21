//
// Build and align assembly to input reads to determine coverage percentage.
//

include { BOWTIE2_ALIGN as BT2_TRNS_ALGN } from '../modules/nf-core/bowtie2/align/'
include { BOWTIE2_BUILD AS BT2_TRNS_BLD  } from '../modules/nf-core/bowtie2/build/'

workflow BT2_ALIGN_STATS {
    take:
        assembly
        input_reads
    main:
        ch_versions = Channel.empty()

        BT2_TRNS_BLD

        ch_versions = ch_versions.mix(EGGNOG_TABLE.out.versions)
        ch_versions = ch_versions.mix(SUM_EGGNOG.out.versions)

    emit:
        hits              = EGGNOG_MAPPER.out.hits
        formateggnogtable = EGGNOG_TABLE.out.eggtab
        versions          = ch_versions
        sumtable          = SUM_EGGNOG.out.eggnog_summary

}
