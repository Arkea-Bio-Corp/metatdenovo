//
// Build and align assembly to input reads to determine coverage percentage.
//

include { BOWTIE2_ALIGN as BT2_TRNS_ALGN } from '../../modules/nf-core/bowtie2/align/'
include { BOWTIE2_BUILD as BT2_TRNS_BLD  } from '../../modules/nf-core/bowtie2/build/'

workflow BT2_ALIGN_STATS {
    take:
        assembly
        input_reads
    main:
        ch_versions = Channel.empty()

        BT2_TRNS_BLD(assembly)
        BT2_TRNS_BLD.out.index
            .map { it[1] }
            .set { index_path }

        BT2_TRNS_ALGN(input_reads, index_path, true, false)

        // ch_versions = ch_versions.mix(EGGNOG_TABLE.out.versions)
        // ch_versions = ch_versions.mix(SUM_EGGNOG.out.versions)

    emit:
        aligned         = BT2_TRNS_ALGN.out.collect_log
        // versions          = ch_versions

}
