//
// Build and align assembly to input reads to determine coverage percentage.
//

include { BOWTIE2_ALIGN as BT2_TRNS_ALGN } from '../../modules/nf-core/bowtie2/align/'
include { BOWTIE2_BUILD as BT2_TRNS_BLD  } from '../../modules/nf-core/bowtie2/build/'
include { PLOT_CONTIGS                   } from '../../modules/local/plot_contigs'

workflow ASSEMBLE_STATS {
    take:
        assembly
        input_reads
    main:
        ch_versions = Channel.empty()

        // Bowtie2 stats
        BT2_TRNS_BLD(assembly)
        BT2_TRNS_BLD.out.index
            .map { it[1] }
            .set { index_path }
        BT2_TRNS_ALGN(input_reads, index_path, true, false)

        // Plot contig distribution
        PLOT_CONTIGS(assembly)

        ch_versions = ch_versions.mix(BT2_TRNS_ALGN.out.versions)
        ch_versions = ch_versions.mix(BT2_TRNS_BLD.out.versions)

    emit:
        aligned         = BT2_TRNS_ALGN.out.collect_log
        versions        = ch_versions

}
