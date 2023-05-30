/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMetatdenovo.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// set an empty multiqc channel
ch_multiqc_files = Channel.empty()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: local
//
include { COLLECT_STATS                    } from '../modules/local/collect_stats'
include { UNPIGZ as UNPIGZ_CONTIGS         } from '../modules/local/unpigz'
include { UNPIGZ as UNPIGZ_GFF             } from '../modules/local/unpigz'
include { HMMER_HMMSCAN                    } from '../modules/local/hmmscan/main'
include { EGGNOG_MAPPER                    } from '../modules/local/eggnog/mapper'
include { EGGNOG_DOWNLOAD                  } from '../modules/local/eggnog/download'

//include { MERGE_TABLES                     } from '../modules/local/merge_summary_tables'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK     } from '../subworkflows/local/input_check'
include { BT2_ALIGN       } from './bowtie_align'
include { CDHITEST        } from './cd_hit_est'
include { DEDUPE          } from './dedupe'
include { MAPPY           } from './eggnog'
include { HMMERTIME       } from './hmmscan'
include { KRAKEN_ID       } from './kraken2'
include { SALMONY         } from './salmon'
include { RRNA_REMOVE     } from './sortmerna'
include { LONGORF_PREDICT } from './transdecoder'
include { TRIMMYTRIM      } from './trim_galore'
include { TRINITY_TRIN    } from './trinity'

//
// SUBWORKFLOW: Consisting of local modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules (mostly)
//
include { CAT_FASTQ 	          	  } from '../modules/nf-core/cat/fastq/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { BBMAP_DEDUPE                } from '../modules/nf-core/bbmap/dedupe/main'
include { BOWTIE2_ALIGN               } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_BUILD               } from '../modules/nf-core/bowtie2/build/main'
include { CDHIT_CDHIT                 } from '../modules/nf-core/cdhit/main'
include { KRAKEN2_KRAKEN2             } from '../modules/nf-core/kraken2/main'
include { SALMON_INDEX                } from '../modules/nf-core/salmon/index/main'
include { SALMON_QUANT                } from '../modules/nf-core/salmon/quant/main'
include { SORTMERNA                   } from '../modules/nf-core/sortmerna/main'
include { TRANSDECODER_LONGORF        } from '../modules/nf-core/transdecoder/longorf/main'
include { TRANSDECODER_PREDICT        } from '../modules/nf-core/transdecoder/predict/main'
include { TRIMGALORE                  } from '../modules/nf-core/trimgalore/main'
include { TRINITY                     } from '../modules/nf-core/trinity/main'

//
// SUBWORKFLOWS: Installed directly from nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow METATDENOVO {

    ch_versions = Channel.empty()

    // STEP 0: Validate input
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            new_id = meta.id - ~/_T\d+/
            [ meta + [id: new_id], fastq ]
    }
    .groupTuple()
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Step 1 FastQC
    //
    FASTQC (
        ch_fastq[0]
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // Step 2 Multi QC of raw reads

    // Step 3 Trim Galore!
    // TRIMMYTRIM()

    // 
    // Step 3a FastQC & MultiQC again to compared trimmed reads
    //

    // Step 4
    // Remove host sequences, bowtie2 align to Bos taurus
    // BT2_ALIGN()

    // Step 5 
    // rRNA remove (sortmerna)
    // RRNA_REMOVE()

    // Step 6
    // Deduplication with Dedupe
    // DEDUPE()

    // Step 7
    // Filter by taxa with Kraken2
    // KRAKEN_ID()

    // Step 8
    // Merge reads, normalize, and assemble with Trinity
    // TRINITY_TRIN()

    // Step 9a 
    // concatenate multiple assemblies
    // CAT_CAT() module? should we make a subworkflow for this?

    // Step 9
    // Clustering with CD-HIT-EST to remove redundancies
    // CDHITEST()

    // Step 10
    // ORF prediction and translation to peptide seqs with Transdecoder
    // LONGORF_PREDICT()

    // Step 11 --> important!! use reads from after step 5 here
    // Quantification (salmon, rsem, or bbmap)
    // SALMONY()

    // Step 12 
    // Functional annotation with eggnog-mapper
    // MAPPY()

    // Step 13
    // Functional annotation with hmmscan
    // HMMERTIME()

    // Step 14
    // Kraken2 taxonomical annotation of contigs
    // KRAKEN_ID()

    // Step 15
    // MODULE: Collect statistics from mapping analysis
    //
    // ch_collect_stats
    //     .map { [ it[0], it[1], it[2], it[3], it[4], it[5], [] ] }
    //     .set { ch_collect_stats }

    // COLLECT_STATS(ch_collect_stats)
    // ch_versions     = ch_versions.mix(COLLECT_STATS.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // Step 2
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMetatdenovo.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMetatdenovo.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS_CDS.out.summary.collect{it[1]}.ifEmpty([]))


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END :^)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
