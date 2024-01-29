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
include { COUNTS_PLOT                      } from '../modules/local/counts_plot'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK     } from '../subworkflows/local/input_check'

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

include { CAT_FASTQ 	          	       } from '../modules/nf-core/cat/fastq/'
include { FASTQC as PRE_TRIM_FQC           } from '../modules/nf-core/fastqc/'
include { FASTQC as POST_MERGE_FQC         } from '../modules/nf-core/fastqc/'
include { MULTIQC                          } from '../modules/nf-core/multiqc/'
include { CUSTOM_DUMPSOFTWAREVERSIONS      } from '../modules/nf-core/custom/dumpsoftwareversions/'
include { CUSTOM_DUMPCOUNTS                } from '../modules/nf-core/custom/dumpcounts/'
include { CUSTOM_DUMPLOGS as TRM_LOGS      } from '../modules/nf-core/custom/dumplogs/'
include { CUSTOM_DUMPLOGS as SMR_LOGS      } from '../modules/nf-core/custom/dumplogs/'
include { CUSTOM_DUMPLOGS as KR2_LOGS      } from '../modules/nf-core/custom/dumplogs/'
include { CUSTOM_DUMPLOGS as BT2_LOGS      } from '../modules/nf-core/custom/dumplogs/'
include { BOWTIE2_ALIGN                    } from '../modules/nf-core/bowtie2/align/'
include { BBMAP_DEDUPE                     } from '../modules/nf-core/bbmap/dedupe/'
include { BBMAP_REPAIR                     } from '../modules/nf-core/bbmap/repair/'
include { BBMAP_REFORMAT                   } from '../modules/nf-core/bbmap/reformat/'
include { BBMAP_MERGE                      } from '../modules/nf-core/bbmap/merge/'
include { CDHIT_CDHIT                      } from '../modules/nf-core/cdhit/'
include { KRAKEN2_KRAKEN2 as KRKN_ARCH     } from '../modules/nf-core/kraken2/'
include { KRAKEN2_KRAKEN2 as KRKN_NO_ARCH  } from '../modules/nf-core/kraken2/'
include { SALMON_INDEX                     } from '../modules/nf-core/salmon/index/'
include { SALMON_QUANT                     } from '../modules/nf-core/salmon/quant/'
include { SORTMERNA                        } from '../modules/nf-core/sortmerna/'
include { TRANSDECODER_LONGORF             } from '../modules/nf-core/transdecoder/longorf/'
include { TRANSDECODER_PREDICT             } from '../modules/nf-core/transdecoder/predict/'
include { TRIMGALORE                       } from '../modules/nf-core/trimgalore/'
include { TRIMMOMATIC                      } from '../modules/nf-core/trimmomatic'
include { TRINITY                          } from '../modules/nf-core/trinity/'
include { SOAP_DENOVO_TRANS                } from '../modules/local/soap-denovo-trans'
include { MEGAHIT                          } from '../modules/nf-core/megahit/'
include { TRANS_ABYSS                      } from '../modules/local/transabyss'
include { PLOT_CONTIGS                     } from '../modules/local/plot_contigs'
include { TRANSRATE                        } from '../modules/local/transrate'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow METATDENOVO {

    ch_versions = Channel.empty()

    // 
    // Clustering with CD-HIT-EST to remove redundancies
    // 
    // Channel
    //     .of(params.assembled_contigs)
    //     .map{ [[id: "megahit", single_end: true], it]}
    //     .set { contig_map }
    // println contig_map

    Channel
        .fromPath(params.transdecoder_pep)
        .splitFasta(by: 50000, file: true, compress: false)
        .map { [[id: "megahit", single_end: true], it] }
        .set { split_assemb }

    // TRANSRATE(contig_map)
    // PLOT_CONTIGS(contig_map, "mega🥺hit")

    // 
    // ORF prediction and translation to peptide seqs with Transdecoder
    // 
    // tdecoder_folder = TRANSDECODER_LONGORF(split_assemb).folder
    // ch_versions = ch_versions.mix(TRANSDECODER_LONGORF.out.versions)
    // TRANSDECODER_PREDICT(split_assemb, tdecoder_folder)
    // ch_versions = ch_versions.mix(TRANSDECODER_PREDICT.out.versions)

    // 
    // Functional annotation with eggnog-mapper
    // 
    eggdbchoice = ["diamond", "hmmer", "novel_fams"]
    // eggnog_ch = Channel.fromPath(params.eggnogdir, checkIfExists: true)
    eggnog_ch = Channel.value(file(params.eggnogdir, checkIfExists: true))
    EGGNOG_MAPPER(split_assemb, eggnog_ch, eggdbchoice)
    ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions)

    // 
    // Functional annotation with hmmscan
    // 
    // hmmerdir   = Channel.fromPath(params.hmmdir, checkIfExists: true)
    // hmmerdir   = Channel.value(file(params.hmmdir, checkIfExists: true))
    // hmmerfile  = Channel.value(params.hmmerfile)
    // HMMER_HMMSCAN(split_assemb, hmmerdir, hmmerfile)
    // ch_versions = ch_versions.mix(HMMER_HMMSCAN.out.versions)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END :^o
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
 
