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

include { CAT_FASTQ 	          	        } from '../modules/nf-core/cat/fastq/'
include { FASTQC as PRE_TRIM_FQC            } from '../modules/nf-core/fastqc/'
include { FASTQC as POST_TRIM_FQC           } from '../modules/nf-core/fastqc/'
include { MULTIQC                           } from '../modules/nf-core/multiqc/'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/'
include { BOWTIE2_ALIGN                     } from '../modules/nf-core/bowtie2/align/'
include { BBMAP_DEDUPE                      } from '../modules/nf-core/bbmap/dedupe/'
include { BBMAP_REPAIR                      } from '../modules/nf-core/bbmap/repair/'
include { BBMAP_REFORMAT                    } from '../modules/nf-core/bbmap/reformat/'
include { CDHIT_CDHIT                       } from '../modules/nf-core/cdhit/'
include { KRAKEN2_KRAKEN2 as KRKN_ARCH      } from '../modules/nf-core/kraken2/'
include { KRAKEN2_KRAKEN2 as KRKN_NO_ARCH   } from '../modules/nf-core/kraken2/'
include { SALMON_INDEX                      } from '../modules/nf-core/salmon/index/'
include { SALMON_QUANT                      } from '../modules/nf-core/salmon/quant/'
include { SORTMERNA                         } from '../modules/nf-core/sortmerna/'
include { TRANSDECODER_LONGORF              } from '../modules/nf-core/transdecoder/longorf/'
include { TRANSDECODER_PREDICT              } from '../modules/nf-core/transdecoder/predict/'
include { TRIMGALORE                        } from '../modules/nf-core/trimgalore/'
include { TRIMMOMATIC                       } from '../modules/nf-core/trimmomatic'
include { TRINITY                           } from '../modules/nf-core/trinity/'

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
    PRE_TRIM_FQC (ch_fastq[0])
    ch_versions = ch_versions.mix(PRE_TRIM_FQC.out.versions)

    // Step 2* Multi QC of raw reads
    // * see below

    // Step 3 Trim Galore!
    //
    // TRIMGALORE(ch_fastq[0])
    // ch_versions = ch_versions.mix(TRIMGALORE.out.versions)

    // Step 3 trimmomatic
    //
    TRIMMOMATIC(ch_fastq[0])
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)


    // 
    // Step 3a FastQC & MultiQC again to compared trimmed reads
    //
    POST_TRIM_FQC(TRIMMOMATIC.out.trimmed_reads)
    ch_versions = ch_versions.mix(POST_TRIM_FQC.out.versions)

    // Step 4
    // Remove host sequences, bowtie2 align to Bos taurus
    // 
    index_ch = Channel.fromPath(params.indexdir)
    BOWTIE2_ALIGN(TRIMMOMATIC.out.trimmed_reads, index_ch, true, false)
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

    // SPLIT ~~~~
    BOWTIE2_ALIGN.out.fastq
        .map { [it.get(0), it.get(1)[0], it.get(1)[1]] }
        // TODO: make pe dynamic with meta.single_end
        .splitFastq(by: params.split_size, 
                    pe: true, file: true, 
                    compress: true, decompress: true)
        .map{ [ it.get(0), [it.get(1), it.get(2)] ]}
        .set { split_reads }

    // Step 5 
    // rRNA remove (bowtie)
    // 
    silva_ch = Channel.value(file(params.silva_reference, checkIfExists: true))
    rna_idx  = Channel.value(file(params.rna_idx, checkIfExists: true))
    SORTMERNA(split_reads, silva_ch, rna_idx)
    ch_versions = ch_versions.mix(SORTMERNA.out.versions)


    // Step 6
    // Filter by taxa with Kraken2
    // 
    k2db_ch = Channel.value(file(params.no_archaea_db, checkIfExists: true))
    KRKN_NO_ARCH(SORTMERNA.out.reads, k2db_ch, true, true)
    ch_versions = ch_versions.mix(KRKN_NO_ARCH.out.versions)

    // RECOMBINE ~~~~
    KRKN_NO_ARCH.out.unclassified_reads_fastq
        .groupTuple()
        .map { [it[0], it[1].flatten()] }
        .set { collectedFastqs }
    CAT_FASTQ(collectedFastqs)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    // Step 7
    // Deduplication with Dedupe
    // 
    BBMAP_REPAIR(CAT_FASTQ.out.reads)
    BBMAP_DEDUPE(BBMAP_REPAIR.out.reads)
    BBMAP_REFORMAT(BBMAP_DEDUPE.out.reads)
    ch_versions = ch_versions.mix(BBMAP_DEDUPE.out.versions)

    // Step 8
    // Merge reads, normalize, and assemble with Trinity
    // 
    TRINITY(BBMAP_REFORMAT.out.reads)
    ch_versions = ch_versions.mix(TRINITY.out.versions)

    // Step 9a 
    // concatenate multiple assemblies
    // CAT_CAT() module? should we make a subworkflow for this?

    // Step 9
    // Clustering with CD-HIT-EST to remove redundancies
    // 
    CDHIT_CDHIT(TRINITY.out.transcript_fasta)
    ch_versions = ch_versions.mix(CDHIT_CDHIT.out.versions)

    // Step 10
    // ORF prediction and translation to peptide seqs with Transdecoder
    // 
    tdecoder_folder = TRANSDECODER_LONGORF(CDHIT_CDHIT.out.fasta).folder
    ch_versions = ch_versions.mix(TRANSDECODER_LONGORF.out.versions)
    TRANSDECODER_PREDICT(CDHIT_CDHIT.out.fasta, tdecoder_folder)
    ch_versions = ch_versions.mix(TRANSDECODER_PREDICT.out.versions)


    // Step 11 
    // Quantification w/ salmon
    // 
    salmon_ind = SALMON_INDEX(CAT_FASTQ.out.reads).index
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
    SALMON_QUANT(CAT_FASTQ.out.reads, salmon_ind)   
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)

    // Step 12 
    // Functional annotation with eggnog-mapper
    // 
    eggdbchoice = ["diamond", "mmseqs", "hmmer", "novel_fams"]
    eggnog_ch = Channel.fromPath(params.eggnogdir, checkIfExists: true)
    EGGNOG_MAPPER(TRANSDECODER_PREDICT.out.pep, eggnog_ch, eggdbchoice)
    ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions)

    // Step 13
    // Functional annotation with hmmscan
    // 
    hmmerdir   = Channel.fromPath(params.hmmdir, checkIfExists: true)
    hmmerfile  = Channel.value(params.hmmerfile)
    HMMER_HMMSCAN(TRANSDECODER_PREDICT.out.pep, hmmerdir, hmmerfile)
    ch_versions = ch_versions.mix(HMMER_HMMSCAN.out.versions)

    // Step 14
    // Kraken2 taxonomical annotation of contigs
    // 
    // adjust metamap to switch to single end
    trans_cds = TRANSDECODER_PREDICT.out.cds.map{ [[id: it[0].id, single_end: true], it[1]]} 
    k2_arch_db = Channel.fromPath(params.archaea_db, checkIfExists: true)
    KRKN_ARCH(trans_cds, k2_arch_db, true, true)
    ch_versions = ch_versions.mix(KRKN_ARCH.out.versions)


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
    ch_multiqc_files = ch_multiqc_files.mix(PRE_TRIM_FQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(POST_TRIM_FQC.out.zip.collect{it[1]}.ifEmpty([]))


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
