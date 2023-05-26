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

// If the user supplied hmm files, we will run hmmsearch and then rank the results.
// Create a channel for hmm files.
ch_hmmrs = Channel.empty()
if ( params.hmmdir ) {
    Channel
        .fromPath(params.hmmdir + params.hmmpattern)
        .set { ch_hmmrs }
} else if ( params.hmmfiles ) {
    Channel
        .of( params.hmmfiles.split(',') )
        .map { [ file(it) ] }
        .set { ch_hmmrs }
}


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
//include { MERGE_TABLES                     } from '../modules/local/merge_summary_tables'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

//
// SUBWORKFLOW: Consisting of local modules
//
include { HMMCLASSIFY       } from '../subworkflows/local/hmmclassify'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ 	          	  } from '../modules/nf-core/cat/fastq/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

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


    // Step 3: Trim Galore!
    //


    // Step 4
    // Remove host sequences, bowtie2 align to Bos taurus
    // Subworkflow? bowtie2

    // Step 5
    // rRNA remove (sortmerna)
    //

    // Step 6
    // Deduplication with BBdup
    //

    // Step 7
    // Filter by taxa with Kraken2
    //



    // Step 10
    // Clustering with CD-HIT-EST
    // (for concatenating multiple assemblies)
    //

    // Step 11
    // Call ORFs
    //


    // Step 12
    // Quantification (salmon, rsem, or bbmap)
    //

    // Step 13
    // SUBWORKFLOW: classify ORFs with a set of hmm files
    //

    // Step 14
    // Kraken2 taxonomical annotation
    //

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
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
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
