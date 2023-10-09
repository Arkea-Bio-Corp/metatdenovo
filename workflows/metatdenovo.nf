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
include { CAT_FASTQ as CAT_FASTQ_SALM      } from '../modules/nf-core/cat/fastq/'
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
include { SALMON_MERGE                     } from '../modules/nf-core/salmon/merge/'
include { SORTMERNA                        } from '../modules/nf-core/sortmerna/'
include { TRANSDECODER_LONGORF             } from '../modules/nf-core/transdecoder/longorf/'
include { TRANSDECODER_PREDICT             } from '../modules/nf-core/transdecoder/predict/'
include { SEQKIT_SPLIT2                    } from '../modules/nf-core/seqkit/split2/'
include { TRIMMOMATIC                      } from '../modules/nf-core/trimmomatic'
include { TRINITY                          } from '../modules/nf-core/trinity/'
include { SOAP_DENOVO_TRANS                } from '../modules/local/soap-denovo-trans'
include { MEGAHIT                          } from '../modules/nf-core/megahit/'
include { TRANS_ABYSS                      } from '../modules/local/transabyss'
include { ASSEMBLE_STATS                   } from '../subworkflows/local/bt2_assembly_stats'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow METATDENOVO {

    ch_versions = Channel.empty()
    ch_read_counts = Channel.empty()

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

    // FastQC
    // 
    PRE_TRIM_FQC(ch_fastq[0], "PRE_TRIM")
    ch_versions = ch_versions.mix(PRE_TRIM_FQC.out.versions)


    // Split by split_size
    // ch_fastq[0]
    //     .map { [it.get(0), it.get(1)[0], it.get(1)[1]] }
    //     .splitFastq(by: params.split_size, 
    //                 pe: true, file: false, 
    //                 compress: false, decompress: true)
    //     .map{ [ it.get(0), [it.get(1), it.get(2)] ]}
    //     .set { trim_split }

    SEQKIT_SPLIT2(ch_fastq[0])
    SEQKIT_SPLIT2.out.reads
        .map { meta, reads ->
            int leng = reads.size()/2
            splitreads = reads.collate(leng)
            def combined_arr = []
            for(int i = 0;i<leng;i++) {
                combined_arr.add([splitreads[0][i], splitreads[1][i]])
            }
            [meta, combined_arr]
        }
        .transpose()
        .set { trim_split }

    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)
    //
    // trimmomatic
    //
    adapter_path = Channel.value(file(params.adapter_fa, checkIfExists: true))
    TRIMMOMATIC(trim_split, adapter_path)
    TRM_LOGS(
        TRIMMOMATIC.out.collect_log.collectFile(name: 'collected_logs.txt', newLine: true),
        "Trimmomatic",
        TRIMMOMATIC.out.meta
    )
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)
    ch_read_counts = ch_read_counts.mix(TRIMMOMATIC.out.readcounts)

    // 
    // Remove host sequences, bowtie2 align to Bos taurus
    // 
    index_ch = Channel.value(file(params.indexdir, checkIfExists: true))
    BOWTIE2_ALIGN(TRIMMOMATIC.out.trimmed_reads, index_ch, true, false)
    BT2_LOGS(
        BOWTIE2_ALIGN.out.collect_log.collectFile(name: 'collected_logs.txt', newLine: true),
        "Bowtie2",
        BOWTIE2_ALIGN.out.meta
    )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)
    ch_read_counts = ch_read_counts.mix(BOWTIE2_ALIGN.out.readcounts)

    // 
    // rRNA remove (bowtie)
    // 
    silva_ch = Channel.value(file(params.silva_reference, checkIfExists: true))
    rna_idx  = Channel.value(file(params.rna_idx, checkIfExists: true))
    SORTMERNA(BOWTIE2_ALIGN.out.fastq, silva_ch, rna_idx)
    SMR_LOGS(
        SORTMERNA.out.collect_log.collectFile(name: 'collected_logs.txt', newLine: true),
        "SortMeRNA",
        SORTMERNA.out.meta
    )
    ch_versions = ch_versions.mix(SORTMERNA.out.versions)
    ch_read_counts = ch_read_counts.mix(SORTMERNA.out.readcounts)

    // 
    // Filter by taxa with Kraken2
    // 
    k2db_ch = Channel.value(file(params.no_archaea_db, checkIfExists: true))
    KRKN_NO_ARCH(SORTMERNA.out.reads, k2db_ch, true, true)
    // KR2_LOGS(
    //     KRKN_NO_ARCH.out.collect_log.collectFile(name: 'collected_logs.txt', newLine: true),
    //     "Kraken2_no_archaea",
    //     KRKN_NO_ARCH.out.meta
    // )
    ch_versions = ch_versions.mix(KRKN_NO_ARCH.out.versions)
    ch_read_counts = ch_read_counts.mix(KRKN_NO_ARCH.out.readcounts)

    // RECOMBINE ~~~~
    KRKN_NO_ARCH.out.unclassified_reads_fastq
        .groupTuple()
        .collect ( sample -> sample[1].flatten() )
        .map { [[id: "all_samples", single_end: false], it] }
        .set { collectedFastqs }
    CAT_FASTQ(collectedFastqs)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)
    ch_read_counts = ch_read_counts.mix(CAT_FASTQ.out.readcounts)

    //
    // merge paired reads into single end
    //
    BBMAP_MERGE(CAT_FASTQ.out.reads)
    merged_reads = BBMAP_MERGE.out.merged.map{ [[id: it[0].id, single_end: true], it[1]]}
    ch_versions = ch_versions.mix(BBMAP_MERGE.out.versions)
    ch_read_counts = ch_read_counts.mix(BBMAP_MERGE.out.readcounts)

    // run merged reads through FastQC 
    POST_MERGE_FQC(merged_reads, "POST_MERGE")
    ch_versions = ch_versions.mix(POST_MERGE_FQC.out.versions)

    // 
    // Deduplication with Dedupe
    // 
    BBMAP_DEDUPE(merged_reads)
    ch_versions = ch_versions.mix(BBMAP_DEDUPE.out.versions)
    ch_read_counts = ch_read_counts.mix(BBMAP_DEDUPE.out.readcounts)
    ch_collect_counts = ch_read_counts.collectFile(name: 'collated_counts.csv')
    CUSTOM_DUMPCOUNTS(ch_collect_counts, BBMAP_DEDUPE.out.meta)

    // 
    // Run your assembler of choice
    // 
    if (params.assembler == "Megahit") {
        MEGAHIT(BBMAP_DEDUPE.out.reads)
        ch_versions = ch_versions.mix(MEGAHIT.out.versions)
        assembled_contigs = MEGAHIT.out.contigs
    } else if (params.assembler == "Trans-Abyss") {
        TRANS_ABYSS(BBMAP_DEDUPE.out.reads)
        ch_versions = ch_versions.mix(TRANS_ABYSS.out.versions)
        assembled_contigs = TRANS_ABYSS.out.transcripts
    } else if (params.assembler == "SOAP-DeNovo-Trans") {
        SOAP_DENOVO_TRANS(BBMAP_DEDUPE.out.reads)
        ch_versions = ch_versions.mix(SOAP_DENOVO_TRANS.out.versions)
        assembled_contigs = SOAP_DENOVO_TRANS.out.contigs
    } else { // Trinity case
        TRINITY(BBMAP_DEDUPE.out.reads)
        ch_versions = ch_versions.mix(TRINITY.out.versions)
        assembled_contigs = TRINITY.out.transcript_fasta
    }

    //
    // Run a quick QC check
    // & change meta tag to improve output directory labeling
    assembled_contigs
        .map { [[id: it[0].id, 
                single_end: it[0].single_end, 
                assembler: params.assembler], 
            it[1]] }
        .set { qc_contigs }
    ASSEMBLE_STATS(qc_contigs, BBMAP_DEDUPE.out.reads, params.assembler)

    // Choose to run post assembly steps, can save a bit of time if you only need the 
    // assemblies. Doesn't interact with any MultiQC inputs.
    if (params.run_post_assembly) {
        // 
        // Clustering with CD-HIT-EST to remove redundancies
        // 
        CDHIT_CDHIT(assembled_contigs)
        ch_versions = ch_versions.mix(CDHIT_CDHIT.out.versions)

        // 
        // ORF prediction and translation to peptide seqs with Transdecoder
        // 
        tdecoder_folder = TRANSDECODER_LONGORF(CDHIT_CDHIT.out.fasta).folder
        ch_versions = ch_versions.mix(TRANSDECODER_LONGORF.out.versions)
        TRANSDECODER_PREDICT(CDHIT_CDHIT.out.fasta, tdecoder_folder)
        ch_versions = ch_versions.mix(TRANSDECODER_PREDICT.out.versions)

        // 
        // Quantification w/ salmon
        //
        // need to merge these reads again, by sample
        KRKN_NO_ARCH.out.unclassified_reads_fastq
            .groupTuple()
            .map { [it[0], it[1].flatten()] }
            .set { collectedFastqs }
        CAT_FASTQ_SALM(collectedFastqs)
        ch_versions = ch_versions.mix(CAT_FASTQ_SALM.out.versions)

        salmon_ind = SALMON_INDEX(assembled_contigs).index
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
        SALMON_QUANT(CAT_FASTQ_SALM.out.reads, salmon_ind)   
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)

        SALMON_QUANT.out.quants
            .toList()
            .transpose()
            .toList()
            .set { quant_list }
        SALMON_MERGE(quant_list)
        ch_versions = ch_versions.mix(SALMON_MERGE.out.versions)

        // 
        // Functional annotation with eggnog-mapper
        // 
        eggdbchoice = ["diamond", "mmseqs", "hmmer", "novel_fams"]
        eggnog_ch = Channel.fromPath(params.eggnogdir, checkIfExists: true)
        EGGNOG_MAPPER(TRANSDECODER_PREDICT.out.pep, eggnog_ch, eggdbchoice)
        ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions)

        // 
        // Functional annotation with hmmscan
        // 
        hmmerdir   = Channel.fromPath(params.hmmdir, checkIfExists: true)
        hmmerfile  = Channel.value(params.hmmerfile)
        HMMER_HMMSCAN(TRANSDECODER_PREDICT.out.pep, hmmerdir, hmmerfile)
        ch_versions = ch_versions.mix(HMMER_HMMSCAN.out.versions)

        // 
        // Kraken2 taxonomical annotation of contigs
        // 
        // adjust metamap to switch to single end
        trans_cds = TRANSDECODER_PREDICT.out.cds.map{ [[id: it[0].id, single_end: true], it[1]]} 
        k2_arch_db = Channel.fromPath(params.archaea_db, checkIfExists: true)
        KRKN_ARCH(trans_cds, k2_arch_db, true, true)
        ch_versions = ch_versions.mix(KRKN_ARCH.out.versions)
    }


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml'),
        BBMAP_DEDUPE.out.meta
    )

    // 
    // MultiQC
    //
    workflow_summary    = WorkflowMetatdenovo.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    COUNTS_PLOT(CUSTOM_DUMPCOUNTS.out.readcounts)
    ch_versions = ch_versions.mix(COUNTS_PLOT.out.versions)

    methods_description    = WorkflowMetatdenovo.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(PRE_TRIM_FQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(POST_MERGE_FQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(COUNTS_PLOT.out.counts_png.ifEmpty([]))
    // Plotly output is buggy, skip for now:
    // ch_multiqc_files = ch_multiqc_files.mix(COUNTS_PLOT.out.counts_html.ifEmpty([]))


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        BBMAP_DEDUPE.out.meta
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END :^)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
 