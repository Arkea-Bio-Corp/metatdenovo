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
include { HMMER_HMMSCAN                    } from '../modules/local/hmmscan/main'
include { EGGNOG_MAPPER                    } from '../modules/local/eggnog/mapper'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

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

include { CAT_CAT                      } from '../modules/nf-core/cat/cat/'
include { CAT_CAT as CAT_QUANT         } from '../modules/nf-core/cat/cat/'
include { CAT_FASTQ 	          	   } from '../modules/nf-core/cat/fastq/'
include { MULTIQC                      } from '../modules/nf-core/multiqc/'
include { CUSTOM_DUMPSOFTWAREVERSIONS  } from '../modules/nf-core/custom/dumpsoftwareversions/'
include { CDHIT_CDHIT                  } from '../modules/nf-core/cdhit/'
include { KRAKEN2_KRAKEN2 as KRKN_ARCH } from '../modules/nf-core/kraken2/'
include { SALMON_INDEX                 } from '../modules/nf-core/salmon/index/'
include { SALMON_QUANT                 } from '../modules/nf-core/salmon/quant/'
include { SALMON_MERGE                 } from '../modules/nf-core/salmon/merge/'
include { TRANSDECODER_LONGORF         } from '../modules/nf-core/transdecoder/longorf/'
include { TRANSDECODER_PREDICT         } from '../modules/nf-core/transdecoder/predict/'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def multiqc_report = []


workflow POST_ASSEMBLE_CLUSTER {

    ch_versions = Channel.empty()
    ch_read_counts = Channel.empty()

    // Gather files
    // Sample sheet csv should be organized like so:
    // sample,assembly,reads
    // sample1,/path/to/assembly.fa.gz,/path/to/reads.fq.gz

    Channel.fromPath(ch_input)
        .splitCsv(header: ['sample', 'assembly', 'reads'], skip: 1 )
        .collect(row -> "${row.assembly}")
        .map { [[id: "combined_assembly", single_end: true], it] }
        .set { assembly_list }

    Channel.fromPath(ch_input)
        .splitCsv(header: ['sample', 'assembly', 'reads'], skip: 1 )
        .map { row -> [[id: row.sample, single_end: true], row.reads] }
        .set { reads_list }
    
    // combine gzipped fasta assemblies with CAT 🐈
    CAT_CAT(assembly_list)
    ch_versions = ch_versions.mix(CAT_CAT.out.versions)

    // cluster by sequence similarity
    CDHIT_CDHIT(CAT_CAT.out.file_out)
    ch_versions = ch_versions.mix(CDHIT_CDHIT.out.versions) 

    // ORF prediction
    tdecoder_folder = TRANSDECODER_LONGORF(CDHIT_CDHIT.out.fasta).folder
    ch_versions = ch_versions.mix(TRANSDECODER_LONGORF.out.versions)
    TRANSDECODER_PREDICT(CDHIT_CDHIT.out.fasta, tdecoder_folder)
    ch_versions = ch_versions.mix(TRANSDECODER_PREDICT.out.versions)

    // Quantification w/ salmon
    salmon_ind = SALMON_INDEX(CAT_CAT.out.file_out).index
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
    SALMON_QUANT(reads_list, salmon_ind)   
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)

    SALMON_MERGE(SALMON_QUANT.out.quants.collect())
    ch_versions = ch_versions.mix(SALMON_MERGE.out.versions)

    // Functional annotation with eggnog-mapper
    // eggdbchoice = ["diamond", "mmseqs", "hmmer", "novel_fams"]
    // eggnog_ch = Channel.value(file(params.eggnogdir, checkIfExists: true))
    // EGGNOG_MAPPER(TRANSDECODER_PREDICT.out.pep, eggnog_ch, eggdbchoice)
    // ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions)

    // Functional annotation with hmmscan
    hmmerdir   = Channel.fromPath(params.hmmdir, checkIfExists: true)
    hmmerfile  = Channel.value(params.hmmerfile)
    HMMER_HMMSCAN(TRANSDECODER_PREDICT.out.pep, hmmerdir, hmmerfile)
    ch_versions = ch_versions.mix(HMMER_HMMSCAN.out.versions)
 
    // Kraken2 taxonomical annotation of contigs 
    // adjust metamap to switch to single end
    trans_cds = TRANSDECODER_PREDICT.out.cds.map{ [[id: it[0].id, single_end: true], it[1]]} 
    k2_arch_db = Channel.fromPath(params.archaea_db, checkIfExists: true)
    KRKN_ARCH(trans_cds, k2_arch_db, true, true)
    ch_versions = ch_versions.mix(KRKN_ARCH.out.versions)
}