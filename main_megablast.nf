#!/usr/bin/env nextflow
// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Pipeline version
version = '0.1'

// Input and output direcotries
params.input  = "${launchDir}/test1.fasta"
params.outdir = "${launchDir}/results"

params.method = "megablast"                // other methods not implemented yet

// Taxonomy annotation with BLAST
params.blast_taxdb        = "/Users/sathyajali/Documents/BatchBlaster/blast_db/headtest"
params.blastdb_fasta      = "/Users/sathyajali/Documents/BatchBlaster/blast_db/headtest.fasta"
params.blast_task         = "megablast"   // or "megablast"
params.blast_chunksize    = 4000
params.blast_maxts        = 10
params.blast_hsps         = 1
params.blast_wordsize     = 7
params.blast_evalue       = 0.001
params.blast_reward       = 1
params.blast_penalty      = -1
params.blast_gapopen      = 1
params.blast_gapextend    = 2
params.blast_percidentity = 80
//params.variant_size      = 200

// BLAST-parsing
params.blast_splittax     = true
params.blast_taxcolumns   = "AccID,Kingdom,Phylum,Class,Order,Family,Genus,Species,Function,Dataset"
params.exportexcel        = true


if(params.blast_taxdb){
  bastdb_name = file(params.blast_taxdb).name
  bastdb_dir = file(params.blast_taxdb).parent
}

// Check if input path was provided
if (params.input == false) {
  println( "Please provide the input file with sequences in FASTA format with `--input` parameter.")
  exit(1)
}
// Check if BLAST database was provided
if(params.blast_taxdb == false && params.method == "blast") {
  println( "Please provide the BLAST database with `--blast_taxdb` parameter.")
  exit(1)
}


// Pipeline help message
def helpMsg() {
    log.info"""
    =====================================================================
    BatchBlaster ${version}
    =====================================================================
    
    Pipeline Usage:
    To run the pipeline, enter the following in the command line:
        nextflow run vmikk/batchblaster -r "main" --input ... --outdir ...
    
    If running on HPC, adjust the max heap size of the Java VM and specify the profile:
        export NXF_OPTS="-Xms500M -Xmx2G"
        nextflow run vmikk/batchblaster -r "main" -profile cluster,singularity ...

    Options:
    REQUIRED:
        --input               Input file with sequences (FASTA)
        --outdir              The output directory where the results will be saved

    OPTIONAL:
        --blast_taxdb         BLAST database
        --blast_task          Task to execute (`blastn` or `megablast`)
        --blast_chunksize     Number of sequences per analysis chunk
        --blast_maxts         Maximum number of aligned sequences to keep
        --blast_hsps          Maximum number of HSPs per subject sequence to save for each query
        --blast_wordsize      Word size for wordfinder algorithm
        --blast_evalue        Expectation value (E) threshold for saving hits (default, 10)
        --blast_reward        Reward for a nucleotide match
        --blast_penalty       Penalty for a nucleotide mismatch
        --blast_gapopen       Cost to open a gap
        --blast_gapextend     Cost to extend a gap
        --blast_percidentity  Percent identity
        --blast_splittax      Split headers of database records into columns in the output (default, true). In the header, semicolon should be used as a separator
        --blast_taxcolumns    Column names for header parts of the database records (parameter should be defined as a comma-separated list)

    NEXTFLOW-SPECIFIC:
        -qs                   Queue size (max number of processes that can be executed in parallel); e.g., 8
        -resume               Resume the workflow if it was stopped by an error (execution will continue using the cached results)

    """.stripIndent()
}
// Show help msg
if (params.help){
    helpMsg()
    exit(0)
}


// Print the parameters to the console and to the log
log.info """
    =======================================================================
    BatchBlaster ${version}
    =======================================================================
    Input data path: ${params.input}
    Output path:     ${params.outdir}
    Method:          ${params.method}
    Database:        ${params.blast_taxdb}
    """
    .stripIndent()

// log.info """
//         Pipeline info:
//           Pipeline profile:       ${workflow.profile}
//           Config file used:       ${workflow.configFiles}
//           Container engine:       ${workflow.containerEngine}
//         """
//         .stripIndent()
// 
// log.info """
//         Core Nextflow options:
//           launchDir:              ${workflow.launchDir}
//           workDir:                ${workflow.workDir}
//           projectDir:             ${workflow.projectDir}
//         """
//         .stripIndent()
// 
// log.info "======================================================================="
log.info "\n"



// Taxonomy annotation
// use globally-dereplicated sequences
process blast {
    publishDir "results", mode: 'copy'
    label "main_container"

    // publishDir params.outdir, mode: 'copy'
    // cpus 10

    input:
      path input
      path taxdb_dir

    output:
      path "${input.getBaseName()}.m8.gz", emit: blastchunks

    script:

    // Optional parameters
    wordsize  = params.blast_wordsize     ? "-word_size ${params.blast_wordsize}"         : ""
    evalue    = params.blast_evalue       ? "-evalue    ${params.blast_evalue}"           : ""
    reward    = params.blast_reward       ? "-reward    ${params.blast_reward}"           : ""
    penalty   = params.blast_penalty      ? "-penalty   ${params.blast_penalty}"          : ""
    gapopen   = params.blast_gapopen      ? "-gapopen   ${params.blast_gapopen}"          : ""
    gapextend = params.blast_gapextend    ? "-gapextend ${params.blast_gapextend}"        : ""
    percidentity = params.blast_percidentity ? "-perc_identity ${params.blast_percidentity}" : ""
    
    """
    echo -e "Taxonomy annotation with BLAST\n"
    echo -e "Input file: " ${input}
    echo -e "Database  : " ${taxdb_dir}

    blastn \
      -task ${params.blast_task} \
      -query ${input} \
      -db ${taxdb_dir}/${bastdb_name} \
      -strand both \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" \
      -max_target_seqs ${params.blast_maxts} \
      -max_hsps ${params.blast_hsps} \
      -out ${input.getBaseName()}.m8 \
      -num_threads ${task.cpus} \
      ${wordsize} \
      ${evalue} \
      ${reward} \
      ${penalty} \
      ${gapopen} \
      ${gapextend} \
      ${percidentity} \

    ## Compress results
    gzip -7 ${input.getBaseName()}.m8

    echo "..Done"

    """
}


// Aggregate BLAST results
process blast_merge {
    publishDir "results", mode: 'copy'

    label "main_container"

    publishDir params.outdir, mode: 'copy'
    // cpus 1

    input:
      path(input, stageAs: "?/*")

    output:
      path "Blast_hits.m8.gz", emit: m8

    script:
    """

    echo -e "Aggregating BLAST hits\n"
    
    find . -mindepth 1 -maxdepth 2 -name "*.m8.gz" \
      | parallel -j1 "cat {}" \
      >> Blast_hits.m8.gz

    echo "..Done"

    """
}

// Parse BLAST results
process parse_blast {

    label "parse_results"

    publishDir "results", mode: 'copy'
    // cpus 4

    input:
      path m8
      path fasta
      path refdb

    output:
      path "Blast_hits_wide.csv",      emit: Bwide
      path "Blast_hits_long.csv",      emit: Blong
      path "Blast_hits_BestHits.xlsx", emit: Bxlsx, optional: true
      path "Blast_hits_BestHits.txt.gz", emit: Btsv, optional: true

    script:
    """
    echo -e "Parsing BLAST results using container\n"
    
    pip install openpyxl

    python3 /Users/sathyajali/Documents/BatchBlaster/bin/parse_blast.py \
      --m8           ${m8} \
      --fasta        ${fasta} \
      --db           ${refdb} \
      --maxhits      ${params.blast_maxts} \
      ${params.blast_splittax ? "--splittax" : ""} \
      --taxcolumns   "${params.blast_taxcolumns}" \
      --outputprefix Blast_hits \
      --threads      ${task.cpus} \
      ${params.exportexcel ? "--exportexcel" : ""}

    echo "..Done"
    """
}


workflow {

  // Input file with sequences (FASTA)
  ch_inp = Channel.fromPath(params.input)

  // Input file with sequences (FASTA)
  //ch_fasta = Channel.fromPath(params.input)

  // Split FASTA sequences into multiple chunks
  ch_fasta = ch_inp.splitFasta(by: params.blast_chunksize, file:true, compress:false)

  // Run taxonomy annotation
  blast(ch_fasta, blastdb_dir)

  // Aggregate BLAST results
  ch_blasthits = blast.out.blastchunks.collect()
  blast_merge(ch_blasthits)

  // Parse BLAST results
  ch_blastdbfasta = Channel.fromPath(params.blastdb_fasta)
  parse_blast(
    blast_merge.out.m8,
    ch_inp,
    ch_blastdbfasta)

}


// On completion
workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Duration              : ${workflow.duration}"
    println "Execution status      : ${workflow.success ? 'All done!' : 'Failed' }"
}

// On error
workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
