#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.input = "${baseDir}/input.fasta"
params.outdir = "${baseDir}/results"

process blast_remote {

    label 'blast_container'
    input:
        path input
    output:
        tuple val("core_nt"), path("remote_blast_output_core_nt.tsv"), emit: ch_blast_core
        tuple val("tsa_nt"), path("remote_blast_output_tsa_nt.tsv"), emit: ch_blast_tsa
        //tuple path("remote_blast_output_core_nt.tsv"), emit: ch_blast_core
        //tuple path("remote_blast_output_tsa_nt.tsv"), emit: ch_blast_tsa
        

        publishDir "results", mode: 'copy'
    
    script:
    """
    set -e  # Exit on error if any command fails

    echo -e "Running Remote BLAST search using NCBI\n"

    # Run BLAST against core_nt database
    echo "Running BLAST search on core_nt..."
    blastn -query ${input} -db core_nt -remote \
        -out remote_blast_output_core_nt.tsv \
        -outfmt "6 qseqid sseqid stitle sseq pident qcovs" \
        2> core_nt_error.log \
        || { echo "Error: core_nt BLAST failed. Check core_nt_error.log"; exit 1; }

    # Run BLAST against tsa_nt database
    echo "Running BLAST search on tsa_nt..."
    blastn -query ${input} -db tsa_nt -remote \
        -out remote_blast_output_tsa_nt.tsv \
        -outfmt "6 qseqid sseqid stitle sseq pident qcovs" \
        2> tsa_nt_error.log \
        || { echo "Error: tsa_nt BLAST failed. Check tsa_nt_error.log"; exit 1; }

    echo "..Done"
    """
}

process parse_blast_output {
    label 'python_container'

    input:
        path core_nt_file
        path tsa_nt_file

    output:
        path "its.csv", emit: ch_parsed_results

    publishDir "results", mode: 'copy'

    script:
    '''
    set -e  # Exit on error if any command fails

    echo -e 'Parsing Remote BLAST results\\n'

    python3 bin/parse_blast_results.py ${core_nt_file} ${tsa_nt_file} its.csv

    echo '..Parsing Done'
    '''
}

workflow {
    // Define input channel
    ch_inp = Channel.fromPath(params.input)

    // Run BLAST remote process
    blast_remote_out = blast_remote(ch_inp)

}