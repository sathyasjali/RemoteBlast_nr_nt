#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input  = "${baseDir}/test.fasta"
params.outdir = "${baseDir}/results"
params.parseScript = "${baseDir}/bin/tile.py"  // Parameterized path for the Python script
params.condaProfile = "${baseDir}/my_env.yaml"  // Not required if using 'conda' directive


process blast_remote {
    label 'blast_container'
    
    input:
        path input_file

    output:
        // Emit two TSV files as a tuple (for example, for two BLAST searches)
        tuple val("core_nt"), path("remote_blast_output_core_nt.tsv"), emit: ch_blast_core
        tuple val("tsa_nt"),  path("remote_blast_output_tsa_nt.tsv"),  emit: ch_blast_tsa

    publishDir params.outdir, mode: 'copy'
    
    script:
    """
    set -e

    echo -e "Running Remote BLAST search using NCBI\\n"

    # BLAST for core_nt database
    echo "Adding headers and running BLAST search on core_nt..."
    echo -e "qseqid\tsseqid\tstitle\tpident\tqcovs\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tsseq" > remote_blast_output_core_nt.tsv
    blastn -query ${input_file} -db core_nt -remote \\
         -outfmt "6 qseqid sseqid stitle pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore sseq" \\
         >> remote_blast_output_core_nt.tsv \\
         2> core_nt_error.log || { echo "Error: core_nt BLAST failed"; exit 1; }

    # BLAST for tsa_nt database
    echo "Running BLAST search on tsa_nt (with header)..."
    echo -e "qseqid\tsseqid\tstitle\tpident\tqcovs\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tsseq" > remote_blast_output_tsa_nt.tsv
    blastn -query ${input_file} -db tsa_nt -remote \\
    -outfmt "6 qseqid sseqid stitle pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore sseq" \\
         >> remote_blast_output_tsa_nt.tsv \\
         2> tsa_nt_error.log || { echo "Error: tsa_nt BLAST failed"; exit 1; }

    echo "..Done"
    """
}

process combineBlast {

    publishDir params.outdir, mode: 'copy'
    label "blast_container"

    // Input: pass the output directory containing tsv files.
    input:
      path(input, stageAs: "/*") // use this for 2 levels down "?/*" or "/*" for 1 level

    output:
      path "Blast_hits.tsv", emit: tsv

  
    script:
    """
    echo -e "Aggregating BLAST hits\n"
    
    # Extract header from the core_nt file and write it to Blast_hits.tsv
    echo "Extracting header from results/remote_blast_output_core_nt.tsv..."
    head -n 1 results/remote_blast_output_core_nt.tsv > Blast_hits.tsv

    # Now, for all TSV files in the results folder, skip the header and concatenate the rest.
    echo "Appending TSV contents (skipping duplicate headers)..."
    ls results/*nt.tsv | parallel -j1 "tail -n +2 {}" >> Blast_hits.tsv
    """
}

process parse_blast_output {
    conda 'my_env.yaml'  // This will activate the environment defined in my_env.yaml located in the base directory

    input:
        path tsv_file

    output:
        path "its.csv", emit: ch_parsed_results

    publishDir params.outdir, mode: 'copy'

    script:
    """
    set -e
    echo -e "Parsing and tiling BLAST results\n"
    
    # If you need to activate a conda profile, you can source it if required.
    # source ${params.condaProfile}  (optional if already activated by Nextflow via conda)
    
    python3 ${params.parseScript} ${tsv_file} its.csv
    
    echo "..Parsing Done"
    """
}


workflow {
    // Create an input channel from the FASTA file.
    // This channel will be used to pass the input file to the blast_remote process.
    ch_inp = Channel.fromPath(params.input)
    
    // Run the blast_remote process (this will publish TSV files into the results folder).
    remote = blast_remote(ch_inp)
    
     // Run combineBlast to merge the TSV files. This produces Blast_hits.tsv.
    merged_tsv = combineBlast( file(params.outdir) )
    merged_tsv.view()
    
    // Parse the merged TSV file to generate tiled output as a CSV file.
    parse_blast_output(merged_tsv)
    
}