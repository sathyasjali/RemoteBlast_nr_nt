params.input       = "${baseDir}/test.fasta"
params.outdir      = "${baseDir}/results"
params.parseScript = "${baseDir}/bin/parse_blast_results.py"  // Parameterized path for the Python script

process blast_remote {

    label 'blast_container'

    input:
        path input_file

    output:
        tuple val("core_nt"), path("remote_blast_output_core_nt.tsv"), emit: ch_blast_core

    publishDir params.outdir, mode: 'copy'
    
    script:
    """
    set -e  # Exit on error if any command fails

    echo -e "Running Remote BLAST search using NCBI\\n"

    # Run BLAST against core_nt database
    echo "Adding headers and running BLAST search on core_nt..."
    echo -e "qseqid\tsseqid\tstitle\tpident\tqcovs\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tsseq" > remote_blast_output_core_nt.tsv
    blastn -query ${input_file} -db core_nt -remote \
        -outfmt "6 qseqid sseqid stitle pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore sseq" \
        >> remote_blast_output_core_nt.tsv \
        2> core_nt_error.log \
        || { echo "Error: core_nt BLAST failed. Check core_nt_error.log"; exit 1; }

    echo "..Done"
    """
}

process parse_blast_output {
    label 'python_container'

    input:
        path core_nt_file

    output:
        path "its.csv", emit: ch_parsed_results

    publishDir params.outdir, mode: 'copy'

    script:
    def parse_script = params.parseScript  // Define the variable in Groovy scope
    """
    set -e  # Exit on error if any command fails

    echo -e "Parsing Remote BLAST results\n"
    echo -e "Activating Conda environment: bio_env\n"

    source /opt/homebrew/Caskroom/miniconda/base/etc/profile.d/conda.sh  # Ensure Conda is available
    conda activate /opt/homebrew/Caskroom/miniconda/base/envs/bio_env  # Activate Conda

    python3 ${parse_script} ${core_nt_file} its.csv  # Correctly interpolated parse_script

    echo "..Parsing Done"
    """
}


workflow {
    // Define input channel from the FASTA file
    ch_inp = Channel.fromPath(params.input)

    // Run the BLAST remote process
    blast_remote_out = blast_remote(ch_inp)

    // Debugging: Check emitted output from blast_remote
    blast_remote_out.view()

    // Extract only the file path from the tuple (correcting the extraction issue)
    ch_blast_core = blast_remote_out.map { it[1] }

    // Pass the extracted file path to the parsing process
    parse_blast_output(ch_blast_core)
}
