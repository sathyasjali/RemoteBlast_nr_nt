#!/usr/bin/env python3

i#!/usr/bin/env python3

import pandas as pd
import sys

def parse_blast_results(core_nt_file, output_csv):
    """
    Parses BLAST output from the core_nt database, filters alignments shorter than 200 bp,
    and extracts the required columns.
    
    Parameters:
    - core_nt_file: Path to the BLAST output file for the core_nt database.
    - output_csv: Output CSV file to save the filtered results.
    """
    # Define column headers based on BLAST outfmt 6
    columns = ["qseqid", "sseqid", "stitle", "pident", "qcovs", "length", "mismatch", 
               "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sseq"]

    # Read the BLAST results for core_nt
    df_core = pd.read_csv(core_nt_file, sep="\t", names=columns)

    # Filter alignments with length less than 200 bp
    df_filtered = df_core[df_core["length"] >= 200]

    # Extract the required columns
    df_final = df_filtered[["qseqid", "sseqid", "stitle", "sseq"]]

    # Save the filtered results to a CSV file
    df_final.to_csv(output_csv, index=False)

    print(f"Filtered BLAST results saved to {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python parse_blast_results.py <core_nt_file> <output_csv>")
        sys.exit(1)

    core_nt_file = sys.argv[1]
    output_csv = sys.argv[2]
    
    parse_blast_results(core_nt_file, output_csv)
