#!/usr/bin/env python3

import pandas as pd
import sys

def parse_blast_results(core_nt_file, output_csv):
    """
    Parses BLAST output from the core_nt database, filters alignments shorter than 200 bp,
    and extracts the required columns with descriptive labels.

    Parameters:
    - core_nt_file: Path to the BLAST output file for the core_nt database.
    - output_csv: Output CSV file to save the filtered results.
    """
    # Define column headers based on BLAST outfmt 6
    columns = [
        "Query_ID", "Subject_ID", "Subject_Title", "Percent_Identity", 
        "Query_Coverage", "Alignment_Length", "Mismatches", "Gap_Openings", 
        "Query_Start", "Query_End", "Subject_Start", "Subject_End", "E_Value", 
        "Bit_Score", "Subject_Sequence"
    ]

    # Read the BLAST results and skip the first row (header)
    df_core = pd.read_csv(core_nt_file, sep="\t", names=columns, skiprows=1)

    # Convert "Alignment_Length" column to integer
    df_core["Alignment_Length"] = pd.to_numeric(df_core["Alignment_Length"], errors="coerce")

    # Drop rows where "Alignment_Length" could not be converted
    df_core = df_core.dropna(subset=["Alignment_Length"])

    # Convert back to integer
    df_core["Alignment_Length"] = df_core["Alignment_Length"].astype(int)

    # Filter alignments with length less than 200 bp
    df_filtered = df_core[df_core["Alignment_Length"] >= 200]

    # Extract the required columns with descriptive names
    df_final = df_filtered[["Query_ID", "Subject_ID", "Subject_Title", "Subject_Sequence"]]

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
