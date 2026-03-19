#!/usr/bin/env python3

import argparse
import pandas as pd
import os

def parse_blast(m8_file, fasta_file, refdb, max_hits, split_tax, tax_columns, output_prefix, threads, export_excel):
    """
    Parses a BLAST output file (.m8) and generates summary reports.
    """
    print(f"Parsing BLAST results from {m8_file}...")

    # Define column names based on BLAST output format 6
    columns = [
        "query_id", "subject_id", "percent_identity", "alignment_length", "mismatches",
        "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score", "sequence"
    ]

    # Read the BLAST output file
    df = pd.read_csv(m8_file, sep="\t", names=columns, header=None)

    # Compute query coverage
    df["query_coverage"] = (df["alignment_length"] / (df["q_end"] - df["q_start"] + 1)) * 100

    # Placeholder description (since not available in BLAST output)
    df["description"] = "N/A"

    # Select required columns
    parsed_df = df[["query_id", "subject_id", "description", "percent_identity", "query_coverage", "evalue", "sequence"]]

    # Save results in different formats
    wide_output = f"{output_prefix}_wide.csv"
    long_output = f"{output_prefix}_long.csv"
    best_hits_xlsx = f"{output_prefix}_BestHits.xlsx"
    best_hits_txt = f"{output_prefix}_BestHits.txt.gz"

    # Save in wide format
    parsed_df.to_csv(wide_output, index=False)
    
    # Save in long format
    long_df = parsed_df.melt(id_vars=["query_id", "subject_id"], value_vars=["percent_identity", "query_coverage", "evalue", "sequence"])
    long_df.to_csv(long_output, index=False)

    # Save best hits (sorted by highest identity)
    best_hits = parsed_df.sort_values(by=["percent_identity"], ascending=False).groupby("query_id").head(max_hits)
    
    # Save as Excel (if requested)
    if export_excel:
        best_hits.to_excel(best_hits_xlsx, index=False)

    # Save as compressed text file
    best_hits.to_csv(best_hits_txt, sep="\t", index=False, compression="gzip")

    print(f"Parsed BLAST results saved:\n - {wide_output}\n - {long_output}\n - {best_hits_xlsx} (optional)\n - {best_hits_txt} (optional)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse BLAST results and generate summary reports.")
    
    parser.add_argument("--m8", required=True, help="BLAST output file (.m8)")
    parser.add_argument("--fasta", required=True, help="Query FASTA file")
    parser.add_argument("--db", required=True, help="Reference database")
    parser.add_argument("--maxhits", type=int, default=10, help="Maximum number of hits to retain")
    parser.add_argument("--splittax", action="store_true", help="Split taxonomy data")
    parser.add_argument("--taxcolumns", type=str, default="", help="Taxonomy columns to extract")
    parser.add_argument("--outputprefix", default="Blast_hits", help="Output file prefix")
    parser.add_argument("--threads", type=int, default=4, help="Number of processing threads")
    parser.add_argument("--exportexcel", action="store_true", help="Export results as Excel file")

    args = parser.parse_args()

    parse_blast(
        m8_file=args.m8,
        fasta_file=args.fasta,
        refdb=args.db,
        max_hits=args.maxhits,
        split_tax=args.splittax,
        tax_columns=args.taxcolumns,
        output_prefix=args.outputprefix,
        threads=args.threads,
        export_excel=args.exportexcel
    )
