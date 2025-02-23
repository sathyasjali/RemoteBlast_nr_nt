#!/usr/bin/env python
import csv
import sys

def parse_blast_tsv(input_tsv, output_csv):
    with open(input_tsv, 'r') as tsvfile, open(output_csv, 'w', newline='') as csvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        writer = csv.writer(csvfile)
        
        # Write header
        header = ["query_id", "subject_id", "percentage_identity", "alignment_length",
                  "mismatches", "gap_opens", "query_start", "query_end", "subject_start",
                  "subject_end", "e_value", "bit_score", "aligned_sequence"]
        writer.writerow(header)

        # Write data rows, excluding sequences with alignment length < 200bp or aligned sequence length < 200bp
        for row in reader:
            if row and int(row[3]) >= 200 and len(row[12]) >= 200:  # Check alignment_length and sseq
                writer.writerow(row)

if __name__ == "__main__":
    input_tsv = sys.argv[1]
    output_csv = sys.argv[2]
    parse_blast_tsv(input_tsv, output_csv)
