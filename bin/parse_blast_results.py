import pandas as pd
import argparse

def parse_blast_results(input_file, output_file):
    """
    Parses the BLAST results CSV and extracts required columns into a new CSV file.

    :param input_file: Path to the input CSV file
    :param output_file: Path to the output formatted CSV file
    """
    # Define essential columns required in the final output
    required_columns = [
        "query", "accession", "title", "sequence", "identity", "alignment_length",
        "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "e_value", "bit_score"
    ]

    # Read the CSV file
    try:
        df = pd.read_csv(input_file)

        # Ensure all required columns exist
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            print(f"Warning: Missing columns in input file: {missing_columns}. Proceeding with available data.")

        # Select only existing required columns
        df_selected = df[[col for col in required_columns if col in df.columns]]

        # Save the formatted file
        df_selected.to_csv(output_file, index=False)
        print(f"Successfully saved parsed BLAST results to {output_file}")

    except Exception as e:
        print(f"Error processing file: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse BLAST result CSV file and format output")
    parser.add_argument("--input", required=True, help="Path to the input CSV file")
    parser.add_argument("--output", required=True, help="Path to the output CSV file")
    
    args = parser.parse_args()
    
    parse_blast_results(args.input, args.output)
