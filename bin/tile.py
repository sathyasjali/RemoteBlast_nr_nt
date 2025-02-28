#!/usr/bin/env python3
import sys
import pandas as pd

# Tiling parameters
TILE_LENGTH = 200
OVERLAP = 75
STEP = TILE_LENGTH - OVERLAP  # 125 bp step ensures 75 bp overlap

def parse_and_tile_blast_results(tsv_file, output_csv):
    """
    Reads a BLAST TSV file, filters for alignments with length >= 200 bp, 
    and for each row generates tiling windows from the 'Subject_Sequence' column.
    
    Each window is TILE_LENGTH bp long and each tile overlaps the previous one by OVERLAP bp.
    All tiles for a given row are written in separate columns (Tile_1, Tile_2, etc.)
    in the output CSV.
    
    The output CSV will have the following columns:
      - Query_ID, Subject_ID, Subject_Title, Tile_1, Tile_2, Tile_3, ...
    """
    # Define expected column names (adjust if needed)
    columns = [
        "Query_ID", "Subject_ID", "Subject_Title", "Percent_Identity", 
        "Query_Coverage", "Alignment_Length", "Mismatches", "Gap_Openings", 
        "Query_Start", "Query_End", "Subject_Start", "Subject_End", "E_Value", 
        "Bit_Score", "Subject_Sequence"
    ]
    
    # Read the TSV file (assuming the first line is the header)
    df = pd.read_csv(tsv_file, sep="\t", header=0, names=columns)
    
    # Convert Alignment_Length to numeric, drop rows with invalid values, and filter for >=200 bp
    df["Alignment_Length"] = pd.to_numeric(df["Alignment_Length"], errors="coerce")
    df = df.dropna(subset=["Alignment_Length"])
    df["Alignment_Length"] = df["Alignment_Length"].astype(int)
    df_filtered = df[df["Alignment_Length"] >= 200]
    
    aggregated_rows = []
    max_tiles = 0

    for _, row in df_filtered.iterrows():
        seq = row["Subject_Sequence"]
        seq_len = len(seq)
        
        if seq_len < TILE_LENGTH:
            # If the sequence is shorter than 200bp, output it as a single tile.
            tiles = [seq]
        else:
            # Generate overlapping tiles with a step of STEP (125 bp)
            tiles = [seq[i:i+TILE_LENGTH] for i in range(0, seq_len - TILE_LENGTH + 1, STEP)]
        
        if len(tiles) > max_tiles:
            max_tiles = len(tiles)
        
        # Create a dictionary with the BLAST identifiers and the tiles
        row_dict = {
            "Query_ID": row["Query_ID"],
            "Subject_ID": row["Subject_ID"],
            "Subject_Title": row["Subject_Title"]
        }
        for idx, tile in enumerate(tiles, start=1):
            row_dict[f"Tile_{idx}"] = tile
        
        aggregated_rows.append(row_dict)
    
    # Create a list of all tile column names up to the maximum number encountered
    tile_cols = [f"Tile_{i}" for i in range(1, max_tiles+1)]
    base_cols = ["Query_ID", "Subject_ID", "Subject_Title"]
    all_cols = base_cols + tile_cols

    # Create the final DataFrame and ensure all expected tile columns are present
    df_tiles = pd.DataFrame(aggregated_rows)
    for col in tile_cols:
        if col not in df_tiles.columns:
            df_tiles[col] = ""
    df_tiles = df_tiles[all_cols]  # Reorder columns

    # Save the output as a CSV file
    df_tiles.to_csv(output_csv, index=False)
    print(f"Filtered and tiled BLAST results saved to {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python3 parse_and_tile_blast_results.py <input_tsv> <output_csv>")
    input_tsv = sys.argv[1]
    output_csv = sys.argv[2]
    parse_and_tile_blast_results(input_tsv, output_csv)
