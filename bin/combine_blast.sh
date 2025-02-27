# Use the first argument as the input directory.
input_dir=${1:?Input directory required}

echo "Listing all files in the current working directory:"
ls -lR .

echo "Input directory structure (from \$input_dir):"
find "$input_dir" -type f

echo "Finding and concatenating TSV files:"

# Find all TSV files directly inside the input directory and concatenate them.
find "\$input_dir" -mindepth 2 -maxdepth 2 -name "*.tsv" -print \
  | tee found_tsv_files.txt \
  | parallel -j1 "cat {}" \
  >> Blast_hits.tsv

echo "List of found TSV files:"
cat found_tsv_files.txt

echo "Contents of Blast_hits.tsv:"
cat Blast_hits.tsv