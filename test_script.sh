#!/bin/bash -ue

echo "Listing all files in the current working directory:"
ls -lR .

echo "Input directory structure (from results):"
find results -type f

echo "Finding and concatenating TSV files:"

# Since both files are in 'dir1' (which is one level deep inside results),
# we use -mindepth 1 and -maxdepth 1 to search inside the specified directory.
find results -mindepth 1 -maxdepth 1 -name "*.tsv" -print \
  | tee found_tsv_files.txt \
  | parallel -j1 "cat {}" \
  >> Blast_hits.tsv

echo "List of found TSV files:"
cat found_tsv_files.txt

echo "Contents of Blast_hits.tsv:"
cat Blast_hits.tsv
