## RemoteBlast_nr_nt

Nextflow pipeline for remote BLAST against NCBI nr/nt databases.

### Setup

- Nextflow (>= 21.04.0)
- Docker
- Python 3.8+ (pandas, biopython)

### Run

```bash
nextflow run blast.nf
nextflow run blast.nf -resume
```

### Pipeline

1. `blast_remote` — BLASTn against `core_nt` and `tsa_nt`
2. `combineBlast` — merges results into `Blast_hits.tsv`
3. `parse_blast_output` — filters and tiles sequences into `its.csv`

### Output

- `Blast_hits.tsv` — combined BLAST results
- `its.csv` — parsed and tiled sequences

### Troubleshooting

- Process hangs — check internet connection
- Missing results — check error logs in `work/` directory
- Docker issues — ensure Docker is running, try `-resume`
