
# Remote BLAST Analysis Pipeline

## Overview
This repository contains a **Nextflow** pipeline for running **Remote BLAST (Basic Local Alignment Search Tool)** analysis against the **NCBI nr/nt database**. The pipeline utilizes a **Docker container** for BLAST execution (and/or Conda for parsing) and supports automated parsing and tiling of BLAST output results. Additionally, a YAML configuration file (`my_env.yaml`) is provided to define the Conda environment for the parsing step.

## Repository Structure

```plaintext
.
├── blast.nf                   # Main Nextflow pipeline script
├── my_env.yaml                # Conda environment definition file
├── README.md                  # Project documentation
├── bin
│   └── tile.py                # Python script for parsing and tiling BLAST results
├── data
│   └── test.fasta             # Example input FASTA file
└── results                    # Output directory (populated by the pipeline)



## Prerequisites

- **Nextflow** (`>= 21.04.0`)
- **Docker** (for containerized execution)
- **Python 3.8+** (for results parsing)
- **Conda**  (parse script enabled using conda)

- **Nextflow (version 24.10.4 or later)**
- **Conda**  
  The pipeline uses a Conda environment defined in `my_env.yaml` for the parsing step. Ensure Conda is installed and accessible.

- **BLAST+**  
  The pipeline assumes that BLAST commands (e.g. `blastn`) are available in your system or container.

- **GNU Parallel**  
  Used in the `combineBlast` process for concatenating TSV files.

## Configuration

The pipeline parameters are defined in the Nextflow script (`blast.nf`). Key parameters include:

- **params.input**: Path to the input FASTA file.
- **params.outdir**: Output directory where results will be stored.
- **params.parseScript**: Path to the Python parsing/tiling script (`bin/tile.py`).
- **params.condaProfile**: Path to the Conda initialization script (e.g. `${baseDir}/my_env.yaml` in this case, though not required when using the `conda` directive).

## Features
- Executes **BLASTn** remotely using **NCBI BLAST API**
- Processes **FASTA input files** and retrieves sequence alignments
- Stores results in `.tsv` format and reformats them into `.csv`
- Supports **Docker-based execution** for reproducibility
- Utilizes **Nextflow DSL2** for workflow management


## Installation
```{r, eval=FALSE}
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/
```

Follow installation instructions from [Docker's official website](https://docs.docker.com/get-docker/).

## Usage
### Running the Pipeline
```{r, eval=FALSE}
nextflow run blast.nf 
```

### Expected Output
- `remote_blast_output_core_nt.tsv`: BLAST results for **core_nt database**
- `remote_blast_output_tsa_nt.tsv`: BLAST results for **tsa_nt database**
- `Blast_hits.csv`: combined results with from `core_nt` and `tsa_nt` remote blast

## Configuration
The pipeline uses **Nextflow configuration profiles**. The default **Docker profile** is defined in `nextflow.config`:
```{r, eval=FALSE}
docker.fixOwnership = true
docker {
    enabled = true
    runOptions = '--platform linux/amd64'
}
process {
    withLabel: 'blast_container' {
        container = 'ncbi/blast:2.16.0'
        cpus = 16
        memory = '16GB'
    }
}
```
  
## Workflow
- **blast.nf**  
The main Nextflow script. It defines the following processes:
### 1️⃣ **BLAST Execution (`blast_remote`)**
- Runs **BLASTn** remotely against `core_nt` and `tsa_nt` databases on an input fasta file
- Outputs results in `.tsv` format
- The `core_nt` search writes a header into its TSV file.
- The `tsa_nt` search is run with a header as well.

### 2️⃣ ** Combines tsa_nt and core_nt results  (`combineBlast`)** 
Merges the TSV outputs by extracting the header from the `core_nt` file and concatenating the contents (skipping duplicate headers) from both TSV files into a merged file (`Blast_hits.tsv`).
- if the `combineBlast` process fails `test_script.sh` can be used to combine the `tsv` files
- make sure `test_script.sh` excutable by running `chmod +x test_script.sh` then run `./test_script.sh`


### 3️⃣ ** Parse Result  (`parse_results`)**
- **parse_blast_output**: Uses a Python script (located in `bin/tile.py`) to parse the merged BLAST TSV file, filter out alignments shorter than 200 bp, and generate overlapping tiling windows (200 bp each with 75 bp overlap) from the subject sequence. Each tile is written into a separate column (e.g. `Tile_1`, `Tile_2`, …). The final parsed output is saved as a CSV file (`its.csv`).
- Parses BLAST output and ensures required fields (`query`, `accession`, `title`, `sequence`, etc.)
- Outputs final results in `.csv`

- **my_env.yaml**  
  A Conda environment file that defines the required packages (Python 3.8, pandas, biopython, nextflow) for the parsing step.

- **bin/tile.py**  
  A Python script that combines parsing and tiling. It reads the merged BLAST TSV file, filters for alignments ≥200 bp, and generates overlapping tiles from the `Subject_Sequence` column. Each tile is output in its own column. The final output is saved as a CSV file.

- **data/test.fasta**  
  An example input FASTA file. Replace or update this file with your actual query sequences.

- **results/**  
  The output directory where the merged TSV file (`Blast_hits.tsv`) and parsed CSV file (`its.csv`) will be stored.

## Example Input
A sample **FASTA file** (`test.fasta`):
```{r, eval=FALSE}
>test_1
CGTAGCTAGCTAGCTAGCTAGCTAGCTAATGGCACCCGCCTTTTGGATACTCTGCAGATGCAAGACATGCTCTAAACTCAGTCACCTGGATGCCGAGGGAATCAGAGTTTTTGGCCGCGAATGGAACAGAAACTCATCAACCTACAAACAACACCTCAGGATCCCCCCAGCCAAAGTATCAACTTCACGT
```

## Troubleshooting
- **Process hangs?** Ensure you have an active **internet connection** as BLAST runs remotely.
- **Missing results?** Check `core_nt_error.log` and `tsa_nt_error.log` for potential errors.
- **Docker issues?** Ensure Docker is running and try `nextflow run blast.nf -resume`.

## References
- [Nextflow Official Website](https://www.nextflow.io)
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [NCBI BLAST API](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs)

## License
This project is licensed under the **MIT License**. See the full license text [here](https://opensource.org/licenses/MIT).

