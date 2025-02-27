
# Remote BLAST Analysis Pipeline

## Overview
This repository contains a **Nextflow** pipeline for running **Remote BLAST (Basic Local Alignment Search Tool)** analysis against the **NCBI nr/nt database**. The pipeline utilizes a **Docker container** for environment consistency and supports automated parsing of BLAST output results.

## Features
- Executes **BLASTn** remotely using **NCBI BLAST API**
- Processes **FASTA input files** and retrieves sequence alignments
- Stores results in `.tsv` format and reformats them into `.csv`
- Supports **Docker-based execution** for reproducibility
- Utilizes **Nextflow DSL2** for workflow management

## Requirements
- **Nextflow** (`>= 21.04.0`)
- **Docker** (for containerized execution)
- **Python 3.8+** (for results parsing)

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
### 1️⃣ **BLAST Execution (`blast_remote`)**
- Runs **BLASTn** remotely against `core_nt` and `tsa_nt`
- Outputs results in `.tsv` format

### 2️⃣ ** Combines tsa_nt and core_nt results  (`combineBlast`)**
- Outputs results in `.tsv` format
- if the `combineBlast` process fails `test_script.sh` can be used to combine the `tsv` files
- make sure `test_script.sh` excutable by running `chmod +x test_script.sh` then run `./test_script.sh`


### 3️⃣ ** Parse Result  (`parse_results`)**
- Parses BLAST output and ensures required fields (`query`, `accession`, `title`, `sequence`, etc.)
- Outputs final results in `.csv`

## Example Input
A sample **FASTA file** (`test.fasta`):
```{r, eval=FALSE}
>test_1
CGTAGCTAGCTAGCTAGCTAGCTAGCTAATGGCACCCGCCTTTTGGATACTCTGCAGATGCAAGACATGCTCTAAACTCAGTCACCTGGATGCCGAGGGAATCAGAGTTTTTGGCCGCGAATGGAACAGAAACTCATCAACCTACAAACAACACCTCAGGATCCCCCCAGCCAAAGTATCAACTTCACGT
```

## Troubleshooting
- **Process hangs?** Ensure you have an active **internet connection** as BLAST runs remotely.
- **Missing results?** Check `core_nt_error.log` and `tsa_nt_error.log` for potential errors.
- **Docker issues?** Ensure Docker is running and try `nextflow run main.nf -resume`.

## References
- [Nextflow Official Website](https://www.nextflow.io)
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [NCBI BLAST API](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs)

## License
This project is licensed under the **MIT License**. See the full license text [here](https://opensource.org/licenses/MIT).

