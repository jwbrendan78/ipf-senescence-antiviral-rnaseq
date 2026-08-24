# Rubedo RNA-seq workflow

RNA-seq quantification for the Rubedo validation dataset was performed with nf-core/rnaseq using Salmon pseudoalignment.

Key settings:
- nf-core/rnaseq v3.14.0
- Nextflow v25.10.2
- Docker profile
- Salmon pseudoalignment
- No trimming
- No genome alignment
- GENCODE v49 primary assembly FASTA/GTF
- 62 samples analyzed
- Sample `L001003b` excluded because of a known corrupted FASTQ file

## Local setup

Before running the workflow, update `PROJECT_DIR` and `REF_DIR` in
`run_nfcore_rnaseq_salmon_SS_0_240_exclude_corrupted.sh`, or provide them as environment variables.

The repository samplesheet uses relative FASTQ paths under `fastqs/`. FASTQ files should therefore be organized under the configured Rubedo project directory using the directory structure shown in the samplesheet.

GENCODE v49 primary assembly FASTA and GTF files must also be downloaded locally and provided through `REF_DIR`.
