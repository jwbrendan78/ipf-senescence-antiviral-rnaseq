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
