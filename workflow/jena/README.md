# Jena RNA-seq workflow

Raw RNA-seq reads for the primary Jena dataset were processed using the nf-core/rnaseq workflow with Salmon pseudoalignment against the human GRCh38 reference.

Salmon transcript-level quantifications were summarized to gene-level counts in R using `tximport` and the corresponding transcript-to-gene annotation. These gene-level count estimates were used for downstream DESeq2 and WGCNA analyses.

Raw sequencing data (FASTQ files) and processed bulk RNA-seq data are available from GEO under accession **GSE334185**.

## Downstream analysis inputs

The downstream Jena analyses use:

- `Jena_collab_exps_results/txi_rawcounts.rds` — the gene-level `tximport` object generated from Salmon quantifications
- `Jena_collab_exps_results/salmon.merged.gene_counts.tsv` — gene annotation/count information used by downstream analyses
- `metadata/jena_analysis_metadata.csv` — the 27-sample analysis metadata distributed with this repository

The analysis scripts reorder the repository metadata to match the sample order in the `tximport` count matrix and stop if the samples do not match.
