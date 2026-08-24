# Jena RNA-seq workflow

Raw RNA-seq reads for the primary Jena dataset were processed using the nf-core/rnaseq workflow with Salmon pseudoalignment against the human GRCh38 reference.

Salmon transcript-level quantifications were summarized to gene-level counts in R using `tximport` and the corresponding transcript-to-gene annotation. These gene-level count estimates were used for downstream DESeq2 and WGCNA analyses.

Raw sequencing data (FASTQ files) and processed bulk RNA-seq data are available from GEO under accession **GSE334185**.

The downstream R analysis uses the saved `tximport` object (`txi_rawcounts.rds`) together with aligned sample metadata.
