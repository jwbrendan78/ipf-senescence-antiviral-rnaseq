# Jena RNA-seq analyses

Scripts for differential expression, WGCNA, transcription factor inference, pathway enrichment, and related analyses of the primary Jena IPF/HLF RNA-seq dataset used in this study.

Raw sequencing data (FASTQ files) and processed bulk RNA-seq data have been deposited in the Gene Expression Omnibus (GEO) under Series accession **GSE334185**.

The analysis workflow includes:
- Salmon-based RNA-seq quantification
- DESeq2 differential expression analysis
- WGCNA co-expression network analysis
- Gene ontology and pathway enrichment
- Transcription factor activity inference using CollecTRI and decoupleR
