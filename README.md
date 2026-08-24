# ipf-senescence-antiviral-rnaseq

Code and analysis workflows associated with the study:

**Cellular senescence dysregulates antiviral interferon responses in idiopathic pulmonary fibrosis**

This repository contains the computational analyses for the primary Jena bulk RNA-seq dataset and the independent Rubedo validation dataset, including differential expression, pathway enrichment, WGCNA, transcription factor activity inference, and antiviral gene set enrichment analyses.

## Data availability

### Primary Jena RNA-seq dataset

Raw sequencing data (FASTQ files) and processed bulk RNA-seq data are available from the Gene Expression Omnibus (GEO) under accession **GSE334185**.

### Rubedo validation dataset

The independent validation analysis uses RNA-seq data from GEO Series **GSE328392**, part of SuperSeries **GSE328437**.

This dataset was originally reported in:

Hughes et al. (2026). *Uncovering senescent fibroblast heterogeneity connects DNA damage response to idiopathic pulmonary fibrosis*. npj Aging.  
DOI: **10.1038/s41514-026-00388-4**

## Repository structure

```text
ipf-senescence-antiviral-rnaseq/
├── README.md
├── LICENSE
├── requirements.txt
├── metadata/
│   └── README.md
├── scripts/
│   ├── jena/
│   │   ├── README.md
│   │   ├── 01_jena_rnaseq_deseq2_and_downstream.R
│   │   ├── 02_jena_wgcna_preprocessing.R
│   │   ├── 03_jena_wgcna_batch_corrected.ipynb
│   │   └── 04_jena_wgcna_networks.R
│   └── rubedo/
│       ├── README.md
│       └── 01_rubedo_antiviral_gsea_dream.R
└── workflow/
    ├── jena/
    │   └── README.md
    └── rubedo/
        ├── README.md
        ├── run_nfcore_rnaseq_salmon_SS_0_240_exclude_corrupted.sh
        ├── samplesheet_SS_0_240_exclude_corrupted.csv
        └── exclude_corrupted_fastqs.txt

## Analysis overview

### Jena primary RNA-seq dataset

The primary dataset was analyzed using:

- Salmon-based RNA-seq quantification
- DESeq2 differential expression analysis with disease, infection, senescence condition, and batch included in the model
- PCA of variance-stabilized expression values
- Differentially expressed gene overlap and pathway enrichment analyses
- Batch-corrected WGCNA co-expression network analysis
- Module enrichment and hub-gene network visualization
- Transcription factor activity inference using CollecTRI and decoupleR

### Rubedo validation dataset

The independent Rubedo dataset was reanalyzed using a donor-aware mixed-model framework with `dream`, followed by preranked gene set enrichment analysis focused on antiviral and interferon-related pathways.
