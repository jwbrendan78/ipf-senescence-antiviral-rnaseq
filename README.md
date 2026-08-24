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

For the current study, a subset of this dataset was reanalyzed using a donor-aware differential expression framework and gene set enrichment analysis.

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
```

## Analysis overview

### Jena primary RNA-seq dataset

The primary Jena RNA-seq dataset was analyzed using:

- Salmon-based RNA-seq quantification
- DESeq2 differential expression analysis with disease, infection, senescence condition, and batch included in the model
- Principal component analysis of variance-stabilized expression values
- Differentially expressed gene overlap analyses
- Gene ontology and pathway enrichment analyses
- Batch-corrected WGCNA co-expression network analysis
- Module enrichment and hub-gene network visualization
- Transcription factor activity inference using CollecTRI and decoupleR

### Rubedo validation dataset

The independent Rubedo dataset was reanalyzed using a donor-aware mixed-model framework with `dream`, followed by preranked gene set enrichment analysis focused on antiviral and interferon-related pathways.

The analysis uses serum-starved samples at timepoints 0 and 240. Sample `L001003b` was excluded because of a known corrupted FASTQ file, resulting in 62 samples in the RNA-seq analysis.

## Software and dependencies

### R analyses

The R scripts use Bioconductor and CRAN packages including:

- DESeq2
- tximport
- limma
- variancePartition
- fgsea
- clusterProfiler
- decoupleR
- AnnotationDbi
- org.Hs.eg.db
- ggplot2
- WGCNA

Additional packages used by individual analyses are loaded within the corresponding scripts.

### Python / WGCNA

WGCNA module construction and module-trait analyses were performed using PyWGCNA in Python.

Python package versions used for the analysis are provided in:

`requirements.txt`

### RNA-seq processing

For the Rubedo validation dataset, RNA-seq processing was performed using:

- nf-core/rnaseq v3.14.0
- Nextflow v25.10.2
- Docker
- Salmon pseudoalignment
- No read trimming
- No genome alignment
- GENCODE v49 primary assembly FASTA and GTF reference files

The exact Rubedo workflow and samplesheet are provided in `workflow/rubedo/`.

Details of the Jena RNA-seq processing workflow are provided in `workflow/jena/README.md`.

## Running the analyses

### Jena dataset

The Jena analysis scripts are intended to be run in the following order:

1. `scripts/jena/01_jena_rnaseq_deseq2_and_downstream.R`  
   Performs differential expression analysis, PCA, DEG overlap analyses, pathway enrichment, transcription factor activity inference, and related downstream analyses.

2. `scripts/jena/02_jena_wgcna_preprocessing.R`  
   Generates the batch-corrected, variance-filtered expression matrix and corresponding metadata used for WGCNA.

3. `scripts/jena/03_jena_wgcna_batch_corrected.ipynb`  
   Performs PyWGCNA module construction and module-trait relationship analyses.

4. `scripts/jena/04_jena_wgcna_networks.R`  
   Generates downstream module dendrograms, hub-gene networks, and related WGCNA visualizations.

### Rubedo validation dataset

The Rubedo validation analysis is performed with:

`scripts/rubedo/01_rubedo_antiviral_gsea_dream.R`

RNA-seq quantification details, the nf-core/rnaseq command, and the analysis samplesheet are provided in:

`workflow/rubedo/`

## Reproducibility and repository preparation

The analysis scripts in this repository reflect the code used for the reported analyses. Prior to public release, the scripts were minimally cleaned with assistance from ChatGPT to improve readability, documentation, and portability across computing environments.

These repository-preparation changes included removal of machine-specific file paths, consolidation of user-configurable paths, clarification of comments, and documentation of required inputs. The underlying statistical models, analysis parameters, contrasts, thresholds, and computational methods were not intentionally altered during this cleanup.

Because the original analyses were performed across specific local R, Python, Nextflow, and filesystem environments, users reproducing the analyses may need to modify project paths, reference-file locations, software environments, or other system-specific settings for their own computing setup.

## License

This repository is available under the MIT License. See `LICENSE` for details.
