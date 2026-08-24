# Dataset metadata

This directory contains the analysis-ready metadata tables corresponding to the RNA-seq datasets used in this repository. The GEO records remain the authoritative public archive for the sequencing and processed data.

## Primary Jena RNA-seq dataset

`jena_analysis_metadata.csv` contains the 27 samples used in the batch-adjusted Jena RNA-seq analyses after exclusion of `IPF4_Qm3`. Column names are retained from the analysis metadata so that they correspond directly to the R and PyWGCNA workflows.

Raw sequencing data (FASTQ files) and processed bulk RNA-seq data are available from GEO under accession **GSE334185**.

## Rubedo validation dataset

`rubedo_validation_metadata.csv` contains the 62 serum-starved (SS) samples at timepoints 0 and 240 that were included in the Rubedo nf-core/rnaseq quantification workflow. Sample `L001003b` is not included because its FASTQ file was known to be corrupted.

The validation data are from GEO Series **GSE328392**, part of SuperSeries **GSE328437**, originally reported in Hughes et al. (2026), *Uncovering senescent fibroblast heterogeneity connects DNA damage response to idiopathic pulmonary fibrosis*, npj Aging. DOI: **10.1038/s41514-026-00388-4**.

Only analysis-relevant Rubedo fields are included here; internal storage, location, and collection-tracking fields from the laboratory workbook were omitted.
