#!/usr/bin/env bash
set -euo pipefail

cd /mnt/d/Rubedo_E001

nextflow run nf-core/rnaseq -r 3.14.0 -profile docker \
  --input /mnt/d/Rubedo_E001/samplesheets/samplesheet_local_SS_0_240_exclude_corrupted.csv \
  --outdir /mnt/d/Rubedo_E001/results/rnaseq_salmon_fastqc_notrim_SS_0_240_exclude_corrupted \
  --pseudo_aligner salmon \
  --skip_alignment \
  --skip_trimming \
  --fasta /mnt/c/Users/jwbre/Documents/Campisi/refs/hg38_gencode_v49/GRCh38.primary_assembly.genome.fa.gz \
  --gtf /mnt/c/Users/jwbre/Documents/Campisi/refs/hg38_gencode_v49/gencode.v49.primary_assembly.annotation.gtf.gz \
  --max_cpus 8 \
  --max_memory 20.GB \
  -work-dir /mnt/d/Rubedo_E001/work \
  -resume
