#!/usr/bin/env bash
set -euo pipefail

# Rubedo RNA-seq quantification workflow
#
# Original analysis:
#   nf-core/rnaseq v3.14.0
#   Nextflow v25.10.2
#   Docker profile
#   Salmon pseudoalignment
#   No trimming
#   No genome alignment
#   GENCODE v49 primary assembly FASTA/GTF
#
# Before running, update PROJECT_DIR and REF_DIR below, or set them as
# environment variables with the same names.

PROJECT_DIR="${PROJECT_DIR:-PATH/TO/RUBEDO_PROJECT_DIRECTORY}"
REF_DIR="${REF_DIR:-PATH/TO/GENCODE_V49_REFERENCE_DIRECTORY}"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

SAMPLESHEET="${SAMPLESHEET:-${SCRIPT_DIR}/samplesheet_SS_0_240_exclude_corrupted.csv}"
OUTDIR="${OUTDIR:-${PROJECT_DIR}/results/rnaseq_salmon_fastqc_notrim_SS_0_240_exclude_corrupted}"
FASTA="${FASTA:-${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz}"
GTF="${GTF:-${REF_DIR}/gencode.v49.primary_assembly.annotation.gtf.gz}"
WORK_DIR="${WORK_DIR:-${PROJECT_DIR}/work}"

if [[ ! -d "${PROJECT_DIR}" ]]; then
  echo "ERROR: PROJECT_DIR does not exist: ${PROJECT_DIR}" >&2
  echo "Update PROJECT_DIR at the top of this script or set it as an environment variable." >&2
  exit 1
fi

if [[ ! -f "${SAMPLESHEET}" ]]; then
  echo "ERROR: Samplesheet not found: ${SAMPLESHEET}" >&2
  exit 1
fi

if [[ ! -f "${FASTA}" ]]; then
  echo "ERROR: FASTA not found: ${FASTA}" >&2
  exit 1
fi

if [[ ! -f "${GTF}" ]]; then
  echo "ERROR: GTF not found: ${GTF}" >&2
  exit 1
fi

cd "${PROJECT_DIR}"

nextflow run nf-core/rnaseq -r 3.14.0 -profile docker \
  --input "${SAMPLESHEET}" \
  --outdir "${OUTDIR}" \
  --pseudo_aligner salmon \
  --skip_alignment \
  --skip_trimming \
  --fasta "${FASTA}" \
  --gtf "${GTF}" \
  --max_cpus 8 \
  --max_memory 20.GB \
  -work-dir "${WORK_DIR}" \
  -resume
