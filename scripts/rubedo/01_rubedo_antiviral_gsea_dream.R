# ============================================================
# Part 2: New-dataset antiviral GSEA analysis only
# ============================================================
#
# Public dataset: GEO GSE328392 (SuperSeries GSE328437)
# Original dataset publication: Hughes et al. (2026), npj Aging,
# DOI: 10.1038/s41514-026-00388-4
#
# Purpose
# -------
# Analyze the larger Rubedo RNA-seq dataset at baseline and
# generate antiviral/innate-immune GSEA outputs that can be
# compared directly with the original Jena revision analysis.
#
# Primary comparison
# ------------------
# IPF versus non-PF fibroblasts under:
#   serum = SS
#   timepoint = 0
#
# Replicate handling
# ------------------
# Wells ending in a/b are repeated measurements from the same donor
# and were generated across distinct experimental batches. The primary
# analysis therefore retains all wells, adjusts for batch as a fixed
# effect, and models donor as a random intercept using dream.
#
# For preranked GSEA, genes are ranked by dream's signed standardized
# statistic (z.std), which is comparable across genes despite
# gene-specific degrees of freedom in the mixed model.
#
# Outputs
# -------
# 1. Batch-adjusted, donor-aware dream baseline result
# 2. Complete Hallmark, GO:BP, and Reactome preranked GSEA
# 3. Antiviral/innate-immune candidate tables
# 4. Antiviral wheel plot matching the original plot structure
# 5. Optional old-versus-new baseline pathway comparison table
#
# Run this script from top to bottom in a fresh R session.
# ============================================================

# ============================================================
# 0) USER CONFIGURATION
# ============================================================

# Set this to the local directory containing the Rubedo validation dataset.
project_dir <- "PATH/TO/RUBEDO_PROJECT_DIRECTORY"

pipeline_dir <- file.path(
  project_dir,
  "nextflow"
)

metadata_file <- file.path(
  project_dir,
  "nextflow",
  "samplesheets",
  "Rubedo_metadata.xlsx"
)

metadata_sheet <- "Rubedo_metadata_combined"

output_root <- file.path(
  project_dir,
  "analysis",
  "Part2_new_dataset_GSEA_antiviral"
)

# Baseline subset used for the validation comparison
serum_value <- "SS"
timepoint_value <- 0

# Expected composition based on the experimental design.
# These values are checked and reported but do not force failure.
expected_baseline_wells <- 32
expected_baseline_donors <- 16
expected_nonpf_donors <- 7
expected_ipf_donors <- 9

# Reuse the exact MSigDB gene-set objects from the original
# analysis whenever they are available. This maximizes direct
# comparability and avoids MSigDB-version drift.
# Optional: point this to the GSEA_antiviral output directory from the Jena
# analysis to reuse the exact saved MSigDB gene-set objects. If the files are
# absent, the script falls back to the current MSigDB collections via msigdbr.
old_gsea_root <- "PATH/TO/JENA_GSEA_ANTIVIRAL_RESULTS"

old_geneset_dir <- file.path(
  old_gsea_root,
  "gene_sets"
)

old_gsea_results_rds <- file.path(
  old_gsea_root,
  "tables",
  "GSEA_all_collections_all_contrasts.rds"
)

# ============================================================
# 1) PACKAGES
# ============================================================

required_packages <- c(
  "readxl",
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "edgeR",
  "limma",
  "variancePartition",
  "BiocParallel",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "fgsea",
  "msigdbr",
  "ggplot2",
  "svglite",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages before continuing: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(BiocParallel)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(svglite)
  library(scales)
})

# ============================================================
# 2) OUTPUT DIRECTORIES
# ============================================================

deg_dir <- file.path(output_root, "dream_baseline")
gsea_dir <- file.path(output_root, "GSEA_antiviral")
gsea_table_dir <- file.path(gsea_dir, "tables")
gsea_geneset_dir <- file.path(gsea_dir, "gene_sets")
gsea_rank_dir <- file.path(gsea_dir, "ranked_lists")
gsea_plot_input_dir <- file.path(gsea_dir, "plot_inputs")
gsea_wheel_dir <- file.path(gsea_dir, "wheel_plots")
gsea_wheel_png_dir <- file.path(gsea_wheel_dir, "PNG")
gsea_wheel_svg_dir <- file.path(gsea_wheel_dir, "SVG")
comparison_dir <- file.path(gsea_dir, "old_vs_new_comparison")

for (this_dir in c(
  output_root,
  deg_dir,
  gsea_dir,
  gsea_table_dir,
  gsea_geneset_dir,
  gsea_rank_dir,
  gsea_plot_input_dir,
  gsea_wheel_dir,
  gsea_wheel_png_dir,
  gsea_wheel_svg_dir,
  comparison_dir
)) {
  dir.create(this_dir, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================
# 3) GENERAL HELPERS
# ============================================================

normalize_column_names <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  make.unique(x, sep = "_")
}

find_column <- function(df, candidates, required = TRUE) {
  found <- candidates[candidates %in% names(df)]

  if (length(found) == 0) {
    if (required) {
      stop(
        "Could not find a required metadata column. Tried: ",
        paste(candidates, collapse = ", "),
        "\nAvailable columns: ",
        paste(names(df), collapse = ", ")
      )
    }
    return(NULL)
  }

  found[1]
}

strip_ens_version <- function(x) {
  sub("\\.\\d+$", "", as.character(x))
}

standardize_disease <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  x0 <- gsub("[^a-z0-9]+", "", x0)

  dplyr::case_when(
    x0 %in% c(
      "normal", "healthy", "control", "ctrl", "hlf",
      "nonpf", "nonipf", "noipf", "donorcontrol"
    ) ~ "nonPF",

    x0 %in% c(
      "ipf", "pf", "pulmonaryfibrosis",
      "idiopathicpulmonaryfibrosis"
    ) ~ "IPF",

    TRUE ~ NA_character_
  )
}

map_ids_to_symbols <- function(gene_ids, gene_annot) {
  gene_ids2 <- strip_ens_version(gene_ids)

  annot_clean <- gene_annot %>%
    dplyr::mutate(
      gene_id = strip_ens_version(gene_id),
      gene_name = as.character(gene_name)
    ) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)

  mapped <- tibble::tibble(gene_id = gene_ids2) %>%
    dplyr::left_join(
      annot_clean,
      by = "gene_id"
    )

  missing <- is.na(mapped$gene_name) | mapped$gene_name == ""

  if (any(missing)) {
    fallback_keys <- unique(mapped$gene_id[missing])
    fallback_keys <- fallback_keys[grepl("^ENSG", fallback_keys)]

    if (length(fallback_keys) > 0) {
      fallback <- suppressMessages(
        AnnotationDbi::select(
          org.Hs.eg.db,
          keys = fallback_keys,
          columns = "SYMBOL",
          keytype = "ENSEMBL"
        )
      ) %>%
        dplyr::distinct(ENSEMBL, .keep_all = TRUE)

      mapped <- mapped %>%
        dplyr::left_join(
          fallback,
          by = c("gene_id" = "ENSEMBL")
        ) %>%
        dplyr::mutate(
          gene_name = dplyr::if_else(
            is.na(gene_name) | gene_name == "",
            SYMBOL,
            gene_name
          )
        ) %>%
        dplyr::select(gene_id, gene_name)
    }
  }

  mapped %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)
}

# ============================================================
# 4) LOCATE AND READ NF-CORE SALMON GENE COUNTS
# ============================================================

stopifnot(
  dir.exists(project_dir),
  dir.exists(pipeline_dir),
  file.exists(metadata_file)
)

# Use the final Nextflow result and ignore cached copies under work/.
merged_counts_file <- file.path(
  pipeline_dir,
  "results",
  "rnaseq_salmon_fastqc_notrim_SS_0_240_exclude_corrupted",
  "salmon",
  "salmon.merged.gene_counts.tsv"
)

if (!file.exists(merged_counts_file)) {
  stop(
    "Could not find the final salmon.merged.gene_counts.tsv at:\n",
    merged_counts_file
  )
}

cat(
  "\n================ INPUT FILES ================\n",
  "Merged Salmon gene counts:\n",
  normalizePath(merged_counts_file),
  "\n\nMetadata workbook:\n",
  normalizePath(metadata_file),
  "\n",
  sep = ""
)

counts_tbl <- readr::read_tsv(
  merged_counts_file,
  show_col_types = FALSE,
  progress = FALSE
)

names(counts_tbl) <- normalize_column_names(names(counts_tbl))

gene_id_col <- find_column(
  counts_tbl,
  c("gene_id", "geneid", "ensembl_gene_id")
)

gene_name_col <- find_column(
  counts_tbl,
  c("gene_name", "genename", "gene_symbol", "symbol"),
  required = FALSE
)

if (is.null(gene_name_col)) {
  counts_tbl$gene_name_generated <- NA_character_
  gene_name_col <- "gene_name_generated"
}

gene_annot <- counts_tbl %>%
  dplyr::transmute(
    gene_id = as.character(.data[[gene_id_col]]),
    gene_name = as.character(.data[[gene_name_col]])
  ) %>%
  dplyr::distinct(gene_id, .keep_all = TRUE)

# ============================================================
# 5) READ AND STANDARDIZE METADATA
# ============================================================

meta_raw <- readxl::read_excel(
  metadata_file,
  sheet = metadata_sheet
) %>%
  as.data.frame(stringsAsFactors = FALSE)

names(meta_raw) <- normalize_column_names(names(meta_raw))

sample_col <- find_column(
  meta_raw,
  c(
    "tube_folder_label",
    "folder_label",
    "sample",
    "sample_id",
    "sampleid",
    "library",
    "library_id",
    "library_name"
  )
)

disease_col <- find_column(
  meta_raw,
  c(
    "disease",
    "current_dx",
    "current_diagnosis",
    "dx_at_biopsy",
    "diagnosis",
    "disease_group",
    "group"
  )
)

serum_col <- find_column(
  meta_raw,
  c("serum", "serum_condition")
)

timepoint_col <- find_column(
  meta_raw,
  c("timepoint", "time_point", "time")
)

donor_col <- find_column(
  meta_raw,
  c(
    "cell_label",
    "donor",
    "donor_id",
    "donorid",
    "patient",
    "patient_id",
    "subject",
    "subject_id",
    "individual"
  ),
  required = FALSE
)

batch_col <- find_column(
  meta_raw,
  c(
    "batch",
    "sequencing_batch",
    "seq_batch",
    "experiment",
    "experiment_id",
    "dataset",
    "plate"
  ),
  required = FALSE
)

meta <- meta_raw %>%
  dplyr::mutate(
    sample_id = tolower(trimws(as.character(.data[[sample_col]]))),
    disease = standardize_disease(.data[[disease_col]]),
    serum_standard = toupper(trimws(as.character(.data[[serum_col]]))),
    timepoint_standard = suppressWarnings(
      as.numeric(as.character(.data[[timepoint_col]]))
    )
  )

if (!is.null(donor_col)) {
  meta$donor_id <- trimws(as.character(meta[[donor_col]]))
} else {
  warning(
    "No explicit donor column was found. Donor IDs will be derived ",
    "by removing a terminal a/b from the sample ID."
  )
  meta$donor_id <- sub("[abAB]$", "", meta$sample_id)
}

if (!is.null(batch_col)) {
  meta$batch_candidate <- trimws(as.character(meta[[batch_col]]))
} else {
  meta$batch_candidate <- NA_character_
}

meta <- meta %>%
  dplyr::filter(
    !is.na(sample_id),
    sample_id != ""
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

# ============================================================
# 6) MATCH METADATA TO COUNT COLUMNS AND FILTER BASELINE
# ============================================================

non_sample_count_columns <- c(
  gene_id_col,
  gene_name_col,
  "gene_biotype",
  "gene_type",
  "description"
)

count_sample_columns <- setdiff(
  names(counts_tbl),
  intersect(names(counts_tbl), non_sample_count_columns)
)

matched_samples <- intersect(
  meta$sample_id,
  count_sample_columns
)

cat(
  "\n================ SAMPLE MATCHING ================\n",
  "Metadata rows: ", nrow(meta), "\n",
  "Count sample columns: ", length(count_sample_columns), "\n",
  "Matched samples: ", length(matched_samples), "\n",
  sep = ""
)

if (length(matched_samples) == 0) {
  stop(
    "No metadata sample IDs matched the merged-count columns.\n",
    "First metadata IDs: ",
    paste(utils::head(meta$sample_id, 10), collapse = ", "),
    "\nFirst count columns: ",
    paste(utils::head(count_sample_columns, 10), collapse = ", ")
  )
}

meta_matched <- meta %>%
  dplyr::filter(sample_id %in% matched_samples)

meta_baseline <- meta_matched %>%
  dplyr::filter(
    serum_standard == toupper(serum_value),
    timepoint_standard == timepoint_value,
    disease %in% c("nonPF", "IPF"),
    !is.na(donor_id),
    donor_id != ""
  )

cat(
  "\n================ BASELINE SUBSET ================\n",
  "Serum: ", serum_value, "\n",
  "Timepoint: ", timepoint_value, "\n",
  "Wells retained: ", nrow(meta_baseline), "\n",
  "Donors retained: ", dplyr::n_distinct(meta_baseline$donor_id), "\n",
  sep = ""
)

print(
  meta_baseline %>%
    dplyr::count(disease, name = "wells")
)

print(
  meta_baseline %>%
    dplyr::distinct(donor_id, disease) %>%
    dplyr::count(disease, name = "donors")
)

if (nrow(meta_baseline) != expected_baseline_wells) {
  warning(
    "Expected ", expected_baseline_wells,
    " baseline wells but retained ", nrow(meta_baseline), "."
  )
}

if (dplyr::n_distinct(meta_baseline$donor_id) != expected_baseline_donors) {
  warning(
    "Expected ", expected_baseline_donors,
    " baseline donors but retained ",
    dplyr::n_distinct(meta_baseline$donor_id), "."
  )
}

# Each donor must have one disease label.
donor_disease_check <- meta_baseline %>%
  dplyr::group_by(donor_id) %>%
  dplyr::summarise(
    disease_labels = dplyr::n_distinct(disease),
    .groups = "drop"
  )

stopifnot(all(donor_disease_check$disease_labels == 1))

# ============================================================
# 7) BUILD WELL-LEVEL COUNT MATRIX
# ============================================================

baseline_sample_order <- meta_baseline$sample_id

count_matrix_wells <- counts_tbl %>%
  dplyr::select(
    dplyr::all_of(gene_id_col),
    dplyr::all_of(baseline_sample_order)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE)

rownames(count_matrix_wells) <- count_matrix_wells[[gene_id_col]]
count_matrix_wells[[gene_id_col]] <- NULL

count_matrix_wells <- as.matrix(count_matrix_wells)
storage.mode(count_matrix_wells) <- "numeric"

stopifnot(
  identical(colnames(count_matrix_wells), baseline_sample_order),
  !anyNA(count_matrix_wells),
  all(count_matrix_wells >= 0)
)

# ============================================================
# 8) PREPARE WELL-LEVEL METADATA AND DONOR-LEVEL FILTER COUNTS
# ============================================================

# Keep all 32 wells for the primary mixed-model analysis.
# Each donor contributes two repeated wells from different batches.
meta_well <- meta_baseline %>%
  dplyr::transmute(
    sample_id = as.character(sample_id),
    donor_id = factor(as.character(donor_id)),
    disease = factor(
      as.character(disease),
      levels = c("nonPF", "IPF")
    ),
    batch = factor(trimws(as.character(batch_candidate)))
  ) %>%
  dplyr::arrange(match(sample_id, colnames(count_matrix_wells))) %>%
  as.data.frame(stringsAsFactors = FALSE)

rownames(meta_well) <- meta_well$sample_id

stopifnot(
  identical(rownames(meta_well), colnames(count_matrix_wells)),
  !anyNA(meta_well$donor_id),
  !anyNA(meta_well$disease),
  !anyNA(meta_well$batch),
  nlevels(meta_well$disease) == 2,
  nlevels(meta_well$batch) > 1
)

cat("\n================ WELL-LEVEL MODEL DATA ================\n")
cat("Wells:", nrow(meta_well), "\n")
cat("Donors:", dplyr::n_distinct(meta_well$donor_id), "\n")
print(with(meta_well, table(batch, disease)))

# Create donor-summed counts only for expression filtering and reporting.
# The primary model below still uses the individual well-level columns.
donor_order <- unique(as.character(meta_well$donor_id))

count_matrix_donor <- vapply(
  donor_order,
  function(this_donor) {
    donor_samples <- rownames(meta_well)[
      as.character(meta_well$donor_id) == this_donor
    ]

    rowSums(
      count_matrix_wells[, donor_samples, drop = FALSE]
    )
  },
  FUN.VALUE = numeric(nrow(count_matrix_wells))
)

rownames(count_matrix_donor) <- rownames(count_matrix_wells)
colnames(count_matrix_donor) <- donor_order

donor_meta <- meta_well %>%
  tibble::rownames_to_column("sample_id_rownames") %>%
  dplyr::group_by(donor_id) %>%
  dplyr::summarise(
    disease = dplyr::first(as.character(disease)),
    wells = dplyr::n(),
    batches = paste(
      sort(unique(as.character(batch))),
      collapse = " | "
    ),
    .groups = "drop"
  )

readr::write_csv(
  meta_well %>% tibble::rownames_to_column("sample_rowname"),
  file.path(deg_dir, "well_level_metadata_for_dream.csv")
)

readr::write_csv(
  donor_meta,
  file.path(deg_dir, "donor_level_metadata_summary.csv")
)

# Match the original analysis's detectability rule, but define the
# biological units as donors: >=10 donor-summed counts in >=3 donors.
keep_genes <- rowSums(count_matrix_donor >= 10) >= 3

cat(
  "\nGenes before filtering: ", nrow(count_matrix_wells),
  "\nGenes retained (>=10 counts in >=3 donors): ", sum(keep_genes),
  "\n",
  sep = ""
)

# ============================================================
# 9) VERIFY THE DONOR-AWARE MIXED-MODEL DESIGN
# ============================================================

fixed_design <- stats::model.matrix(
  ~ batch + disease,
  data = meta_well
)

cat("\n================ DREAM DESIGN ================\n")
cat("Formula: ~ batch + disease + (1 | donor_id)\n")
cat("Fixed-effect columns:", ncol(fixed_design), "\n")
cat("Fixed-effect rank:", qr(fixed_design)$rank, "\n")
cat(
  "Fixed-effect design full rank:",
  qr(fixed_design)$rank == ncol(fixed_design),
  "\n"
)

stopifnot(
  qr(fixed_design)$rank == ncol(fixed_design),
  all(table(meta_well$donor_id) == 2)
)

dream_formula <- ~ batch + disease + (1 | donor_id)

writeLines(
  c(
    "Primary method: variancePartition::dream",
    paste0("Formula: ", deparse(dream_formula)),
    "Fixed effect of interest: diseaseIPF (IPF versus non-PF)",
    "Random effect: donor_id",
    "GSEA ranking statistic: dream z.std",
    paste0("Baseline wells: ", nrow(meta_well)),
    paste0("Baseline donors: ", dplyr::n_distinct(meta_well$donor_id)),
    paste0("non-PF donors: ", dplyr::n_distinct(meta_well$donor_id[meta_well$disease == "nonPF"])),
    paste0("IPF donors: ", dplyr::n_distinct(meta_well$donor_id[meta_well$disease == "IPF"])),
    paste0("Genes retained: ", sum(keep_genes))
  ),
  file.path(deg_dir, "analysis_design_summary.txt")
)

# ============================================================
# 10) PRIMARY DREAM ANALYSIS: IPF VERSUS non-PF
# ============================================================

# Salmon gene-level estimated counts can be fractional. The
# limma-voom/dream workflow uses the non-negative numeric counts directly.
dge_new <- edgeR::DGEList(
  counts = count_matrix_wells[keep_genes, , drop = FALSE]
)

dge_new <- edgeR::calcNormFactors(
  dge_new,
  method = "TMM"
)

physical_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))

if (is.na(physical_cores) || physical_cores < 2) {
  dream_workers <- 1L
} else {
  dream_workers <- min(4L, as.integer(physical_cores - 1L))
}

dream_param <- if (dream_workers > 1L) {
  BiocParallel::SnowParam(
    workers = dream_workers,
    type = "SOCK",
    progressbar = TRUE
  )
} else {
  BiocParallel::SerialParam(progressbar = TRUE)
}

cat("dream workers:", dream_workers, "\n")

set.seed(20260729L)

vobj_dream <- variancePartition::voomWithDreamWeights(
  dge_new,
  dream_formula,
  meta_well,
  BPPARAM = dream_param
)

fit_dream <- variancePartition::dream(
  vobj_dream,
  dream_formula,
  meta_well,
  BPPARAM = dream_param
)

fit_dream <- variancePartition::eBayes(fit_dream)

cat("\n================ DREAM COEFFICIENTS ================\n")
print(colnames(fit_dream))

dream_coef <- grep(
  "^diseaseIPF$",
  colnames(fit_dream),
  value = TRUE
)

if (length(dream_coef) != 1) {
  stop(
    "Could not identify exactly one diseaseIPF coefficient. Available: ",
    paste(colnames(fit_dream), collapse = ", ")
  )
}

res_new_baseline <- variancePartition::topTable(
  fit_dream,
  coef = dream_coef,
  number = Inf,
  sort.by = "P"
)

stopifnot(
  "z.std" %in% names(res_new_baseline),
  "adj.P.Val" %in% names(res_new_baseline)
)

annot_for_join <- gene_annot %>%
  dplyr::mutate(
    gene_id_key = strip_ens_version(gene_id)
  ) %>%
  dplyr::distinct(gene_id_key, .keep_all = TRUE) %>%
  dplyr::select(gene_id_key, gene_name)

res_new_baseline_tbl <- res_new_baseline %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::mutate(
    gene_id_key = strip_ens_version(gene_id)
  ) %>%
  dplyr::left_join(
    annot_for_join,
    by = "gene_id_key"
  ) %>%
  dplyr::select(
    gene_id,
    gene_name,
    dplyr::everything(),
    -gene_id_key
  )

readr::write_csv(
  res_new_baseline_tbl,
  file.path(
    deg_dir,
    "dream_new_dataset_IPF_vs_nonPF_SS_time0.csv"
  )
)

saveRDS(
  dge_new,
  file.path(deg_dir, "edgeR_DGEList_well_level.rds")
)

saveRDS(
  vobj_dream,
  file.path(deg_dir, "voomWithDreamWeights_object.rds")
)

saveRDS(
  fit_dream,
  file.path(deg_dir, "dream_fit_IPF_vs_nonPF_SS_time0.rds")
)

saveRDS(
  res_new_baseline,
  file.path(deg_dir, "dream_result_IPF_vs_nonPF_SS_time0.rds")
)

new_deg_summary <- tibble::tibble(
  metric = c(
    "Genes tested",
    "Genes with FDR < 0.05",
    "Positive logFC and FDR < 0.05",
    "Negative logFC and FDR < 0.05"
  ),
  n = c(
    nrow(res_new_baseline),
    sum(res_new_baseline$adj.P.Val < 0.05, na.rm = TRUE),
    sum(
      res_new_baseline$adj.P.Val < 0.05 &
        res_new_baseline$logFC > 0,
      na.rm = TRUE
    ),
    sum(
      res_new_baseline$adj.P.Val < 0.05 &
        res_new_baseline$logFC < 0,
      na.rm = TRUE
    )
  )
)

readr::write_csv(
  new_deg_summary,
  file.path(deg_dir, "dream_new_dataset_summary.csv")
)

cat("\n================ NEW DATASET DREAM SUMMARY ================\n")
print(new_deg_summary)

# ============================================================
# 11) CREATE THE DREAM z.std RANKED LIST
# ============================================================

make_dream_gsea_rank <- function(dream_result, gene_annot) {
  rank_df <- dream_result %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::mutate(
      gene_id = strip_ens_version(gene_id),
      stat = as.numeric(z.std)
    ) %>%
    dplyr::filter(is.finite(stat))

  symbol_map <- map_ids_to_symbols(
    rank_df$gene_id,
    gene_annot
  )

  rank_df <- rank_df %>%
    dplyr::left_join(
      symbol_map,
      by = "gene_id"
    ) %>%
    dplyr::filter(
      !is.na(gene_name),
      gene_name != ""
    ) %>%
    dplyr::group_by(gene_name) %>%
    dplyr::slice_max(
      order_by = abs(stat),
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(stat))

  ranks <- rank_df$stat
  names(ranks) <- rank_df$gene_name
  ranks <- sort(ranks, decreasing = TRUE)

  list(
    ranks = ranks,
    table = rank_df %>%
      dplyr::select(
        gene_symbol = gene_name,
        gene_id,
        stat
      )
  )
}

new_rank_result <- make_dream_gsea_rank(
  dream_result = res_new_baseline,
  gene_annot = gene_annot
)

new_gsea_ranks <- list(
  new_IPF_vs_nonPF_SS_time0 = new_rank_result$ranks
)

readr::write_csv(
  new_rank_result$table,
  file.path(
    gsea_rank_dir,
    "ranked_genes_new_IPF_vs_nonPF_SS_time0_dream_zstd.csv"
  )
)

saveRDS(
  new_gsea_ranks,
  file.path(
    gsea_rank_dir,
    "GSEA_new_dataset_ranked_list_dream_zstd.rds"
  )
)

cat(
  "\nRanked genes for GSEA: ",
  length(new_rank_result$ranks),
  "\nMaximum dream z.std: ",
  max(new_rank_result$ranks),
  "\nMinimum dream z.std: ",
  min(new_rank_result$ranks),
  "\n",
  sep = ""
)

# ============================================================
# 12) LOAD THE SAME MSIGDB COLLECTIONS AS THE OLD ANALYSIS
# ============================================================

old_hallmark_rds <- file.path(
  old_geneset_dir,
  "MSigDB_Hallmark_sets.rds"
)

old_gobp_rds <- file.path(
  old_geneset_dir,
  "MSigDB_GO_BP_sets.rds"
)

old_reactome_rds <- file.path(
  old_geneset_dir,
  "MSigDB_Reactome_sets.rds"
)

reuse_old_gene_sets <- all(
  file.exists(c(
    old_hallmark_rds,
    old_gobp_rds,
    old_reactome_rds
  ))
)

if (reuse_old_gene_sets) {
  hallmark_sets <- readRDS(old_hallmark_rds)
  gobp_sets <- readRDS(old_gobp_rds)
  reactome_sets <- readRDS(old_reactome_rds)

  file.copy(
    old_hallmark_rds,
    file.path(gsea_geneset_dir, basename(old_hallmark_rds)),
    overwrite = TRUE
  )

  file.copy(
    old_gobp_rds,
    file.path(gsea_geneset_dir, basename(old_gobp_rds)),
    overwrite = TRUE
  )

  file.copy(
    old_reactome_rds,
    file.path(gsea_geneset_dir, basename(old_reactome_rds)),
    overwrite = TRUE
  )

  writeLines(
    c(
      "Gene sets reused from the original Jena revision analysis.",
      normalizePath(old_geneset_dir)
    ),
    file.path(gsea_geneset_dir, "gene_set_source.txt")
  )

  cat("\nReused the original saved MSigDB gene-set objects.\n")

} else {
  warning(
    "Original saved gene-set RDS files were not found. ",
    "Retrieving the current MSigDB collections with msigdbr."
  )

  get_msigdb_collection <- function(collection, subcollection = NULL) {
    args <- list(
      db_species = "HS",
      species = "Homo sapiens",
      collection = collection
    )

    if (!is.null(subcollection)) {
      args$subcollection <- subcollection
    }

    tbl <- do.call(msigdbr::msigdbr, args)

    tbl %>%
      dplyr::select(gs_name, gene_symbol, dplyr::everything()) %>%
      dplyr::filter(
        !is.na(gs_name),
        gs_name != "",
        !is.na(gene_symbol),
        gene_symbol != ""
      ) %>%
      dplyr::distinct(gs_name, gene_symbol, .keep_all = TRUE)
  }

  hallmark_tbl <- get_msigdb_collection("H")

  gobp_tbl <- get_msigdb_collection(
    "C5",
    "GO:BP"
  )

  reactome_tbl <- get_msigdb_collection(
    "C2",
    "CP:REACTOME"
  )

  hallmark_sets <- split(
    hallmark_tbl$gene_symbol,
    hallmark_tbl$gs_name
  )

  gobp_sets <- split(
    gobp_tbl$gene_symbol,
    gobp_tbl$gs_name
  )

  reactome_sets <- split(
    reactome_tbl$gene_symbol,
    reactome_tbl$gs_name
  )

  hallmark_sets <- lapply(hallmark_sets, unique)
  gobp_sets <- lapply(gobp_sets, unique)
  reactome_sets <- lapply(reactome_sets, unique)

  saveRDS(
    hallmark_sets,
    file.path(gsea_geneset_dir, "MSigDB_Hallmark_sets.rds")
  )

  saveRDS(
    gobp_sets,
    file.path(gsea_geneset_dir, "MSigDB_GO_BP_sets.rds")
  )

  saveRDS(
    reactome_sets,
    file.path(gsea_geneset_dir, "MSigDB_Reactome_sets.rds")
  )

  writeLines(
    c(
      "Gene sets retrieved with msigdbr during this run.",
      paste0("msigdbr version: ", packageVersion("msigdbr"))
    ),
    file.path(gsea_geneset_dir, "gene_set_source.txt")
  )
}

cat(
  "\n================ GENE-SET COUNTS ================\n",
  "Hallmark pathways: ", length(hallmark_sets), "\n",
  "GO:BP pathways: ", length(gobp_sets), "\n",
  "Reactome pathways: ", length(reactome_sets), "\n",
  sep = ""
)

# ============================================================
# 13) RUN FULL PRERANKED GSEA
# ============================================================

gsea_collections <- list(
  Hallmark = hallmark_sets,
  GO_BP = gobp_sets,
  Reactome = reactome_sets
)

gsea_all_results <- list()
gsea_summary_list <- list()

analysis_number <- 0L

for (contrast_name in names(new_gsea_ranks)) {
  ranks <- new_gsea_ranks[[contrast_name]]

  stopifnot(
    is.numeric(ranks),
    !is.null(names(ranks)),
    all(is.finite(ranks))
  )

  for (collection_name in names(gsea_collections)) {
    analysis_number <- analysis_number + 1L

    cat(
      "\nRunning: ",
      contrast_name,
      " | ",
      collection_name,
      "\n",
      sep = ""
    )

    set.seed(20260729L + analysis_number)

    fgsea_result <- fgsea::fgseaMultilevel(
      pathways = gsea_collections[[collection_name]],
      stats = ranks,
      minSize = 15,
      maxSize = 500,
      eps = 0,
      scoreType = "std",
      nproc = 1
    )

    fgsea_result <- as.data.frame(fgsea_result) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        direction = dplyr::case_when(
          NES > 0 ~ "Positive",
          NES < 0 ~ "Negative",
          TRUE ~ "Neutral"
        )
      ) %>%
      dplyr::arrange(padj, pval)

    result_name <- paste(
      contrast_name,
      collection_name,
      sep = "__"
    )

    gsea_all_results[[result_name]] <- fgsea_result

    collection_output_dir <- file.path(
      gsea_table_dir,
      collection_name
    )

    dir.create(
      collection_output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    fgsea_csv <- fgsea_result %>%
      dplyr::mutate(
        leadingEdge = vapply(
          leadingEdge,
          function(x) paste(x, collapse = ";"),
          character(1)
        )
      )

    readr::write_csv(
      fgsea_csv,
      file.path(
        collection_output_dir,
        paste0(
          "GSEA_",
          contrast_name,
          "_",
          collection_name,
          "_all.csv"
        )
      )
    )

    significant_result <- fgsea_csv %>%
      dplyr::filter(
        !is.na(padj),
        padj < 0.05
      )

    readr::write_csv(
      significant_result,
      file.path(
        collection_output_dir,
        paste0(
          "GSEA_",
          contrast_name,
          "_",
          collection_name,
          "_FDR05.csv"
        )
      )
    )

    gsea_summary_list[[result_name]] <- tibble::tibble(
      contrast = contrast_name,
      collection = collection_name,
      tested_pathways = nrow(fgsea_result),
      significant_FDR05 = nrow(significant_result),
      positive_FDR05 = sum(
        significant_result$NES > 0,
        na.rm = TRUE
      ),
      negative_FDR05 = sum(
        significant_result$NES < 0,
        na.rm = TRUE
      )
    )

    cat(
      "Tested: ", nrow(fgsea_result),
      " | FDR < 0.05: ", nrow(significant_result),
      "\n",
      sep = ""
    )
  }
}

saveRDS(
  gsea_all_results,
  file.path(
    gsea_table_dir,
    "GSEA_new_dataset_all_collections.rds"
  )
)

gsea_run_summary <- dplyr::bind_rows(gsea_summary_list)

readr::write_csv(
  gsea_run_summary,
  file.path(
    gsea_table_dir,
    "GSEA_new_dataset_run_summary.csv"
  )
)

cat("\n================ GSEA RUN SUMMARY ================\n")
print(gsea_run_summary, width = Inf)

# ============================================================
# 14) COMBINE RESULTS AND EXTRACT ANTIVIRAL CANDIDATES
# ============================================================

gsea_plot_all <- dplyr::bind_rows(
  lapply(
    names(gsea_all_results),
    function(result_name) {
      name_parts <- strsplit(
        result_name,
        "__",
        fixed = TRUE
      )[[1]]

      tibble::as_tibble(
        gsea_all_results[[result_name]]
      ) %>%
        dplyr::mutate(
          contrast = name_parts[1],
          collection = name_parts[2]
        )
    }
  )
)

antiviral_terms <- c(
  "INTERFERON",
  "ANTIVIR",
  "VIRAL",
  "VIRUS",
  "RIG_I",
  "RIGI",
  "MDA5",
  "IFIH1",
  "MAVS",
  "DDX58",
  "INNATE_IMMUNE",
  "PATTERN_RECOGNITION",
  "TYPE_I_INTERFERON",
  "TYPE_II_INTERFERON",
  "ISG",
  "OAS",
  "PKR",
  "TLR3",
  "TLR7",
  "TLR8",
  "TLR9",
  "STING",
  "TMEM173"
)

antiviral_pattern <- paste(
  antiviral_terms,
  collapse = "|"
)

gsea_antiviral_candidates <- gsea_plot_all %>%
  dplyr::filter(
    grepl(
      antiviral_pattern,
      pathway,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::mutate(
    minus_log10_FDR = -log10(
      pmax(padj, .Machine$double.xmin)
    ),
    enrichment_direction = dplyr::case_when(
      NES > 0 ~ "Positive",
      NES < 0 ~ "Negative",
      TRUE ~ "Neutral"
    ),
    pathway_label = pathway %>%
      gsub(
        "^(HALLMARK_|GOBP_|GO_|REACTOME_)",
        "",
        .
      ) %>%
      gsub("_", " ", .)
  ) %>%
  dplyr::arrange(
    padj,
    dplyr::desc(abs(NES))
  )

gsea_antiviral_candidates_csv <- gsea_antiviral_candidates %>%
  dplyr::mutate(
    leadingEdge = vapply(
      leadingEdge,
      function(x) paste(x, collapse = ";"),
      character(1)
    )
  )

readr::write_csv(
  gsea_antiviral_candidates_csv,
  file.path(
    gsea_plot_input_dir,
    "GSEA_new_dataset_antiviral_pathway_candidates.csv"
  )
)

cat("\n================ TOP ANTIVIRAL CANDIDATES ================\n")

print(
  gsea_antiviral_candidates %>%
    dplyr::select(
      collection,
      pathway,
      NES,
      padj,
      size
    ) %>%
    dplyr::slice_head(n = 30),
  n = 30,
  width = Inf
)

# ============================================================
# 15) ANTIVIRAL WHEEL PLOT MATCHING THE ORIGINAL STRUCTURE
# ============================================================

wheel_pathways <- tibble::tribble(
  ~pathway, ~pathway_label, ~pathway_order,

  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "Interferon alpha response",
  1,

  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "Interferon gamma response",
  2,

  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
  "IFN-alpha/beta signaling",
  3,

  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "IFN-gamma signaling",
  4,

  "REACTOME_PKR_MEDIATED_SIGNALING",
  "PKR-mediated signaling",
  5,

  "GOBP_RESPONSE_TO_TYPE_I_INTERFERON",
  "Response to type I IFN",
  6,

  "GOBP_RESPONSE_TO_TYPE_II_INTERFERON",
  "Response to type II IFN",
  7,

  "GOBP_RESPONSE_TO_VIRUS",
  "Response to virus",
  8,

  "GOBP_DEFENSE_RESPONSE_TO_VIRUS",
  "Defense response to virus",
  9,

  "GOBP_NEGATIVE_REGULATION_OF_VIRAL_PROCESS",
  "Negative regulation of viral process",
  10,

  "GOBP_VIRAL_GENOME_REPLICATION",
  "Viral genome replication",
  11,

  "GOBP_REGULATION_OF_INNATE_IMMUNE_RESPONSE",
  "Regulation of innate immunity",
  12
)

new_contrast_name <- "new_IPF_vs_nonPF_SS_time0"

new_wheel_data <- gsea_plot_all %>%
  dplyr::filter(
    contrast == new_contrast_name
  ) %>%
  dplyr::inner_join(
    wheel_pathways,
    by = "pathway"
  ) %>%
  dplyr::mutate(
    pathway_label_wrapped = stringr::str_wrap(
      pathway_label,
      width = 22
    ),

    pathway_label_wrapped = factor(
      pathway_label_wrapped,
      levels = stringr::str_wrap(
        wheel_pathways$pathway_label,
        width = 22
      )
    ),

    abs_NES = abs(NES),

    direction_symbol = dplyr::if_else(
      NES >= 0,
      "+",
      "\u2212"
    ),

    enrichment_direction = dplyr::if_else(
      NES >= 0,
      "Positive NES",
      "Negative NES"
    ),

    minus_log10_FDR = -log10(
      pmax(padj, .Machine$double.xmin)
    ),

    minus_log10_FDR_plot = pmin(
      minus_log10_FDR,
      50
    )
  ) %>%
  dplyr::arrange(pathway_order)

if (nrow(new_wheel_data) != nrow(wheel_pathways)) {
  missing_wheel_pathways <- setdiff(
    wheel_pathways$pathway,
    new_wheel_data$pathway
  )

  stop(
    "The fixed antiviral wheel panel is incomplete. Missing: ",
    paste(missing_wheel_pathways, collapse = ", ")
  )
}

new_wheel_data_csv <- new_wheel_data %>%
  dplyr::mutate(
    leadingEdge = vapply(
      leadingEdge,
      function(x) paste(x, collapse = ";"),
      character(1)
    )
  )

readr::write_csv(
  new_wheel_data_csv,
  file.path(
    gsea_wheel_dir,
    "GSEA_new_dataset_antiviral_wheel_plot_data.csv"
  )
)

# The original plots used a -log10(FDR) cap of 50.
wheel_fill_max <- 50

# Preserve the original visual scale when possible.
wheel_y_max <- 4

if (max(new_wheel_data$abs_NES, na.rm = TRUE) > wheel_y_max - 0.25) {
  wheel_y_max <- ceiling(
    max(new_wheel_data$abs_NES, na.rm = TRUE) + 0.5
  )

  warning(
    "The new-data NES values exceed the original y-scale. ",
    "The wheel y-limit was increased to ", wheel_y_max, "."
  )
}

new_wheel_plot <- ggplot2::ggplot(
  new_wheel_data,
  ggplot2::aes(
    x = pathway_label_wrapped,
    y = abs_NES,
    fill = minus_log10_FDR_plot
  )
) +
  ggplot2::geom_col(
    width = 0.82,
    color = "grey25",
    linewidth = 0.4
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      y = abs_NES + 0.18,
      label = direction_symbol
    ),
    size = 7,
    fontface = "bold"
  ) +
  ggplot2::coord_polar(
    start = -pi / 12,
    clip = "off"
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, wheel_y_max),
    breaks = seq(0, wheel_y_max, by = 1),
    expand = c(0, 0)
  ) +
  ggplot2::scale_fill_viridis_c(
    option = "C",
    limits = c(0, wheel_fill_max),
    name = expression(-log[10]("FDR")),
    guide = ggplot2::guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(6.5, "cm"),
      barheight = grid::unit(0.55, "cm")
    )
  ) +
  ggplot2::labs(
    title = "New dataset: baseline IPF vs non-PF",
    subtitle = paste0(
      "SS, timepoint 0; batch-adjusted donor-aware dream",
      "\nBar length = |NES|; +/\u2212 indicates enrichment direction",
      "\nColor scale capped at \u2212log10(FDR) = 50"
    ),
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 18) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 25,
      hjust = 0.5,
      margin = ggplot2::margin(b = 8)
    ),

    plot.subtitle = ggplot2::element_text(
      size = 16,
      hjust = 0.5,
      lineheight = 1.08,
      margin = ggplot2::margin(b = 18)
    ),

    axis.text.x = ggplot2::element_text(
      size = 13,
      face = "bold",
      lineheight = 0.92
    ),

    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),

    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),

    legend.position = "bottom",
    legend.direction = "horizontal",

    legend.title = ggplot2::element_text(
      face = "bold",
      size = 17
    ),

    legend.text = ggplot2::element_text(size = 15),
    legend.margin = ggplot2::margin(t = 14),

    plot.margin = ggplot2::margin(
      40,
      85,
      40,
      85
    )
  )

ggplot2::ggsave(
  filename = file.path(
    gsea_wheel_png_dir,
    "GSEA_wheel_new_IPF_vs_nonPF_SS_time0.png"
  ),
  plot = new_wheel_plot,
  width = 11,
  height = 11.5,
  dpi = 300,
  bg = "white"
)

svglite::svglite(
  file = file.path(
    gsea_wheel_svg_dir,
    "GSEA_wheel_new_IPF_vs_nonPF_SS_time0.svg"
  ),
  width = 11,
  height = 11.5,
  bg = "white"
)

print(new_wheel_plot)
grDevices::dev.off()

cat(
  "\nNew-dataset antiviral wheel plot completed.\n",
  "Output directory:\n",
  normalizePath(gsea_wheel_dir),
  "\n",
  sep = ""
)

# ============================================================
# 16) OPTIONAL DIRECT COMPARISON WITH THE OLD BASELINE GSEA
# ============================================================

if (file.exists(old_gsea_results_rds)) {
  old_gsea_results <- readRDS(old_gsea_results_rds)

  old_result_names <- c(
    "IPF_vs_HLF_QUI_mock__Hallmark",
    "IPF_vs_HLF_QUI_mock__GO_BP",
    "IPF_vs_HLF_QUI_mock__Reactome"
  )

  if (all(old_result_names %in% names(old_gsea_results))) {
    old_baseline_gsea <- dplyr::bind_rows(
      lapply(
        old_result_names,
        function(this_name) {
          parts <- strsplit(this_name, "__", fixed = TRUE)[[1]]

          tibble::as_tibble(
            old_gsea_results[[this_name]]
          ) %>%
            dplyr::mutate(collection = parts[2])
        }
      )
    ) %>%
      dplyr::filter(pathway %in% wheel_pathways$pathway) %>%
      dplyr::select(
        collection,
        pathway,
        old_NES = NES,
        old_padj = padj,
        old_size = size
      )

    new_baseline_gsea <- new_wheel_data %>%
      dplyr::select(
        collection,
        pathway,
        pathway_label,
        pathway_order,
        new_NES = NES,
        new_padj = padj,
        new_size = size
      )

    old_vs_new_baseline <- wheel_pathways %>%
      dplyr::left_join(
        old_baseline_gsea,
        by = "pathway"
      ) %>%
      dplyr::left_join(
        new_baseline_gsea,
        by = c(
          "pathway",
          "pathway_label",
          "pathway_order"
        ),
        suffix = c("_old", "_new")
      ) %>%
      dplyr::mutate(
        old_direction = dplyr::case_when(
          old_NES > 0 ~ "Positive",
          old_NES < 0 ~ "Negative",
          TRUE ~ NA_character_
        ),
        new_direction = dplyr::case_when(
          new_NES > 0 ~ "Positive",
          new_NES < 0 ~ "Negative",
          TRUE ~ NA_character_
        ),
        direction_concordant = old_direction == new_direction
      ) %>%
      dplyr::arrange(pathway_order)

    readr::write_csv(
      old_vs_new_baseline,
      file.path(
        comparison_dir,
        "GSEA_antiviral_old_vs_new_baseline_comparison.csv"
      )
    )

    cat(
      "\n================ OLD VS NEW BASELINE ================\n"
    )

    print(
      old_vs_new_baseline %>%
        dplyr::select(
          pathway_label,
          old_NES,
          old_padj,
          new_NES,
          new_padj,
          direction_concordant
        ),
      n = Inf,
      width = Inf
    )
  } else {
    warning(
      "The old GSEA RDS was found, but the expected baseline result names ",
      "were not all present. The old-versus-new table was skipped."
    )
  }
} else {
  warning(
    "Old GSEA result RDS not found. The new analysis completed, ",
    "but the automatic old-versus-new comparison table was skipped."
  )
}

# ============================================================
# 17) SAVE ANALYSIS OBJECTS AND SESSION INFORMATION
# ============================================================

saveRDS(
  list(
    metadata_well_level = meta_well,
    metadata_donor_level = donor_meta,
    count_matrix_well = count_matrix_wells,
    count_matrix_donor_for_filtering = count_matrix_donor,
    dge = dge_new,
    voom_object = vobj_dream,
    dream_fit = fit_dream,
    dream_result = res_new_baseline,
    ranks = new_gsea_ranks,
    gsea_results = gsea_all_results,
    antiviral_candidates = gsea_antiviral_candidates,
    antiviral_wheel_data = new_wheel_data
  ),
  file.path(
    output_root,
    "Part2_new_dataset_GSEA_analysis_objects.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_root, "sessionInfo.txt")
)

cat(
  "\n============================================================\n",
  "PART 2 NEW-DATASET ANTIVIRAL GSEA COMPLETED\n",
  "============================================================\n",
  "Main output directory:\n",
  normalizePath(output_root),
  "\n",
  sep = ""
)
