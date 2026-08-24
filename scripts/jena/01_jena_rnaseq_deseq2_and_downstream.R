# ----------------------------
# Setup
# ----------------------------

# Set this to the local directory containing the Jena RNA-seq analysis files.
# Processed sequencing data are available from GEO: GSE334185.
base_dir <- "PATH/TO/JENA_PROJECT_DIRECTORY"

if (!dir.exists(base_dir)) {
  stop(
    "Project directory not found. Update 'base_dir' at the top of this script ",
    "to the local directory containing the Jena analysis files."
  )
}

setwd(base_dir)

# ============================================================
# Revision analysis output locations
# ============================================================

revision_dir <- file.path(
  base_dir,
  "Jena_collab_revisions"
)

revision_root_dir <- file.path(
  revision_dir,
  "results_revision"
)

# Reproduced outputs from the original script
revision_results_dir <- file.path(
  revision_root_dir,
  "Original_analysis_reproduced"
)

# New analyses added for the revision
reviewer_results_dir <- file.path(
  revision_root_dir,
  "Reviewer_additions"
)

dir.create(
  revision_results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  reviewer_results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat(
  "Original-analysis outputs:\n",
  normalizePath(revision_results_dir),
  "\n\n"
)

cat(
  "Reviewer-addition outputs:\n",
  normalizePath(reviewer_results_dir),
  "\n"
)



# .libPaths(c("C:/Users/jwbre/Rlibs_clean", .libPaths()))
library(tximport)
library(readr)
library(DESeq2)
library(dplyr)


# one time transcript to gene integer fixing method -----------------------
#one time making counts file from transcript abundances if the counts files have integers
# drop_sample <- "IPF4_Qm3"
# 
# # meta already created earlier; just be sure it's trimmed + factorized
# meta <- meta %>% filter(tube_folder_label != drop_sample)
# rownames(meta) <- meta$tube_folder_label

# salmon_dir <- "E:/Jena_collab_exps/results/rnaseq_salmon_fastqc_notrim/salmon"
# 
# files <- file.path(salmon_dir, meta$tube_folder_label, "quant.sf")
# names(files) <- meta$tube_folder_label
# 
# # verify paths
# stopifnot(all(file.exists(files)))
# 
# tx2gene <- read_tsv("Jena_collab_exps_results/tx2gene.tsv",
#                     col_names = c("TXNAME","GENEID"),
#                     show_col_types = FALSE)
# 
# txi <- tximport(files, type="salmon", tx2gene=tx2gene, countsFromAbundance="no")
# 
# # SAVE so you don't need SSD later
# saveRDS(txi, file = "Jena_collab_exps_results/txi_rawcounts.rds")
# saveRDS(meta, file = "Jena_collab_exps_results/meta_aligned.rds")


# important imports -------------------------------------------------------

txi  <- readRDS("Jena_collab_exps_results/txi_rawcounts.rds")
meta <- readRDS("Jena_collab_exps_results/meta_aligned.rds")



# ============================================================
# Revision checkpoint 1:
# Verify inputs, experimental groups, and original DESeq2 model
# ============================================================

# library(DESeq2)
# 
# # Reload the unchanged original analysis inputs
# txi  <- readRDS("Jena_collab_exps_results/txi_rawcounts.rds")
# meta <- readRDS("Jena_collab_exps_results/meta_aligned.rds")

# ------------------------------------------------------------
# 1) Basic input checks
# ------------------------------------------------------------

stopifnot(
  "counts" %in% names(txi),
  all(c("batch", "disease", "infection", "condition") %in% colnames(meta)),
  identical(rownames(meta), colnames(txi$counts)),
  anyDuplicated(rownames(meta)) == 0,
  anyDuplicated(colnames(txi$counts)) == 0
)

design_vars <- c("batch", "disease", "infection", "condition")

stopifnot(
  all(complete.cases(meta[, design_vars]))
)

cat("\n================ INPUT CHECKS ================\n")
cat("Genes in tximport object:", nrow(txi$counts), "\n")
cat("Samples in tximport object:", ncol(txi$counts), "\n")
cat("Samples in metadata:", nrow(meta), "\n")
cat("Metadata and count-column order identical: ",
    identical(rownames(meta), colnames(txi$counts)), "\n",
    sep = "")

# ------------------------------------------------------------
# 2) Re-create the original reference levels
# ------------------------------------------------------------

meta$batch     <- factor(meta$batch)
meta$condition <- relevel(factor(meta$condition), ref = "CTL")
meta$disease   <- relevel(factor(meta$disease), ref = "normal")
meta$infection <- relevel(factor(meta$infection), ref = "mock")

cat("\n================ REFERENCE LEVELS ================\n")
cat("Batch levels:", paste(levels(meta$batch), collapse = ", "), "\n")
cat("Disease levels:", paste(levels(meta$disease), collapse = ", "), "\n")
cat("Infection levels:", paste(levels(meta$infection), collapse = ", "), "\n")
cat("Condition levels:", paste(levels(meta$condition), collapse = ", "), "\n")

# ------------------------------------------------------------
# 3) Show group and batch composition
# ------------------------------------------------------------

cat("\n================ BIOLOGICAL GROUP COUNTS ================\n")
print(with(meta, table(disease, infection, condition)))

cat("\n================ DISEASE BY BATCH ================\n")
print(with(meta, table(disease, batch)))

cat("\n================ CONDITION BY BATCH ================\n")
print(with(meta, table(condition, batch)))

if ("cell_label" %in% colnames(meta)) {
  cat("\n================ SAMPLES PER DONOR ================\n")
  print(sort(table(meta$cell_label)))
}

# ------------------------------------------------------------
# 4) Verify that the design matrix is full rank
# ------------------------------------------------------------

design_formula <- ~ batch + disease * infection * condition

design_matrix <- model.matrix(
  design_formula,
  data = as.data.frame(meta)
)

design_rank <- qr(design_matrix)$rank
design_columns <- ncol(design_matrix)

cat("\n================ DESIGN MATRIX ================\n")
cat("Model:", deparse(design_formula), "\n")
cat("Number of design columns:", design_columns, "\n")
cat("Design-matrix rank:", design_rank, "\n")
cat("Full rank:", design_rank == design_columns, "\n")

stopifnot(design_rank == design_columns)
# The DESeq2 model is fitted once later in the original
# "filtering and deseq" section.

cat("\nInput and design verification completed successfully.\n")



# ============================================================
# Revision checkpoint 2:
# Inspect original QUI/mock baseline samples
# ============================================================

baseline_meta <- meta[
  meta$condition == "CTL" &
    meta$infection == "mock",
  ,
  drop = FALSE
]

baseline_meta <- droplevels(baseline_meta)

cat("\n================ BASELINE SAMPLE COUNT ================\n")
cat("Total QUI/mock samples:", nrow(baseline_meta), "\n")

cat("\n================ BASELINE DISEASE COUNTS ================\n")
print(table(baseline_meta$disease))

cat("\n================ BASELINE DISEASE BY BATCH ================\n")
print(with(baseline_meta, table(disease, batch)))

cat("\n================ BASELINE SAMPLE DETAILS ================\n")

baseline_details <- data.frame(
  sample = rownames(baseline_meta),
  donor = if ("cell_label" %in% colnames(baseline_meta)) {
    as.character(baseline_meta$cell_label)
  } else {
    NA_character_
  },
  disease = as.character(baseline_meta$disease),
  batch = as.character(baseline_meta$batch),
  stringsAsFactors = FALSE
)

print(baseline_details, row.names = FALSE)



# ============================================================
# GEO processed file: tximport gene-level estimated counts
# ============================================================

# gene_annot <- readr::read_tsv(
#   "Jena_collab_exps_results/salmon.merged.gene_counts.tsv",
#   show_col_types = FALSE
# ) %>%
#   dplyr::select(gene_id, gene_name) %>%
#   dplyr::distinct()
# 
# counts_out <- txi$counts %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column("gene_id") %>%
#   dplyr::left_join(gene_annot, by = "gene_id") %>%
#   dplyr::relocate(gene_name, .after = gene_id)
# 
# readr::write_tsv(
#   counts_out,
#   file.path(geo_root, "Jena_collab_tximport_gene_estimated_counts.tsv")
# )

# ============================================================
# GEO helper output: raw files + key characteristics
# ============================================================

# geo_root <- file.path(revision_results_dir, "GEO_output")
# dir.create(geo_root, recursive = TRUE, showWarnings = FALSE)
# 
# # Make sure metadata is aligned to tximport counts
# stopifnot(identical(rownames(meta), colnames(txi$counts)))
# 
# meta_all <- meta %>%
#   mutate(
#     batch     = as.character(batch),
#     condition = as.character(condition),
#     infection = as.character(infection),
#     disease   = as.character(disease),
#     donor     = as.character(cell_label)
#   )
# 
# geo_raw_files <- meta_all %>%
#   mutate(
#     condition_geo = dplyr::recode(condition, "CTL" = "QUI", "IR" = "SEN"),
#     disease_geo   = dplyr::recode(disease, "normal" = "healthy"),
#     serum_clean   = ifelse(is.na(serum) | serum == "", "not provided", as.character(serum)),
#     sex_clean     = ifelse(is.na(sex) | sex == "", "not provided", as.character(sex)),
#     age_clean     = ifelse(is.na(age) | age == "", "not provided", as.character(age))
#   ) %>%
#   rowwise() %>%
#   transmute(
#     title = paste(
#       c(
#         disease_geo,
#         infection,
#         condition_geo,
#         paste("donor", donor),
#         batch
#       )[c(
#         disease_geo,
#         infection,
#         condition_geo,
#         paste("donor", donor),
#         batch
#       ) != ""],
#       collapse = ", "
#     ),
#     
#     description = paste0(
#       "Donor ID: ", donor,
#       "; disease: ", disease_geo,
#       "; infection: ", infection,
#       "; senescence treatment: ", condition_geo,
#       "; batch: ", batch
#     ),
#     
#     `*raw file` = paste0(tube_folder_label, "_1.fq.gz"),
#     `raw file`  = paste0(tube_folder_label, "_2.fq.gz"),
#     
#     `characteristics: donor` = donor,
#     `characteristics: disease` = disease_geo,
#     `characteristics: infection` = infection,
#     `characteristics: senescence treatment` = condition_geo,
#     `characteristics: batch` = batch,
#     
#     # Optional: keep only if you want these in GEO
#     `characteristics: serum` = serum_clean,
#     `characteristics: sex` = sex_clean,
#     `characteristics: age` = age_clean
#   ) %>%
#   ungroup()
# 
# out_geo_raw_path <- file.path(geo_root, "geo_output_raw_files.tsv")
# readr::write_tsv(geo_raw_files, out_geo_raw_path)
# 
# cat("Wrote GEO raw file helper TSV:\n", out_geo_raw_path, "\n")

# Batch-aware script checks ------------------------------------------------
# This script assumes meta_aligned.rds has already been rebuilt from the
# updated metadata CSV and contains the batch column.
stopifnot("batch" %in% colnames(meta))
stopifnot(identical(rownames(meta), colnames(txi$counts)))
meta$batch <- factor(meta$batch)
cat("Batch counts in meta_aligned.rds:\n")
print(table(meta$batch, useNA = "ifany"))


# PCAs --------------------------------------------------------------------


pca_dir <- file.path(revision_results_dir, "PCA")
pca_png_dir <- file.path(pca_dir, "PNG")
pca_svg_dir <- file.path(pca_dir, "SVG")
dir.create(pca_png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pca_svg_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(dplyr)
library(svglite)

# helper: save both PNG and SVG with the same stem
save_plot_dual <- function(plot, out_stem, png_dir, svg_dir,
                           width = 9, height = 8, dpi = 300) {
  ggsave(
    filename = file.path(png_dir, paste0(out_stem, ".png")),
    plot = plot,
    width = width, height = height, dpi = dpi
  )
  
  svglite::svglite(
    file = file.path(svg_dir, paste0(out_stem, ".svg")),
    width = width, height = height
  )
  print(plot)
  dev.off()
}

dds_pca <- DESeqDataSetFromTximport(txi, colData=meta, design=~1)
vsd <- vst(dds_pca, blind=TRUE)

pca_df <- plotPCA(
  vsd,
  intgroup = c("infection", "disease", "condition", "batch"),
  returnData = TRUE
) %>%
  mutate(
    disease = recode(disease, "normal" = "healthy"),
    condition = recode(condition, "CTL" = "QUI", "IR" = "SEN"),
    batch = factor(batch)
  )

percentVar <- round(100 * attr(pca_df, "percentVar"))

p_all <- ggplot(
  pca_df,
  aes(PC1, PC2,
      color = infection,
      shape = condition,
      label = disease)
) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(
    size = 3,
    max.overlaps = 50,
    box.padding = 0.4,
    point.padding = 0.3,
    show.legend = FALSE
  ) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(
    "PCA of RNA-seq Samples",
    subtitle = "Variance-stabilized counts (tximport → DESeq2)"
  ) +
  theme_classic()

p_all
save_plot_dual(
  plot = p_all,
  out_stem = "PCA_infection_disease_condition",
  png_dir = pca_png_dir,
  svg_dir = pca_svg_dir,
  width = 9, height = 8, dpi = 300
)

p2 <- ggplot(
  pca_df,
  aes(PC1, PC2,
      color = condition,
      shape = disease)
) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(
    "PCA of RNA-seq Samples",
    subtitle = "Variance-stabilized counts (tximport → DESeq2)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

p2
save_plot_dual(
  plot = p2,
  out_stem = "PCA_disease_condition",
  png_dir = pca_png_dir,
  svg_dir = pca_svg_dir,
  width = 9, height = 8, dpi = 300
)

# New PCA: batch + condition + infection ----------------------------------

p_batch <- ggplot(
  pca_df,
  aes(PC1, PC2,
      color = batch,
      shape = condition,
      label = infection)
) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(
    size = 3,
    max.overlaps = 50,
    box.padding = 0.4,
    point.padding = 0.3,
    show.legend = FALSE
  ) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(
    "PCA of RNA-seq Samples by Batch",
    subtitle = "Color = batch; shape = condition; label = infection"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

p_batch

save_plot_dual(
  plot = p_batch,
  out_stem = "PCA_batch_condition_infection",
  png_dir = pca_png_dir,
  svg_dir = pca_svg_dir,
  width = 9, height = 8, dpi = 300
)


# Batch-corrected PCA visualization only ----------------------------------

library(limma)

# Preserve the biological design while removing batch for visualization
design_bc <- model.matrix(
  ~ disease * infection * condition,
  data = as.data.frame(colData(vsd))
)

# Remove batch effect from VST assay for PCA visualization only
vst_mat_batch_corrected <- limma::removeBatchEffect(
  assay(vsd),
  batch = colData(vsd)$batch,
  design = design_bc
)

# Run PCA manually on batch-corrected VST matrix
pca_bc <- prcomp(t(vst_mat_batch_corrected))

percentVar_bc <- round(
  100 * (pca_bc$sdev^2 / sum(pca_bc$sdev^2))
)

pca_bc_df <- as.data.frame(pca_bc$x[, 1:2]) %>%
  mutate(
    name = rownames(.)
  ) %>%
  left_join(
    as.data.frame(colData(vsd)) %>%
      mutate(name = rownames(.)),
    by = "name"
  ) %>%
  mutate(
    disease = recode(disease, "normal" = "healthy"),
    condition = recode(condition, "CTL" = "QUI", "IR" = "SEN"),
    batch = factor(batch)
  )


p_bc_infection <- ggplot(
  pca_bc_df,
  aes(PC1, PC2,
      color = infection,
      shape = condition,
      label = disease)
) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(
    size = 3,
    max.overlaps = 50,
    box.padding = 0.4,
    point.padding = 0.3,
    show.legend = FALSE
  ) +
  xlab(paste0("PC1: ", percentVar_bc[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_bc[2], "% variance")) +
  ggtitle(
    "Batch-corrected PCA of RNA-seq Samples",
    subtitle = "Batch removed from VST matrix; color = infection"
  ) +
  theme_classic()

p_bc_infection

save_plot_dual(
  plot = p_bc_infection,
  out_stem = "PCA_batch_corrected_infection_disease_condition",
  png_dir = pca_png_dir,
  svg_dir = pca_svg_dir,
  width = 9, height = 8, dpi = 300
)


p_bc_batch_check <- ggplot(
  pca_bc_df,
  aes(PC1, PC2,
      color = batch,
      shape = condition,
      label = infection)
) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(
    size = 3,
    max.overlaps = 50,
    box.padding = 0.4,
    point.padding = 0.3,
    show.legend = FALSE
  ) +
  xlab(paste0("PC1: ", percentVar_bc[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_bc[2], "% variance")) +
  ggtitle(
    "Batch-corrected PCA of RNA-seq Samples",
    subtitle = "Batch removed from VST matrix; color = batch"
  ) +
  theme_classic()

p_bc_batch_check

save_plot_dual(
  plot = p_bc_batch_check,
  out_stem = "PCA_batch_corrected_batch_condition_infection",
  png_dir = pca_png_dir,
  svg_dir = pca_svg_dir,
  width = 9, height = 8, dpi = 300
)


# filter info section -----------------------------------------------------


#filter step based on low count diagnostics after running deseq
# --- Low-count diagnostics (before DESeq) ---
#   > cts <- counts(dds_3way)
# >
#   > # how many samples have at least X counts per gene?
#   > nsamp_ge10 <- rowSums(cts >= 10)
# > nsamp_ge5  <- rowSums(cts >= 5)
# > nsamp_ge20 <- rowSums(cts >= 20)
# >
#   > summary(nsamp_ge10)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.000   0.000   0.000   5.682   3.000  27.000
# >
#   > # How many genes would be kept under common rules?
#   > c(
#     +     keep_ge10_in_2  = sum(nsamp_ge10 >= 2),
#     +     keep_ge10_in_4  = sum(nsamp_ge10 >= 4),
#     +     keep_ge10_in_6  = sum(nsamp_ge10 >= 6),
#     +     keep_ge10_in_8  = sum(nsamp_ge10 >= 8),
#     +     keep_ge5_in_6   = sum(nsamp_ge5  >= 6),
#     +     keep_ge20_in_6  = sum(nsamp_ge20 >= 6)
#     + )
# keep_ge10_in_2 keep_ge10_in_4 keep_ge10_in_6 keep_ge10_in_8  keep_ge5_in_6 keep_ge20_in_6
# 20839          19144          18227          17543          20855          16041

# Keep genes that are detectably expressed (≥10 counts) in at least one full biological group
  #your smallest group size is 6
# requiring ≥6 samples with ≥10 counts ensures:
#   the gene is consistently expressed, not a single-sample artifact
# the gene could plausibly contribute to within-group comparisons

# with(meta, table(disease, infection, condition))
# , , condition = CTL
#
# infection
# disease  mock IAV
# normal    3   3
# IPF       3   4
#
# , , condition = IR
#
# infection
# disease  mock IAV
# normal    3   3
# IPF       4   4

#smallest group size is actually 3 so change the group size to 3


# filtering and deseq -----------------------------------------------------

#deseq with interaction variable since we are asking if disease plays a role on infection
meta$batch     <- factor(meta$batch)
meta$condition <- relevel(factor(meta$condition), ref = "CTL")
meta$disease   <- relevel(factor(meta$disease), ref = "normal")
meta$infection <- relevel(factor(meta$infection), ref = "mock")

dds_3way <- DESeqDataSetFromTximport(
  txi,
  colData = meta,
  design = ~ batch + disease * infection * condition
)

keep <- rowSums(counts(dds_3way) >= 10) >= 3
dds_3way <- dds_3way[keep, ]

dds_3way <- DESeq(dds_3way)

# Identify any genes whose DESeq2 coefficient fitting did not converge
nonconverged_genes <- rownames(dds_3way)[
  !is.na(mcols(dds_3way)$betaConv) &
    !mcols(dds_3way)$betaConv
]

cat(
  "\nNumber of non-converged genes:",
  length(nonconverged_genes),
  "\n"
)

print(nonconverged_genes)

resultsNames(dds_3way)

nm <- resultsNames(dds_3way)
nm

deg_dir <- file.path(revision_results_dir, "DESeq2", "DEGs_3way_senescence")
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

library(readr)
library(dplyr)

gene_annot <- read_tsv(
  "Jena_collab_exps_results/salmon.merged.gene_counts.tsv",
  show_col_types = FALSE
) %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::distinct()

# helper: figure-friendly gene labels
clean_gene_label <- function(gene_name, gene_id) {
  out <- ifelse(!is.na(gene_name) & gene_name != "", gene_name, NA_character_)
  out[grepl("^ENSG", out)] <- NA_character_
  out
}


# ============================================================
# Revision analysis 1:
# Original-dataset baseline disease comparison
#
# Contrast:
#   IPF QUI mock vs non-PF QUI mock
#
# Because mock and CTL are the reference levels, this is the
# disease_IPF_vs_normal coefficient from the full 3-way model.
#
# Positive log2FC = higher expression in IPF
# Negative log2FC = lower expression in IPF
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ashr)
  library(dplyr)
  library(readr)
  library(tibble)
})

# ------------------------------------------------------------
# 1) Output directory
# ------------------------------------------------------------

baseline_deg_dir <- file.path(
  reviewer_results_dir,
  "DESeq2",
  "Baseline_IPF_vs_nonPF_QUI_mock"
)

dir.create(
  baseline_deg_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2) Confirm coefficient exists
# ------------------------------------------------------------

baseline_coef <- "disease_IPF_vs_normal"

stopifnot(
  baseline_coef %in% resultsNames(dds_3way)
)

cat("\n================ BASELINE CONTRAST ================\n")
cat("Coefficient:", baseline_coef, "\n")
cat("Comparison: IPF QUI mock vs non-PF QUI mock\n")
cat("Positive log2FC means higher expression in IPF.\n")

# ------------------------------------------------------------
# 3) Extract the raw DESeq2 result
# ------------------------------------------------------------

res_baseline_raw <- results(
  dds_3way,
  name = baseline_coef,
  alpha = 0.05
)

# ------------------------------------------------------------
# 4) Shrink log2 fold changes using the same ashr approach
#    already used in the original analysis
# ------------------------------------------------------------

res_baseline_shr <- lfcShrink(
  dds_3way,
  coef = baseline_coef,
  type = "ashr"
)

# ------------------------------------------------------------
# 5) Load gene annotation
# ------------------------------------------------------------

# gene_annot <- readr::read_tsv(
#   "Jena_collab_exps_results/salmon.merged.gene_counts.tsv",
#   show_col_types = FALSE
# ) %>%
#   dplyr::select(gene_id, gene_name) %>%
#   dplyr::distinct(gene_id, .keep_all = TRUE)

# ------------------------------------------------------------
# 6) Convert results to annotated tables
# ------------------------------------------------------------

result_to_table <- function(res_object, gene_annot) {
  
  as.data.frame(res_object) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(
      gene_annot,
      by = "gene_id"
    ) %>%
    dplyr::relocate(
      gene_name,
      .after = gene_id
    ) %>%
    dplyr::arrange(
      is.na(padj),
      padj
    )
}

baseline_raw_tbl <- result_to_table(
  res_baseline_raw,
  gene_annot
)

baseline_shr_tbl <- result_to_table(
  res_baseline_shr,
  gene_annot
)

# ------------------------------------------------------------
# 7) Save tables and DESeq2 result objects
# ------------------------------------------------------------

readr::write_csv(
  baseline_raw_tbl,
  file.path(
    baseline_deg_dir,
    "DEG_IPF_QUI_mock_vs_nonPF_QUI_mock_raw.csv"
  )
)

readr::write_csv(
  baseline_shr_tbl,
  file.path(
    baseline_deg_dir,
    "DEG_IPF_QUI_mock_vs_nonPF_QUI_mock_shrunk_ashr.csv"
  )
)

saveRDS(
  res_baseline_raw,
  file.path(
    baseline_deg_dir,
    "res_IPF_QUI_mock_vs_nonPF_QUI_mock_raw.rds"
  )
)

saveRDS(
  res_baseline_shr,
  file.path(
    baseline_deg_dir,
    "res_IPF_QUI_mock_vs_nonPF_QUI_mock_shrunk_ashr.rds"
  )
)

# ------------------------------------------------------------
# 8) Summarize the result
#
# Statistical significance comes from the DESeq2 padj.
# Effect-size classification uses the shrunken log2FC.
# ------------------------------------------------------------

baseline_sig_tbl <- baseline_shr_tbl %>%
  dplyr::mutate(
    significant = !is.na(padj) & padj < 0.05,
    direction = dplyr::case_when(
      significant & log2FoldChange >= 0.5  ~ "Higher in IPF",
      significant & log2FoldChange <= -0.5 ~ "Lower in IPF",
      significant                         ~ "Significant, |log2FC| < 0.5",
      TRUE                                ~ "Not significant"
    )
  )

baseline_summary <- tibble::tibble(
  metric = c(
    "Genes tested with non-NA adjusted p-value",
    "Genes with padj < 0.05",
    "Higher in IPF: padj < 0.05 and log2FC >= 0.5",
    "Lower in IPF: padj < 0.05 and log2FC <= -0.5"
  ),
  n = c(
    sum(!is.na(baseline_sig_tbl$padj)),
    sum(baseline_sig_tbl$significant),
    sum(baseline_sig_tbl$direction == "Higher in IPF"),
    sum(baseline_sig_tbl$direction == "Lower in IPF")
  )
)

readr::write_csv(
  baseline_summary,
  file.path(
    baseline_deg_dir,
    "baseline_DEG_summary.csv"
  )
)

cat("\n================ BASELINE DEG SUMMARY ================\n")
print(baseline_summary)

# ------------------------------------------------------------
# 9) Print the top 15 genes by adjusted p-value
# ------------------------------------------------------------

baseline_top15 <- baseline_shr_tbl %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::select(
    gene_id,
    gene_name,
    baseMean,
    log2FoldChange,
    lfcSE,
    pvalue,
    padj
  ) %>%
  dplyr::slice_head(n = 15)

cat("\n================ TOP 15 BASELINE GENES ================\n")
print(
  as.data.frame(baseline_top15),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10) Check the previously identified non-converged gene
# ------------------------------------------------------------

cat("\n================ NON-CONVERGED GENE RESULT ================\n")

print(
  baseline_raw_tbl %>%
    dplyr::filter(gene_id %in% nonconverged_genes)
)

cat(
  "\nBaseline result files written to:\n",
  normalizePath(baseline_deg_dir),
  "\n"
)



###############################################################################
## DESeq2 3-way analysis output
## Model: disease * infection * condition
##
## Reference levels (baseline):
##   disease   = normal
##   infection = mock
##   condition = CTL
##
## Naming convention used below:
##   "<context>: A vs B"
##   where context explicitly states disease + infection (or condition)
###############################################################################

library(DESeq2)
library(dplyr)
library(readr)

## ---------------------------------------------------------------------------
## 1) Output directories
## ---------------------------------------------------------------------------

deg_dir <- file.path(revision_results_dir, "DESeq2", "DEGs_3way_clear_labels")
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 2) Gene annotation (Ensembl → gene symbol)
## ---------------------------------------------------------------------------

if (!exists("gene_annot")) {
  gene_annot <- read_tsv(
    "Jena_collab_exps_results/salmon.merged.gene_counts.tsv",
    show_col_types = FALSE
  ) %>%
    select(gene_id, gene_name) %>%
    distinct()
}

## ---------------------------------------------------------------------------
## 3) Helper function: write DESeq2 results to CSV
## ---------------------------------------------------------------------------

write_res_csv <- function(res, filename, gene_annot = NULL) {
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  
  if (!is.null(gene_annot)) {
    df <- left_join(df, gene_annot, by = "gene_id")
  }
  
  col_order <- c("gene_id", "gene_name",
                 setdiff(colnames(df), c("gene_id", "gene_name")))
  df <- df[, col_order[col_order %in% colnames(df)]]
  
  df <- df %>% arrange(is.na(padj), padj)
  write.csv(df, filename, row.names = FALSE)
}


###############################################################################
## SECTION A: Senescence effect (IR vs CTL) within each disease + infection
###############################################################################

# normal, mock: IR vs CTL
res_normal_mock_IR_vs_CTL <- results(
  dds_3way,
  name = "condition_IR_vs_CTL"
)

# normal, IAV: IR vs CTL
res_normal_IAV_IR_vs_CTL <- results(
  dds_3way,
  contrast = list(c("condition_IR_vs_CTL", "infectionIAV.conditionIR"))
)

# IPF, mock: IR vs CTL
res_IPF_mock_IR_vs_CTL <- results(
  dds_3way,
  contrast = list(c("condition_IR_vs_CTL", "diseaseIPF.conditionIR"))
)

# IPF, IAV: IR vs CTL
res_IPF_IAV_IR_vs_CTL <- results(
  dds_3way,
  contrast = list(c("condition_IR_vs_CTL",
                    "diseaseIPF.conditionIR",
                    "infectionIAV.conditionIR",
                    "diseaseIPF.infectionIAV.conditionIR"))
)

###############################################################################
## SECTION B: Infection effect (IAV vs mock) within each disease + condition
###############################################################################

# normal CTL: IAV vs mock
res_normal_CTL_IAV_vs_mock <- results(
  dds_3way,
  name = "infection_IAV_vs_mock"
)

# normal IR: IAV vs mock
res_normal_IR_IAV_vs_mock <- results(
  dds_3way,
  contrast = list(c("infection_IAV_vs_mock", "infectionIAV.conditionIR"))
)

# IPF CTL: IAV vs mock
res_IPF_CTL_IAV_vs_mock <- results(
  dds_3way,
  contrast = list(c("infection_IAV_vs_mock", "diseaseIPF.infectionIAV"))
)

# IPF IR: IAV vs mock
res_IPF_IR_IAV_vs_mock <- results(
  dds_3way,
  contrast = list(c("infection_IAV_vs_mock",
                    "diseaseIPF.infectionIAV",
                    "infectionIAV.conditionIR",
                    "diseaseIPF.infectionIAV.conditionIR"))
)

###############################################################################
## SECTION C: Interaction / comparison-of-comparisons
###############################################################################

# Does infection change senescence effect in normal?
res_normal_infection_modifies_senescence <- results(
  dds_3way,
  name = "infectionIAV.conditionIR"
)

# Does infection change senescence effect in IPF?
res_IPF_infection_modifies_senescence <- results(
  dds_3way,
  contrast = list(c("infectionIAV.conditionIR",
                    "diseaseIPF.infectionIAV.conditionIR"))
)

# Does disease change how infection modifies senescence? (3-way)
res_3way_disease_modifies_infection_senescence <- results(
  dds_3way,
  name = "diseaseIPF.infectionIAV.conditionIR"
)

###############################################################################
## SECTION D: Collect all results in a clearly named list
###############################################################################

res_list <- list(
  
  ## Senescence (IR vs CTL)
  normal_mock_IR_vs_CTL = res_normal_mock_IR_vs_CTL,
  normal_IAV_IR_vs_CTL  = res_normal_IAV_IR_vs_CTL,
  IPF_mock_IR_vs_CTL    = res_IPF_mock_IR_vs_CTL,
  IPF_IAV_IR_vs_CTL     = res_IPF_IAV_IR_vs_CTL,
  
  ## Infection (IAV vs mock)
  normal_CTL_IAV_vs_mock = res_normal_CTL_IAV_vs_mock,
  normal_IR_IAV_vs_mock  = res_normal_IR_IAV_vs_mock,
  IPF_CTL_IAV_vs_mock    = res_IPF_CTL_IAV_vs_mock,
  IPF_IR_IAV_vs_mock     = res_IPF_IR_IAV_vs_mock,
  
  ## Interaction biology
  normal_infection_modifies_senescence = res_normal_infection_modifies_senescence,
  IPF_infection_modifies_senescence    = res_IPF_infection_modifies_senescence,
  disease_modifies_infection_senescence = res_3way_disease_modifies_infection_senescence
)

saveRDS(res_list,
        file = file.path(deg_dir, "res_3way_clear_labels.rds"))

###############################################################################
# LFC shrinkage (recommended for downstream plots/ranking)
###############################################################################
library(ashr)

shrink_res <- function(dds, res_obj) {
  lfcShrink(dds, res = res_obj, type = "ashr")
}

# Build a parallel list of SHRUNK results (same names as res_list)
res_list_shr <- lapply(res_list, function(r) shrink_res(dds_3way, r))

# Save both versions
saveRDS(res_list,     file = file.path(deg_dir, "res_3way_clear_labels_raw.rds"))
saveRDS(res_list_shr, file = file.path(deg_dir, "res_3way_clear_labels_shrunk_ashr.rds"))



###############################################################################
## SECTION E: Write CSVs (explicit A vs B filenames)
###############################################################################

# --- Senescence ---
write_res_csv(res_list_shr$normal_mock_IR_vs_CTL,
              file.path(deg_dir, "DEG_normal_mock_IR_vs_normal_mock_CTL.csv"), # normal IR mock vs normal CTL mock
              gene_annot)

write_res_csv(res_list_shr$normal_IAV_IR_vs_CTL,
              file.path(deg_dir, "DEG_normal_IAV_IR_vs_normal_IAV_CTL.csv"), # normal IR IAV vs normal CTL IAV
              gene_annot)

write_res_csv(res_list_shr$IPF_mock_IR_vs_CTL,
              file.path(deg_dir, "DEG_IPF_mock_IR_vs_IPF_mock_CTL.csv"), # IPF IR mock vs IPF CTL mock
              gene_annot)

write_res_csv(res_list_shr$IPF_IAV_IR_vs_CTL,
              file.path(deg_dir, "DEG_IPF_IAV_IR_vs_IPF_IAV_CTL.csv"), # IPF IR IAV vs IPF CTL IAV
              gene_annot)

# # --- Infection ---
write_res_csv(res_list_shr$normal_CTL_IAV_vs_mock,
              file.path(deg_dir, "DEG_normal_CTL_IAV_vs_normal_CTL_mock.csv"), # normal CTL IAV vs normal CTL mock
              gene_annot)

write_res_csv(res_list_shr$normal_IR_IAV_vs_mock,
              file.path(deg_dir, "DEG_normal_IR_IAV_vs_normal_IR_mock.csv"), # normal IR IAV vs normal IR mock
              gene_annot)

write_res_csv(res_list_shr$IPF_CTL_IAV_vs_mock,
              file.path(deg_dir, "DEG_IPF_CTL_IAV_vs_IPF_CTL_mock.csv"), # IPF CTL IAV vs IPF CTL mock
              gene_annot)

write_res_csv(res_list_shr$IPF_IR_IAV_vs_mock,
              file.path(deg_dir, "DEG_IPF_IR_IAV_vs_IPF_IR_mock.csv"), # IPF IR IAV vs IPF IR mock
              gene_annot)

# --- Interaction / modulation ---
write_res_csv(res_list_shr$normal_infection_modifies_senescence,
              file.path(deg_dir, "DEG_normal_infection_modifies_senescence.csv"),
              gene_annot)

write_res_csv(res_list_shr$IPF_infection_modifies_senescence,
              file.path(deg_dir, "DEG_IPF_infection_modifies_senescence.csv"),
              gene_annot)

write_res_csv(res_list_shr$disease_modifies_infection_senescence,
              file.path(deg_dir, "DEG_disease_modifies_infection_senescence_3way.csv"),
              gene_annot)

###############################################################################
## SECTION F: Quick summary (how many significant genes per contrast)
###############################################################################

sig_counts <- sapply(res_list_shr, function(x) sum(x$padj < 0.05, na.rm = TRUE))
sig_counts <- sort(sig_counts, decreasing = TRUE)
print(sig_counts)



# volcano plots -----------------------------------------------------------


#volcano plots

library(ggplot2)
library(dplyr)
library(ggrepel)
library(svglite)

# output folders
vol_dir <- file.path(revision_results_dir, "DESeq2", "Volcano")
vol_png_dir <- file.path(vol_dir, "PNG")
vol_svg_dir <- file.path(vol_dir, "SVG")
dir.create(vol_png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(vol_svg_dir, recursive = TRUE, showWarnings = FALSE)

# helper: save both PNG and SVG with same stem
save_plot_dual <- function(plot, out_stem, png_dir, svg_dir,
                           width = 7.5, height = 6.5, dpi = 300) {
  ggsave(
    filename = file.path(png_dir, paste0(out_stem, ".png")),
    plot = plot,
    width = width, height = height, dpi = dpi
  )
  
  svglite::svglite(
    file = file.path(svg_dir, paste0(out_stem, ".svg")),
    width = width, height = height
  )
  print(plot)
  dev.off()
}

# helper: pretty labels for titles/file-facing text
pretty_cmp_label <- function(x) {
  x |>
    gsub("^normal", "healthy", x = _) |>
    gsub("^IPF", "IPF", x = _) |>
    gsub("CTL", "QUI", x = _) |>
    gsub("IR", "SEN", x = _)
}

make_volcano_like_image <- function(res, title, out_stem,
                                    png_dir, svg_dir,
                                    gene_annot = NULL,
                                    p_col = "padj",
                                    p_cut = 0.05,
                                    lfc_cut = 0.5,
                                    label_n = 10,
                                    xlim = NULL) {
  
  if (!missing(p_col) && p_col != "padj") {
    warning("make_volcano_like_image now uses padj only; ignoring p_col = ", p_col)
  }
  
  df <- as.data.frame(res) %>%
    mutate(gene_id = rownames(res))
  
  if (!is.null(gene_annot)) {
    df <- left_join(df, gene_annot, by = "gene_id")
  } else {
    df$gene_name <- NA_character_
  }
  
  pvals <- df$padj
  pvals[pvals == 0] <- .Machine$double.xmin
  
  df <- df %>%
    mutate(
      p_to_plot = pvals,
      neglog10p = -log10(p_to_plot),
      sig_class = case_when(
        !is.na(p_to_plot) & p_to_plot < p_cut & log2FoldChange <= -lfc_cut ~ "Down",
        !is.na(p_to_plot) & p_to_plot < p_cut & log2FoldChange >=  lfc_cut ~ "Up",
        TRUE ~ "Not Sig"
      ),
      label = clean_gene_label(gene_name, gene_id)
    )
  
  label_df <- df %>%
    filter(sig_class != "Not Sig", !is.na(p_to_plot), !is.na(label), label != "") %>%
    arrange(p_to_plot) %>%
    head(label_n)
  
  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10p, color = sig_class)) +
    geom_point(size = 2, alpha = 0.9) +
    scale_color_manual(
      values = c("Down" = "dodgerblue3",
                 "Not Sig" = "grey75",
                 "Up" = "firebrick3")
    ) +
    labs(
      title = title,
      x = "log2FC",
      y = "-log10(padj)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 15)
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = label),
      size = 6.25,
      show.legend = FALSE,
      max.overlaps = 50
    )
  
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim)
  }
  
  save_plot_dual(
    plot = p,
    out_stem = out_stem,
    png_dir = png_dir,
    svg_dir = svg_dir,
    width = 7.5,
    height = 6.5,
    dpi = 300
  )
  
  p
}


# ============================================================
# Revision volcano:
# Baseline IPF QUI mock vs HLF QUI mock
#
# Positive log2FC = higher in IPF
# Negative log2FC = lower in IPF
# ============================================================

baseline_vol_dir <- file.path(
  reviewer_results_dir,
  "DESeq2",
  "Volcano",
  "Baseline_IPF_vs_HLF_QUI_mock"
)

baseline_vol_png_dir <- file.path(baseline_vol_dir, "PNG")
baseline_vol_svg_dir <- file.path(baseline_vol_dir, "SVG")

dir.create(
  baseline_vol_png_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  baseline_vol_svg_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  exists("res_baseline_shr"),
  exists("gene_annot"),
  exists("make_volcano_like_image")
)

# Full log2FC range
p_baseline_full <- make_volcano_like_image(
  res = res_baseline_shr,
  title = "IPF vs HLF, QUI mock (full range)",
  out_stem = "volcano_IPF_vs_HLF_QUI_mock_full",
  png_dir = baseline_vol_png_dir,
  svg_dir = baseline_vol_svg_dir,
  gene_annot = gene_annot,
  p_cut = 0.05,
  lfc_cut = 0.5,
  label_n = 10
)

# Capped display range, matching the original analysis
p_baseline_capped <- make_volcano_like_image(
  res = res_baseline_shr,
  title = "IPF vs HLF, QUI mock (|log2FC| <= 5)",
  out_stem = "volcano_IPF_vs_HLF_QUI_mock_capped",
  png_dir = baseline_vol_png_dir,
  svg_dir = baseline_vol_svg_dir,
  gene_annot = gene_annot,
  p_cut = 0.05,
  lfc_cut = 0.5,
  label_n = 10,
  xlim = c(-5, 5)
)

cat(
  "\nBaseline volcano plots written to:\n",
  normalizePath(baseline_vol_dir),
  "\n"
)



# loop for comparisons ----------------------------------------------------

comparisons <- c(
  "mock_IR_vs_CTL",
  "IAV_IR_vs_CTL",
  "CTL_IAV_vs_mock",
  "IR_IAV_vs_mock"
)

groups <- c("normal", "IPF")

for (grp in groups) {
  for (cmp in comparisons) {
    
    res_name <- paste0(grp, "_", cmp)
    res_obj  <- res_list_shr[[res_name]]
    
    grp_label <- ifelse(grp == "normal", "healthy", grp)
    cmp_label <- cmp |>
      gsub("CTL", "QUI", x = _) |>
      gsub("IR", "SEN", x = _)
    
    # full range
    make_volcano_like_image(
      res = res_obj,
      title = paste0(grp_label, ", ", gsub("_", " ", cmp_label), " (full range)"),
      out_stem = paste0("volcano_", pretty_cmp_label(res_name), "_full"),
      png_dir = vol_png_dir,
      svg_dir = vol_svg_dir,
      gene_annot = gene_annot,
      label_n = 10
    )
    
    # capped range
    make_volcano_like_image(
      res = res_obj,
      title = paste0(grp_label, ", ", gsub("_", " ", cmp_label), " (|log2FC| <= 5)"),
      out_stem = paste0("volcano_", pretty_cmp_label(res_name), "_capped"),
      png_dir = vol_png_dir,
      svg_dir = vol_svg_dir,
      gene_annot = gene_annot,
      label_n = 10,
      xlim = c(-5, 5)
    )
  }
}


# venn diagram loop -------------------------------------------------------

library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(ggvenn)
library(patchwork)
library(svglite)

# Output directory
venn_dir <- file.path(revision_results_dir, "DESeq2", "Venn")
dir.create(venn_dir, recursive = TRUE, showWarnings = FALSE)

# helper: save ggplot as both PNG and SVG
save_plot_dual <- function(plot, out_stem, png_dir, svg_dir,
                           width = 6, height = 5, dpi = 300) {
  ggsave(
    filename = file.path(png_dir, paste0(out_stem, ".png")),
    plot = plot,
    width = width, height = height, dpi = dpi
  )
  
  svglite::svglite(
    file = file.path(svg_dir, paste0(out_stem, ".svg")),
    width = width, height = height
  )
  print(plot)
  dev.off()
}

# helper: rename labels for display / filenames
pretty_cmp_label <- function(x) {
  x <- gsub("^normal", "healthy", x)
  x <- gsub("^Normal", "Healthy", x)
  x <- gsub("CTL", "QUI", x)
  x <- gsub("IR", "SEN", x)
  x
}

# Helper: get DEG gene sets
get_deg_set <- function(res,
                        gene_annot = NULL,
                        direction = c("up", "down"),
                        p_col = "padj",
                        p_cut = 0.05,
                        lfc_cut = 0.5) {
  direction <- match.arg(direction)
  
  df <- as.data.frame(res) %>%
    mutate(gene_id = rownames(res))
  
  if (!is.null(gene_annot)) {
    df <- df %>%
      left_join(gene_annot, by = "gene_id") %>%
      mutate(gene_label = ifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id))
  } else {
    df <- df %>%
      mutate(gene_label = gene_id)
  }
  
  df <- df %>%
    filter(
      !is.na(.data[[p_col]]),
      .data[[p_col]] < p_cut,
      !is.na(log2FoldChange)
    )
  
  if (direction == "up") {
    df <- df %>% filter(log2FoldChange >= lfc_cut)
  } else {
    df <- df %>% filter(log2FoldChange <= -lfc_cut)
  }
  
  unique(df$gene_label)
}


# Helper: 2-set ggvenn with % labels, no axes
make_ggvenn2 <- function(A, B,
                         A_name = "Healthy", B_name = "IPF",
                         fill = c("#b2182b", "#f4a6a6"),
                         title = "") {
  
  venn_data <- setNames(
    list(unique(A), unique(B)),
    c(A_name, B_name)
  )
  
  ggvenn(
    data = venn_data,
    fill_color = fill,
    stroke_size = 0.6,
    set_name_size = 6,
    text_size = 5,
    show_percentage = TRUE,
    show_elements = FALSE
  ) +
    ggtitle(title) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.margin = margin(8, 8, 8, 8)
    )
}


# Thresholds (match volcano)
p_cut   <- 0.05
lfc_cut <- 0.5

A_name <- "Healthy"
B_name <- "IPF"

comparisons <- c(
  "mock_IR_vs_CTL",
  "IAV_IR_vs_CTL",
  "CTL_IAV_vs_mock",
  "IR_IAV_vs_mock"
)

# optional: collect a master summary across all comparisons
master_counts <- list()

for (cmp in comparisons) {
  
  normal_name <- paste0("normal_", cmp)
  ipf_name    <- paste0("IPF_", cmp)
  
  res_normal <- res_list_shr[[normal_name]]
  res_ipf    <- res_list_shr[[ipf_name]]
  
  if (is.null(res_normal) || is.null(res_ipf)) {
    warning("Skipping ", cmp, " because one of these is missing: ",
            normal_name, " or ", ipf_name)
    next
  }
  
  # subfolder per comparison
  cmp_label <- pretty_cmp_label(cmp)
  cmp_dir <- file.path(venn_dir, cmp_label)
  cmp_png_dir <- file.path(cmp_dir, "PNG")
  cmp_svg_dir <- file.path(cmp_dir, "SVG")
  dir.create(cmp_png_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cmp_svg_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  # Build gene sets
  up_normal <- get_deg_set(res_normal, gene_annot, direction = "up",
                           p_cut = p_cut, lfc_cut = lfc_cut)
  up_ipf    <- get_deg_set(res_ipf, gene_annot, direction = "up",
                           p_cut = p_cut, lfc_cut = lfc_cut)
  
  down_normal <- get_deg_set(res_normal, gene_annot, direction = "down",
                             p_cut = p_cut, lfc_cut = lfc_cut)
  down_ipf    <- get_deg_set(res_ipf, gene_annot, direction = "down",
                             p_cut = p_cut, lfc_cut = lfc_cut)
  
  # ensure clean unique character vectors
  up_normal   <- unique(as.character(up_normal))
  up_ipf      <- unique(as.character(up_ipf))
  down_normal <- unique(as.character(down_normal))
  down_ipf    <- unique(as.character(down_ipf))
  
  up_overlap   <- intersect(up_normal, up_ipf)
  down_overlap <- intersect(down_normal, down_ipf)
  
  
  # Export CSVs
  write_csv(tibble(gene = up_normal),
            file.path(cmp_dir, paste0("UP_healthy_", cmp_label, ".csv")))
  write_csv(tibble(gene = up_ipf),
            file.path(cmp_dir, paste0("UP_IPF_", cmp_label, ".csv")))
  write_csv(tibble(gene = up_overlap),
            file.path(cmp_dir, paste0("UP_overlap_healthy_vs_IPF_", cmp_label, ".csv")))
  
  write_csv(tibble(gene = down_normal),
            file.path(cmp_dir, paste0("DOWN_healthy_", cmp_label, ".csv")))
  write_csv(tibble(gene = down_ipf),
            file.path(cmp_dir, paste0("DOWN_IPF_", cmp_label, ".csv")))
  write_csv(tibble(gene = down_overlap),
            file.path(cmp_dir, paste0("DOWN_overlap_healthy_vs_IPF_", cmp_label, ".csv")))
  
  # per-comparison summary CSV
  counts_tbl <- tibble(
    comparison = cmp_label,
    set = c("UP Healthy", "UP IPF", "UP overlap",
            "DOWN Healthy", "DOWN IPF", "DOWN overlap"),
    n   = c(length(up_normal), length(up_ipf), length(up_overlap),
            length(down_normal), length(down_ipf), length(down_overlap))
  )
  write_csv(counts_tbl, file.path(cmp_dir, paste0("venn_counts_summary_", cmp_label, ".csv")))
  print(counts_tbl)
  
  master_counts[[cmp_label]] <- counts_tbl
  
  
  # Titles
  title_cmp <- gsub("_", " ", cmp_label)
  
  # UP venn (red)
  p_up <- make_ggvenn2(
    up_normal, up_ipf,
    A_name = A_name,
    B_name = B_name,
    fill = c("#b2182b", "#f4a6a6"),
    title = paste0("Upregulated DEGs (", title_cmp, ")")
  )
  
  save_plot_dual(
    plot = p_up,
    out_stem = paste0("Venn_UP_", cmp_label),
    png_dir = cmp_png_dir,
    svg_dir = cmp_svg_dir,
    width = 6, height = 5, dpi = 300
  )
  
  # DOWN venn (blue)
  p_down <- make_ggvenn2(
    down_normal, down_ipf,
    A_name = A_name,
    B_name = B_name,
    fill = c("#08519c", "#bdd7e7"),
    title = paste0("Downregulated DEGs (", title_cmp, ")")
  )
  
  save_plot_dual(
    plot = p_down,
    out_stem = paste0("Venn_DOWN_", cmp_label),
    png_dir = cmp_png_dir,
    svg_dir = cmp_svg_dir,
    width = 6, height = 5, dpi = 300
  )
}


# Optional: master summary across all comparisons
if (length(master_counts) > 0) {
  master_counts_tbl <- bind_rows(master_counts)
  write_csv(master_counts_tbl, file.path(venn_dir, "venn_counts_summary_ALL_comparisons.csv"))
}


#venn diagram between comparisons

library(readr)
library(dplyr)
library(tibble)
library(ggvenn)
library(ggplot2)
library(svglite)

venn_dir <- file.path(revision_results_dir, "DESeq2", "Venn")


# helper: save ggplot as both PNG and SVG
save_plot_dual <- function(plot, out_stem, png_dir, svg_dir,
                           width = 7, height = 5.5, dpi = 300) {
  ggsave(
    filename = file.path(png_dir, paste0(out_stem, ".png")),
    plot = plot,
    width = width, height = height, dpi = dpi
  )
  
  svglite::svglite(
    file = file.path(svg_dir, paste0(out_stem, ".svg")),
    width = width, height = height
  )
  print(plot)
  dev.off()
}

# helper: rename labels for display / filenames
pretty_cmp_label <- function(x) {
  x <- gsub("^normal", "healthy", x)
  x <- gsub("^Normal", "Healthy", x)
  x <- gsub("CTL", "QUI", x)
  x <- gsub("IR", "SEN", x)
  x
}

# Helpers
read_gene_csv <- function(path) {
  read_csv(path, show_col_types = FALSE)$gene |>
    unique() |>
    as.character()
}

# Build "only" gene set for a comparison from the saved CSVs
# group_only: "IPF" or "Healthy"
# direction: "UP" or "DOWN"
get_group_only_set <- function(cmp, group_only = c("IPF", "Healthy"), direction = c("UP", "DOWN")) {
  group_only <- match.arg(group_only)
  direction  <- match.arg(direction)
  
  cmp_label <- pretty_cmp_label(cmp)
  cmp_dir <- file.path(venn_dir, cmp_label)
  
  # full sets
  set_healthy <- read_gene_csv(file.path(cmp_dir, paste0(direction, "_healthy_", cmp_label, ".csv")))
  set_ipf     <- read_gene_csv(file.path(cmp_dir, paste0(direction, "_IPF_", cmp_label, ".csv")))
  
  # overlap (Healthy ∩ IPF)
  overlap <- read_gene_csv(file.path(cmp_dir, paste0(direction, "_overlap_healthy_vs_IPF_", cmp_label, ".csv")))
  
  if (group_only == "IPF") {
    setdiff(set_ipf, overlap)
  } else {
    setdiff(set_healthy, overlap)
  }
}

make_cross_venn_and_exports <- function(cmp_A, cmp_B,
                                        group_only = c("IPF", "Healthy"),
                                        direction = c("UP", "DOWN"),
                                        fill = c("#b2182b", "#f4a6a6")) {
  
  group_only <- match.arg(group_only)
  direction  <- match.arg(direction)
  
  cmp_A_label <- pretty_cmp_label(cmp_A)
  cmp_B_label <- pretty_cmp_label(cmp_B)
  
  # Get the two sets being compared
  A <- get_group_only_set(cmp_A, group_only = group_only, direction = direction)
  B <- get_group_only_set(cmp_B, group_only = group_only, direction = direction)
  
  # Bubbles
  left_only  <- setdiff(A, B)
  overlap    <- intersect(A, B)
  right_only <- setdiff(B, A)
  
  # Output folder
  out_dir <- file.path(
    venn_dir,
    "cross",
    paste0(group_only, "_", direction),
    paste0(cmp_A_label, "_vs_", cmp_B_label)
  )
  out_png_dir <- file.path(out_dir, "PNG")
  out_svg_dir <- file.path(out_dir, "SVG")
  dir.create(out_png_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_svg_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save bubble gene lists
  write_csv(
    tibble(gene = sort(left_only)),
    file.path(out_dir, "left_only.csv")
  )
  
  write_csv(
    tibble(gene = sort(overlap)),
    file.path(out_dir, "overlap.csv")
  )
  
  write_csv(
    tibble(gene = sort(right_only)),
    file.path(out_dir, "right_only.csv")
  )
  
  # Counts summary
  counts_tbl <- tibble(
    bubble = c("left_only", "overlap", "right_only"),
    n      = c(length(left_only), length(overlap), length(right_only))
  )
  write_csv(counts_tbl, file.path(out_dir, "counts_summary.csv"))
  
  # Make Venn plot
  venn_data <- setNames(
    list(unique(A), unique(B)),
    c(paste0(group_only, "-only ", direction, ": ", cmp_A_label),
      paste0(group_only, "-only ", direction, ": ", cmp_B_label))
  )
  
  title <- paste0(group_only, "-only ", direction, " DEGs: ", cmp_A_label, " vs ", cmp_B_label)
  
  p <- ggvenn(
    venn_data,
    fill_color = fill,
    stroke_size = 0.6,
    set_name_size = 5,
    text_size = 5,
    show_percentage = TRUE,
    show_elements = FALSE
  ) +
    ggtitle(title) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
    )
  
  save_plot_dual(
    plot = p,
    out_stem = paste0("Venn_", group_only, "only_", direction, "_", cmp_A_label, "_vs_", cmp_B_label),
    png_dir = out_png_dir,
    svg_dir = out_svg_dir,
    width = 7, height = 5.5, dpi = 300
  )
  
  message("Saved: ", out_dir)
  return(invisible(list(plot = p, counts = counts_tbl, out_dir = out_dir)))
}


pairs <- list(
  c("CTL_IAV_vs_mock", "IR_IAV_vs_mock"),
  c("mock_IR_vs_CTL",  "IAV_IR_vs_CTL")
)

for (pair in pairs) {
  cmp_A <- pair[1]
  cmp_B <- pair[2]
  
  for (group_only in c("IPF", "Healthy")) {
    for (direction in c("UP", "DOWN")) {
      
      fill <- if (direction == "UP") c("#b2182b", "#f4a6a6") else c("#08519c", "#bdd7e7")
      
      make_cross_venn_and_exports(
        cmp_A, cmp_B,
        group_only = group_only,
        direction  = direction,
        fill       = fill
      )
    }
  }
}



###############################################################################
# Two heatmaps:
#  (1) Overlap only (Normal ∩ IPF): Up then Down
#  (2) Non-overlap only: Normal-only and IPF-only, Up block then Down block
###############################################################################



# heatmap loop ------------------------------------------------------------

# heatmap loop ------------------------------------------------------------

# Heatmaps (looped): Shared + Disease-specific (full + capped)
library(dplyr)
library(tibble)
library(pheatmap)
library(svglite)


# Output directory
heatmap_dir <- file.path(revision_results_dir, "DESeq2", "Heatmaps")
dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)


# Comparisons to run
comparisons <- c(
  "mock_IR_vs_CTL",
  "IAV_IR_vs_CTL",
  "CTL_IAV_vs_mock",
  "IR_IAV_vs_mock"
)


# Thresholds + params
padj_cut <- 0.05
lfc_cut  <- 0.5
top_n    <- 10
cap_val  <- 5


# Annotation: enforce 1:1 gene_id -> gene_name
gene_annot1 <- gene_annot %>% distinct(gene_id, .keep_all = TRUE)

# helper: figure-friendly gene labels
clean_gene_label <- function(gene_name, gene_id) {
  out <- ifelse(!is.na(gene_name) & gene_name != "", gene_name, NA_character_)
  out[grepl("^ENSG", out)] <- NA_character_
  out
}

# helper: rename labels for display / filenames
pretty_cmp_label <- function(x) {
  x <- gsub("^normal", "healthy", x)
  x <- gsub("^Normal", "Healthy", x)
  x <- gsub("CTL", "QUI", x)
  x <- gsub("IR", "SEN", x)
  x
}


# Helpers
prep_res_id <- function(res, label, gene_annot1) {
  as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    left_join(gene_annot1, by = "gene_id") %>%
    mutate(
      gene_name = clean_gene_label(gene_name, gene_id),
      contrast = label
    ) %>%
    select(gene_id, gene_name, log2FoldChange, padj)
}

pick_top_ud <- function(df, ids, top_n = 10) {
  up <- df %>%
    filter(gene_id %in% ids, log2FoldChange > 0) %>%
    arrange(padj) %>%
    slice_head(n = top_n)
  
  dn <- df %>%
    filter(gene_id %in% ids, log2FoldChange < 0) %>%
    arrange(padj) %>%
    slice_head(n = top_n)
  
  list(up = up, dn = dn)
}

build_mats <- function(gene_ids, df_norm, df_ipf) {
  
  fc_tbl <- tibble(gene_id = gene_ids) %>%
    left_join(df_norm %>% select(gene_id, log2FoldChange, padj, gene_name), by = "gene_id") %>%
    rename(Healthy = log2FoldChange, padj_healthy = padj, gene_name_healthy = gene_name) %>%
    left_join(df_ipf %>% select(gene_id, log2FoldChange, padj, gene_name), by = "gene_id") %>%
    rename(IPF = log2FoldChange, padj_ipf = padj, gene_name_ipf = gene_name)
  
  row_labels <- fc_tbl %>%
    transmute(
      gene_id,
      label = dplyr::coalesce(gene_name_healthy, gene_name_ipf)
    )
  
  fc_tbl2 <- fc_tbl %>%
    filter(gene_id %in% row_labels$gene_id[row_labels$label %in% row_labels$label[!is.na(row_labels$label)]])
  
  row_labels2 <- row_labels %>%
    filter(!is.na(label), label != "")
  
  fc_mat <- fc_tbl2 %>%
    select(gene_id, Healthy, IPF) %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
  rownames(fc_mat) <- make.unique(row_labels2$label)
  
  padj_mat <- fc_tbl2 %>%
    select(gene_id, padj_healthy, padj_ipf) %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
  
  star_mat <- ifelse(
    padj_mat < 0.001, "***",
    ifelse(padj_mat < 0.01, "**",
           ifelse(padj_mat < 0.05, "*", ""))
  )
  
  star_mat[is.na(padj_mat)] <- ""
  star_mat[is.na(fc_mat)]   <- ""
  
  rownames(star_mat) <- rownames(fc_mat)
  colnames(star_mat) <- c("Healthy", "IPF")
  
  list(fc_mat = fc_mat, star_mat = star_mat)
}

save_pheatmap_dual <- function(mat_to_plot, star_mat, title, out_stem, png_dir, svg_dir,
                               breaks_vec, na_col = "white",
                               width = 5, height = 6.5,
                               fontsize_number = 11,
                               fontsize_row = 11,
                               fontsize_col = 13) {
  
  # PNG
  pheatmap(
    mat_to_plot,
    color  = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    breaks = breaks_vec,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = star_mat,
    number_color = "black",
    fontsize_number = fontsize_number,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    main = title,
    na_col = na_col,
    filename = file.path(png_dir, paste0(out_stem, ".png")),
    width = width, height = height
  )
  
  # SVG
  svglite::svglite(
    file = file.path(svg_dir, paste0(out_stem, ".svg")),
    width = width, height = height
  )
  pheatmap(
    mat_to_plot,
    color  = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    breaks = breaks_vec,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = star_mat,
    number_color = "black",
    fontsize_number = fontsize_number,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    main = title,
    na_col = na_col
  )
  dev.off()
}

plot_save_pheatmap_full <- function(fc_mat, star_mat, title, out_stem, png_dir, svg_dir,
                                    na_col = "white", width = 5, height = 6.5) {
  
  if (nrow(fc_mat) == 0) {
    warning("Skipping heatmap with no labeled genes: ", title)
    return(invisible(NULL))
  }
  
  max_abs <- max(abs(fc_mat), na.rm = TRUE)
  
  save_pheatmap_dual(
    mat_to_plot = fc_mat,
    star_mat = star_mat,
    title = title,
    out_stem = out_stem,
    png_dir = png_dir,
    svg_dir = svg_dir,
    breaks_vec = seq(-max_abs, max_abs, length.out = 101),
    na_col = na_col,
    width = width,
    height = height,
    fontsize_number = 11,
    fontsize_row = 11,
    fontsize_col = 13
  )
  
  message("Saved heatmap to: ", file.path(png_dir, paste0(out_stem, ".png")))
}

plot_save_pheatmap_capped <- function(fc_mat, star_mat, title, out_stem, png_dir, svg_dir,
                                      cap_val = 5, na_col = "white",
                                      width = 5, height = 6.5) {
  
  if (nrow(fc_mat) == 0) {
    warning("Skipping heatmap with no labeled genes: ", title)
    return(invisible(NULL))
  }
  
  fc_cap <- pmax(pmin(fc_mat, cap_val), -cap_val)
  
  save_pheatmap_dual(
    mat_to_plot = fc_cap,
    star_mat = star_mat,
    title = title,
    out_stem = out_stem,
    png_dir = png_dir,
    svg_dir = svg_dir,
    breaks_vec = seq(-cap_val, cap_val, length.out = 101),
    na_col = na_col,
    width = width,
    height = height,
    fontsize_number = 11,
    fontsize_row = 11,
    fontsize_col = 13
  )
  
  message("Saved capped heatmap to: ", file.path(png_dir, paste0(out_stem, ".png")))
}


# Loop through comparisons
for (cmp in comparisons) {
  
  normal_key <- paste0("normal_", cmp)
  ipf_key    <- paste0("IPF_", cmp)
  
  res_norm <- res_list_shr[[normal_key]]
  res_ipf  <- res_list_shr[[ipf_key]]
  
  if (is.null(res_norm) || is.null(res_ipf)) {
    warning("Skipping ", cmp, ": missing ", normal_key, " or ", ipf_key)
    next
  }
  
  cmp_label <- pretty_cmp_label(cmp)
  
  # subfolder per comparison
  cmp_dir <- file.path(heatmap_dir, cmp_label)
  cmp_png_dir <- file.path(cmp_dir, "PNG")
  cmp_svg_dir <- file.path(cmp_dir, "SVG")
  dir.create(cmp_png_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cmp_svg_dir, recursive = TRUE, showWarnings = FALSE)
  
  # prepare tables
  df_norm <- prep_res_id(res_norm, "Healthy", gene_annot1)
  df_ipf  <- prep_res_id(res_ipf,  "IPF",    gene_annot1)
  
  # significant sets
  sig_norm <- df_norm %>%
    filter(!is.na(padj), padj < padj_cut, !is.na(log2FoldChange), abs(log2FoldChange) >= lfc_cut)
  
  sig_ipf <- df_ipf %>%
    filter(!is.na(padj), padj < padj_cut, !is.na(log2FoldChange), abs(log2FoldChange) >= lfc_cut)
  
  overlap_ids <- intersect(sig_norm$gene_id, sig_ipf$gene_id)
  norm_only   <- setdiff(sig_norm$gene_id, sig_ipf$gene_id)
  ipf_only    <- setdiff(sig_ipf$gene_id, sig_norm$gene_id)
  
  # genes for overlap heatmap: overlap Up then overlap Down (ranked by Healthy padj)
  ov_norm <- pick_top_ud(sig_norm, overlap_ids, top_n)
  genes_overlap <- c(ov_norm$up$gene_id, ov_norm$dn$gene_id) %>% unique()
  
  # genes for disease-specific heatmap:
  # Healthy-only Up, IPF-only Up, Healthy-only Down, IPF-only Down
  norm_pick <- pick_top_ud(sig_norm, norm_only, top_n)
  ipf_pick  <- pick_top_ud(sig_ipf,  ipf_only,  top_n)
  
  genes_nonoverlap <- c(
    norm_pick$up$gene_id,
    ipf_pick$up$gene_id,
    norm_pick$dn$gene_id,
    ipf_pick$dn$gene_id
  ) %>% unique()
  
  # matrices
  m_overlap    <- build_mats(genes_overlap,    df_norm, df_ipf)
  m_nonoverlap <- build_mats(genes_nonoverlap, df_norm, df_ipf)
  
  # nice title label
  title_cmp <- gsub("_", " ", cmp_label)
  
  # 1) Shared (OVERLAP) — full + capped
  plot_save_pheatmap_full(
    m_overlap$fc_mat,
    m_overlap$star_mat,
    title = paste0("Shared DEGs: ", title_cmp),
    out_stem = paste0("Heatmap_SHARED_", cmp_label, "_FULL"),
    png_dir = cmp_png_dir,
    svg_dir = cmp_svg_dir,
    width = 5, height = 6.5
  )
  
  plot_save_pheatmap_capped(
    m_overlap$fc_mat,
    m_overlap$star_mat,
    title = paste0("Shared DEGs: ", title_cmp),
    out_stem = paste0("Heatmap_SHARED_", cmp_label, "_CAPPED_", cap_val),
    png_dir = cmp_png_dir,
    svg_dir = cmp_svg_dir,
    cap_val = cap_val,
    width = 5, height = 6.5
  )
  
  # 2) Disease-specific (NON-OVERLAP) — full + capped
  plot_save_pheatmap_full(
    m_nonoverlap$fc_mat,
    m_nonoverlap$star_mat,
    title = paste0("Disease-Specific DEGs: ", title_cmp),
    out_stem = paste0("Heatmap_DISEASE_SPECIFIC_", cmp_label, "_FULL"),
    png_dir = cmp_png_dir,
    svg_dir = cmp_svg_dir,
    width = 5, height = 8
  )
  
  plot_save_pheatmap_capped(
    m_nonoverlap$fc_mat,
    m_nonoverlap$star_mat,
    title = paste0("Disease-Specific DEGs: ", title_cmp),
    out_stem = paste0("Heatmap_DISEASE_SPECIFIC_", cmp_label, "_CAPPED_", cap_val),
    png_dir = cmp_png_dir,
    svg_dir = cmp_svg_dir,
    cap_val = cap_val,
    width = 5, height = 8
  )
}


###############################################################################
# GO enrichment (GO:BP) for 6 sets:
# Shared Up/Down, Normal-unique Up/Down, IPF-unique Up/Down
# Output: ONE faceted dotplot with top 5 GO terms per set
###############################################################################


# Go enrichment loop ------------------------------------------------------

# GO:BP dotplots per comparison (looped)
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(forcats)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(svglite)
})


# Comparisons to run
comparisons <- c(
  "mock_IR_vs_CTL",
  "IAV_IR_vs_CTL",
  "CTL_IAV_vs_mock",
  "IR_IAV_vs_mock"
)


# Parameters
padj_cut    <- 0.05
lfc_cut     <- 0.5
top_k_terms <- 5

# Output dir (subfolders per comparison)
go_dir <- file.path(revision_results_dir, "DESeq2", "GO")
dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)


# helper: save ggplot as both PNG and SVG
save_plot_dual <- function(plot, out_stem, png_dir, svg_dir,
                           width = 12, height = 8, dpi = 300) {
  ggsave(
    filename = file.path(png_dir, paste0(out_stem, ".png")),
    plot = plot,
    width = width, height = height, dpi = dpi
  )
  
  svglite::svglite(
    file = file.path(svg_dir, paste0(out_stem, ".svg")),
    width = width, height = height
  )
  print(plot)
  dev.off()
}

# helper: rename labels for display / filenames
pretty_cmp_label <- function(x) {
  x <- gsub("^normal", "healthy", x)
  x <- gsub("^Normal", "Healthy", x)
  x <- gsub("CTL", "QUI", x)
  x <- gsub("IR", "SEN", x)
  x
}


# Helpers
strip_ensembl_version <- function(x) sub("\\.\\d+$", "", x)

get_sig_ids <- function(res, direction = c("up", "down"), padj_cut = 0.05, lfc_cut = 0.5) {
  direction <- match.arg(direction)
  df <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id")
  
  df <- df %>%
    filter(!is.na(padj), padj < padj_cut,
           !is.na(log2FoldChange),
           abs(log2FoldChange) >= lfc_cut)
  
  if (direction == "up") {
    df <- df %>% filter(log2FoldChange > 0)
  } else {
    df <- df %>% filter(log2FoldChange < 0)
  }
  
  unique(df$gene_id)
}

ens2ent <- function(ens_vec) {
  ens_vec <- strip_ensembl_version(ens_vec)
  ens_vec <- unique(ens_vec[!is.na(ens_vec) & ens_vec != ""])
  
  if (length(ens_vec) == 0) return(character(0))
  
  mapped <- clusterProfiler::bitr(
    ens_vec,
    fromType = "ENSEMBL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  unique(mapped$ENTREZID)
}

run_go_bp <- function(ens_ids, set_name, universe_ent) {
  if (length(ens_ids) < 5) return(NULL)
  gene_ent <- ens2ent(ens_ids)
  if (length(gene_ent) < 5) return(NULL)
  
  out <- enrichGO(
    gene          = gene_ent,
    universe      = universe_ent,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  if (is.null(out) || nrow(as.data.frame(out)) == 0) return(NULL)
  
  as.data.frame(out) %>%
    mutate(Set = set_name)
}

make_go_dotplot <- function(go_df, cmp, top_k_terms = 5) {
  # Top K per category
  go_top <- go_df %>%
    filter(!is.na(p.adjust)) %>%
    group_by(Set) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_k_terms) %>%
    ungroup()
  
  # Order categories on x-axis
  set_order <- c(
    "Healthy Unique Down",
    "Shared Down",
    "IPF Unique Down",
    "Shared Up",
    "Healthy Unique Up",
    "IPF Unique Up"
  )
  
  # Keep only sets that exist
  set_order_present <- intersect(set_order, unique(go_top$Set))
  go_top <- go_top %>% mutate(Set = factor(Set, levels = set_order_present))
  
  # Order GO terms by best significance anywhere in plot
  term_order <- go_top %>%
    group_by(Description) %>%
    summarise(best_p = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(best_p) %>%
    pull(Description)
  
  go_top <- go_top %>%
    mutate(Description = factor(Description, levels = rev(term_order)))
  
  ggplot(go_top, aes(x = Set, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust), alpha = 0.95) +
    scale_size_continuous(name = "Gene Count") +
    scale_color_gradient(name = "Adjusted p-value", low = "red", high = "blue") +
    labs(
      title = paste0("GO biological enrichment: ", cmp),
      x = "Comparison Category",
      y = "GO Term"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      panel.grid.major = element_line(color = "grey85"),
      panel.grid.minor = element_blank()
    )
}


# Loop over comparisons
for (cmp in comparisons) {
  
  normal_key <- paste0("normal_", cmp)
  ipf_key    <- paste0("IPF_", cmp)
  
  res_norm <- res_list_shr[[normal_key]]
  res_ipf  <- res_list_shr[[ipf_key]]
  
  if (is.null(res_norm) || is.null(res_ipf)) {
    warning("Skipping ", cmp, ": missing ", normal_key, " or ", ipf_key)
    next
  }
  
  cmp_label <- pretty_cmp_label(cmp)
  
  # Subfolder per comparison
  cmp_dir <- file.path(go_dir, cmp_label)
  cmp_png_dir <- file.path(cmp_dir, "PNG")
  cmp_svg_dir <- file.path(cmp_dir, "SVG")
  dir.create(cmp_png_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cmp_svg_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Gene universe: union of genes tested in both results
  universe_ens <- unique(c(rownames(res_norm), rownames(res_ipf)))
  universe_ent <- ens2ent(universe_ens)
  
  if (length(universe_ent) < 100) {
    warning("Universe too small after mapping for ", cmp, " (n=", length(universe_ent), "). Skipping.")
    next
  }
  
  # Build DEG sets (ENSEMBL)
  norm_up   <- get_sig_ids(res_norm, "up",   padj_cut, lfc_cut)
  norm_down <- get_sig_ids(res_norm, "down", padj_cut, lfc_cut)
  
  ipf_up    <- get_sig_ids(res_ipf,  "up",   padj_cut, lfc_cut)
  ipf_down  <- get_sig_ids(res_ipf,  "down", padj_cut, lfc_cut)
  
  shared_up   <- intersect(norm_up, ipf_up)
  shared_down <- intersect(norm_down, ipf_down)
  
  norm_unique_up   <- setdiff(norm_up, ipf_up)
  norm_unique_down <- setdiff(norm_down, ipf_down)
  
  ipf_unique_up    <- setdiff(ipf_up, norm_up)
  ipf_unique_down  <- setdiff(ipf_down, norm_down)
  
  # Categories to plot
  gene_sets <- list(
    "Shared Up"           = shared_up,
    "Shared Down"         = shared_down,
    "Healthy Unique Up"   = norm_unique_up,
    "Healthy Unique Down" = norm_unique_down,
    "IPF Unique Up"       = ipf_unique_up,
    "IPF Unique Down"     = ipf_unique_down
  )
  
  # Run GO for each category
  go_list <- lapply(names(gene_sets), function(nm) run_go_bp(gene_sets[[nm]], nm, universe_ent))
  go_df <- bind_rows(go_list)
  
  if (nrow(go_df) == 0) {
    warning("No GO terms passed cutoffs for ", cmp, ". Skipping plot.")
    next
  }
  
  # Plot + save
  p <- make_go_dotplot(go_df, cmp_label, top_k_terms = top_k_terms)
  
  save_plot_dual(
    plot = p,
    out_stem = paste0("GO_BP_top", top_k_terms, "_dotplot_", cmp_label),
    png_dir = cmp_png_dir,
    svg_dir = cmp_svg_dir,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  message("Saved GO dot plot to: ", file.path(cmp_png_dir, paste0("GO_BP_top", top_k_terms, "_dotplot_", cmp_label, ".png")))
  
  print(p)
}


# kegg loop ---------------------------------------------------------------

# KEGG dotplots per comparison (looped)
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(forcats)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(svglite)
})


# Comparisons to run
comparisons <- c(
  "mock_IR_vs_CTL",
  "IAV_IR_vs_CTL",
  "CTL_IAV_vs_mock",
  "IR_IAV_vs_mock"
)


# Parameters
padj_cut    <- 0.05
lfc_cut     <- 0.5
top_k_terms <- 5
organism_kegg <- "hsa"  # human

# Output dir (subfolders per comparison)
kegg_dir <- file.path(revision_results_dir, "DESeq2", "KEGG")
dir.create(kegg_dir, recursive = TRUE, showWarnings = FALSE)


# helper: save ggplot as both PNG and SVG
save_plot_dual <- function(plot, out_stem, png_dir, svg_dir,
                           width = 12, height = 8, dpi = 300) {
  ggsave(
    filename = file.path(png_dir, paste0(out_stem, ".png")),
    plot = plot,
    width = width, height = height, dpi = dpi
  )
  
  svglite::svglite(
    file = file.path(svg_dir, paste0(out_stem, ".svg")),
    width = width, height = height
  )
  print(plot)
  dev.off()
}

# helper: rename labels for display / filenames
pretty_cmp_label <- function(x) {
  x <- gsub("^normal", "healthy", x)
  x <- gsub("^Normal", "Healthy", x)
  x <- gsub("CTL", "QUI", x)
  x <- gsub("IR", "SEN", x)
  x
}

# Helpers (same as GO)
strip_ensembl_version <- function(x) sub("\\.\\d+$", "", x)

get_sig_ids <- function(res, direction = c("up", "down"), padj_cut = 0.05, lfc_cut = 0.5) {
  direction <- match.arg(direction)
  df <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id")
  
  df <- df %>%
    filter(!is.na(padj), padj < padj_cut,
           !is.na(log2FoldChange),
           abs(log2FoldChange) >= lfc_cut)
  
  if (direction == "up") {
    df <- df %>% filter(log2FoldChange > 0)
  } else {
    df <- df %>% filter(log2FoldChange < 0)
  }
  
  unique(df$gene_id)
}

ens2ent <- function(ens_vec) {
  ens_vec <- strip_ensembl_version(ens_vec)
  ens_vec <- unique(ens_vec[!is.na(ens_vec) & ens_vec != ""])
  if (length(ens_vec) == 0) return(character(0))
  
  mapped <- clusterProfiler::bitr(
    ens_vec,
    fromType = "ENSEMBL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  unique(mapped$ENTREZID)
}


# KEGG enrichment
run_kegg <- function(ens_ids, set_name, universe_ent, organism = "hsa") {
  if (length(ens_ids) < 5) return(NULL)
  gene_ent <- ens2ent(ens_ids)
  if (length(gene_ent) < 5) return(NULL)
  
  out <- suppressMessages(
    enrichKEGG(
      gene          = gene_ent,
      universe      = universe_ent,
      organism      = organism,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    )
  )
  
  if (is.null(out) || nrow(as.data.frame(out)) == 0) return(NULL)
  
  as.data.frame(out) %>%
    mutate(Set = set_name)
}


# Plotting
make_kegg_dotplot <- function(kegg_df, cmp, top_k_terms = 5) {
  
  kegg_top <- kegg_df %>%
    filter(!is.na(p.adjust)) %>%
    group_by(Set) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_k_terms) %>%
    ungroup()
  
  set_order <- c(
    "Healthy Unique Down",
    "Shared Down",
    "IPF Unique Down",
    "Shared Up",
    "Healthy Unique Up",
    "IPF Unique Up"
  )
  
  # keep only sets that exist in the current data
  set_order_present <- intersect(set_order, unique(kegg_top$Set))
  kegg_top <- kegg_top %>% mutate(Set = factor(Set, levels = set_order_present))
  
  # order pathways by best p.adjust anywhere
  term_order <- kegg_top %>%
    group_by(Description) %>%
    summarise(best_p = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(best_p) %>%
    pull(Description)
  
  kegg_top <- kegg_top %>%
    mutate(Description = factor(Description, levels = rev(term_order)))
  
  ggplot(kegg_top, aes(x = Set, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust), alpha = 0.95) +
    scale_size_continuous(name = "Gene Count") +
    scale_color_gradient(name = "Adjusted p-value", low = "red", high = "blue") +
    labs(
      title = paste0("KEGG pathway enrichment: ", cmp),
      x = "Comparison Category",
      y = "KEGG Pathway"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      panel.grid.major = element_line(color = "grey85"),
      panel.grid.minor = element_blank()
    )
}


# Loop over comparisons
for (cmp in comparisons) {
  
  normal_key <- paste0("normal_", cmp)
  ipf_key    <- paste0("IPF_", cmp)
  
  res_norm <- res_list_shr[[normal_key]]
  res_ipf  <- res_list_shr[[ipf_key]]
  
  if (is.null(res_norm) || is.null(res_ipf)) {
    warning("Skipping ", cmp, ": missing ", normal_key, " or ", ipf_key)
    next
  }
  
  cmp_label <- pretty_cmp_label(cmp)
  
  # Subfolder per comparison
  cmp_dir <- file.path(kegg_dir, cmp_label)
  cmp_png_dir <- file.path(cmp_dir, "PNG")
  cmp_svg_dir <- file.path(cmp_dir, "SVG")
  dir.create(cmp_png_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cmp_svg_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Universe
  universe_ens <- unique(c(rownames(res_norm), rownames(res_ipf)))
  universe_ent <- ens2ent(universe_ens)
  
  if (length(universe_ent) < 100) {
    warning("Universe too small after mapping for ", cmp, " (n=", length(universe_ent), "). Skipping.")
    next
  }
  
  # Build DEG sets (ENSEMBL)
  norm_up   <- get_sig_ids(res_norm, "up",   padj_cut, lfc_cut)
  norm_down <- get_sig_ids(res_norm, "down", padj_cut, lfc_cut)
  
  ipf_up    <- get_sig_ids(res_ipf,  "up",   padj_cut, lfc_cut)
  ipf_down  <- get_sig_ids(res_ipf,  "down", padj_cut, lfc_cut)
  
  shared_up   <- intersect(norm_up, ipf_up)
  shared_down <- intersect(norm_down, ipf_down)
  
  norm_unique_up   <- setdiff(norm_up, ipf_up)
  norm_unique_down <- setdiff(norm_down, ipf_down)
  
  ipf_unique_up    <- setdiff(ipf_up, norm_up)
  ipf_unique_down  <- setdiff(ipf_down, norm_down)
  
  gene_sets <- list(
    "Shared Up"           = shared_up,
    "Shared Down"         = shared_down,
    "Healthy Unique Up"   = norm_unique_up,
    "Healthy Unique Down" = norm_unique_down,
    "IPF Unique Up"       = ipf_unique_up,
    "IPF Unique Down"     = ipf_unique_down
  )
  
  # Run KEGG for each category
  kegg_list <- lapply(names(gene_sets), function(nm) run_kegg(gene_sets[[nm]], nm, universe_ent, organism = organism_kegg))
  kegg_df <- bind_rows(kegg_list)
  
  if (nrow(kegg_df) == 0) {
    warning("No KEGG pathways passed cutoffs for ", cmp, ". Skipping plot.")
    next
  }
  
  # Plot + save
  p <- make_kegg_dotplot(kegg_df, cmp_label, top_k_terms = top_k_terms)
  
  save_plot_dual(
    plot = p,
    out_stem = paste0("KEGG_top", top_k_terms, "_dotplot_", cmp_label),
    png_dir = cmp_png_dir,
    svg_dir = cmp_svg_dir,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  message("Saved KEGG dot plot to: ", file.path(cmp_png_dir, paste0("KEGG_top", top_k_terms, "_dotplot_", cmp_label, ".png")))
  
  print(p)
}


#TF info
#ULM = Univariate Linear Model (in decoupleR): for each TF, it tests whether the genes in that TF’s regulon (targets, with sign mor) show a coordinated shift in your gene-level statistic (your log2FC or Wald stat). The output is a TF activity score + p-value per TF

# #TF enrichment newer ----------------------------------------------------------

###############################################################################
## CollecTRI + decoupleR ULM TF activity inference + outputs you asked for
## - Barplots colored by activated (red) vs repressed (blue)
## - No TF volcano plots; instead CSV of TF scores/pvals/padj
## - TF-target gene “volcano-like” plots for STAT1 and IRF3 per comparison
###############################################################################

# #TF enrichment newer ----------------------------------------------------------

# build CollecTRI network -------------------------------------------------

# suppressPackageStartupMessages({
#   library(readr)
#   library(dplyr)
#   library(tibble)
#   library(decoupleR)
#   library(ggplot2)
#   library(ggrepel)
#   library(AnnotationDbi)
#   library(org.Hs.eg.db)
# })
# 
# tmp2 <- readr::read_tsv("CollecTRI_source.tsv", skip = 1, show_col_types = FALSE)
# 
# collectri_net <- tmp2 %>%
#   transmute(
#     source = as.character(`Transcription Factor (Associated Gene Name)`),
#     target = as.character(`Target Gene (Associated Gene Name)`),
#     mor = case_when(
#       `[DoRothEA_A] Effect` == "Stimulate" ~  1,
#       `[DoRothEA_A] Effect` == "Inhibit"   ~ -1,
#       TRUE ~ NA_real_
#     )
#   ) %>%
#   filter(
#     !is.na(source), source != "",
#     !is.na(target), target != "",
#     !is.na(mor)
#   ) %>%
#   distinct()
# 
# table(collectri_net$mor)
# head(collectri_net)

# build CollecTRI network -------------------------------------------------

# ============================================================
# CollecTRI network used for the original paper analysis
#
# IMPORTANT:
# This block is retained exactly for reproducibility of the
# Figure 4 TF analysis and is also used for reviewer additions.
# ============================================================

collectri_path <- file.path(base_dir, "CollecTRI_source.tsv")
stopifnot(file.exists(collectri_path))

collectri_raw <- readr::read_tsv(
  collectri_path,
  skip = 1,
  show_col_types = FALSE
) %>%
  janitor::clean_names()

to_mor <- function(x) {
  
  x <- tolower(as.character(x))
  
  dplyr::case_when(
    x %in% c(
      "activation",
      "activate",
      "activator",
      "up",
      "positive",
      "+",
      "+1",
      "pos"
    ) ~ 1,
    
    x %in% c(
      "repression",
      "repress",
      "repressor",
      "down",
      "negative",
      "-",
      "-1",
      "neg"
    ) ~ -1,
    
    TRUE ~ NA_real_
  )
}

collectri_net_paper <- collectri_raw %>%
  dplyr::transmute(
    source = transcription_factor_associated_gene_name,
    target = target_gene_associated_gene_name,
    
    # Original paper-analysis behavior:
    # use TRRUST direction when explicitly recognized;
    # otherwise assign a positive interaction.
    mor = dplyr::coalesce(
      to_mor(trrust_regulation),
      1
    ),
    
    likelihood = 1
  ) %>%
  dplyr::mutate(
    mor = ifelse(is.na(mor), 1, mor),
    mor = ifelse(mor >= 0, 1, -1),
    source = as.character(source),
    target = as.character(target)
  ) %>%
  dplyr::filter(
    source != "",
    target != ""
  ) %>%
  dplyr::distinct(
    source,
    target,
    .keep_all = TRUE
  )

# Preserve the object name expected by the existing TF functions
collectri_net <- collectri_net_paper

cat("Paper-version CollecTRI network loaded\n")
cat("Edges:", nrow(collectri_net), "\n")
cat("TFs:", length(unique(collectri_net$source)), "\n")
cat("Targets:", length(unique(collectri_net$target)), "\n")
print(table(collectri_net$mor))


###############################################################################
## CollecTRI + decoupleR ULM TF activity inference + outputs you asked for
## - Barplots colored by activated (red) vs repressed (blue)
## - No TF volcano plots; instead CSV of TF scores/pvals/padj
## - TF-target gene “volcano-like” plots for STAT1 and IRF3 per comparison
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(decoupleR)
  library(ggplot2)
  library(ggrepel)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(svglite)
})

# Force tidyverse verbs to be the ones used
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
summarize <- dplyr::summarize

## ---------------------------------------------------------------------------
## 0) Helpers
## ---------------------------------------------------------------------------

# helper: save ggplot as both PNG and SVG
save_plot_dual <- function(plot, out_stem, png_dir, svg_dir,
                           width = 8, height = 5.5, dpi = 300) {
  ggsave(
    filename = file.path(png_dir, paste0(out_stem, ".png")),
    plot = plot,
    width = width, height = height, dpi = dpi
  )
  
  svglite::svglite(
    file = file.path(svg_dir, paste0(out_stem, ".svg")),
    width = width, height = height
  )
  print(plot)
  dev.off()
}

# helper: rename labels for display / filenames
pretty_cmp_label <- function(x) {
  x <- gsub("^normal", "healthy", x)
  x <- gsub("^Normal", "Healthy", x)
  x <- gsub("CTL", "QUI", x)
  x <- gsub("IR", "SEN", x)
  x
}

# helper: figure-friendly gene labels
clean_gene_label <- function(gene_name, gene_id) {
  out <- ifelse(!is.na(gene_name) & gene_name != "", gene_name, NA_character_)
  out[grepl("^ENSG", out)] <- NA_character_
  out
}

## ---------------------------------------------------------------------------
## 1) Gene ID -> symbol mapper (uses gene_annot + org.Hs fallback)
## ---------------------------------------------------------------------------

strip_ens_version <- function(x) sub("\\.\\d+$", "", x)

map_ids_to_symbols <- function(gene_ids, gene_annot) {
  gene_ids2 <- strip_ens_version(gene_ids)
  
  m1 <- tibble(gene_id = gene_ids2) %>%
    left_join(
      gene_annot %>% mutate(gene_id = strip_ens_version(gene_id)),
      by = "gene_id"
    )
  
  missing <- is.na(m1$gene_name) | m1$gene_name == ""
  if (any(missing)) {
    fb <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(m1$gene_id[missing]),
      columns = c("SYMBOL"),
      keytype = "ENSEMBL"
    ) %>%
      distinct(ENSEMBL, SYMBOL)
    
    m1 <- m1 %>%
      left_join(fb, by = c("gene_id" = "ENSEMBL")) %>%
      mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", SYMBOL, gene_name)) %>%
      dplyr::select(gene_id, gene_name)
  } else {
    m1 <- m1 %>% dplyr::select(gene_id, gene_name)
  }
  
  m1
}

## ---------------------------------------------------------------------------
## 2) TF activity runner (ULM)
## ---------------------------------------------------------------------------

run_tf_ulm <- function(res_obj, net, gene_annot,
                       score_col = c("stat", "log2FoldChange"),
                       min_targets = 5) {
  
  score_col <- match.arg(score_col)
  
  df <- as.data.frame(res_obj) %>%
    tibble::rownames_to_column("gene_id") %>%
    mutate(gene_id = strip_ens_version(gene_id))
  
  v <- df[[score_col]]
  names(v) <- df$gene_id
  v <- v[is.finite(v)]
  
  map_tbl <- map_ids_to_symbols(names(v), gene_annot)
  
  v_tbl <- tibble(gene_id = names(v), score_gene = as.numeric(v)) %>%
    left_join(map_tbl, by = "gene_id") %>%
    filter(!is.na(gene_name), gene_name != "") %>%
    group_by(gene_name) %>%
    summarize(score_gene = mean(score_gene), .groups = "drop")
  
  m <- matrix(v_tbl$score_gene, ncol = 1)
  rownames(m) <- v_tbl$gene_name
  colnames(m) <- "contrast"
  
  acts <- run_ulm(
    mat = m,
    network = net,
    .source = "source",
    .target = "target",
    .mor = "mor",
    minsize = min_targets
  ) %>%
    mutate(
      padj = p.adjust(p_value, method = "BH"),
      direction = ifelse(score >= 0, "Activated", "Repressed")
    ) %>%
    arrange(p_value)
  
  acts
}

## ---------------------------------------------------------------------------
## 3) Plot helper: TF barplot (red=activated, blue=repressed)
## ---------------------------------------------------------------------------

make_tf_barplot_signed <- function(acts, title, out_stem, png_dir, svg_dir, top_n = 20) {
  
  df <- acts %>%
    filter(is.finite(score)) %>%
    mutate(abs_score = abs(score)) %>%
    arrange(desc(abs_score)) %>%
    head(top_n) %>%
    mutate(
      source = factor(source, levels = source[order(score)]),
      direction = ifelse(score >= 0, "Activated", "Repressed")
    )
  
  p <- ggplot(df, aes(x = source, y = score, fill = direction)) +
    geom_col(color = "black", linewidth = 0.3) +
    coord_flip() +
    scale_fill_manual(values = c(Activated = "firebrick3", Repressed = "dodgerblue3")) +
    labs(title = title, x = "TFs", y = "TF activity score") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 11)
    )
  
  save_plot_dual(
    plot = p,
    out_stem = out_stem,
    png_dir = png_dir,
    svg_dir = svg_dir,
    width = 8,
    height = 5.5,
    dpi = 300
  )
  
  p
}



## ---------------------------------------------------------------------------
## 4) TF-target gene “volcano-like” plot for a given TF
##    (shows DE genes that are targets of TF, colored by consistency with mor)
## ---------------------------------------------------------------------------

make_tf_target_plot <- function(tf, res_obj, net, gene_annot,
                                title, out_stem, png_dir, svg_dir, out_csv,
                                label_n = 25,
                                lfc_col = "log2FoldChange",
                                p_col = "padj") {
  
  tf_net <- net %>%
    dplyr::filter(source == tf) %>%
    dplyr::select(source, target, mor)
  
  if (nrow(tf_net) == 0) {
    warning("No targets found in CollecTRI for TF: ", tf)
    return(invisible(NULL))
  }
  
  deg <- as.data.frame(res_obj) %>%
    tibble::rownames_to_column("gene_id") %>%
    mutate(gene_id = strip_ens_version(gene_id)) %>%
    left_join(gene_annot %>% mutate(gene_id = strip_ens_version(gene_id)), by = "gene_id") %>%
    mutate(
      gene = clean_gene_label(gene_name, gene_id),
      gene_for_csv = ifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id)
    ) %>%
    dplyr::filter(!is.na(.data[[lfc_col]]))
  
  deg2 <- deg %>%
    group_by(gene_id, gene, gene_for_csv) %>%
    summarize(
      logfc   = mean(.data[[lfc_col]], na.rm = TRUE),
      p_value = suppressWarnings(min(pvalue, na.rm = TRUE)),
      padj    = suppressWarnings(min(padj, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      p_value = ifelse(is.infinite(p_value) | is.na(p_value), 1, p_value),
      padj    = ifelse(is.infinite(padj)    | is.na(padj),    1, padj)
    )
  
  df <- tf_net %>%
    inner_join(deg2, by = c("target" = "gene_for_csv"))
  
  if (nrow(df) == 0) {
    warning("No overlap between TF targets and DEG genes for TF: ", tf)
    return(invisible(NULL))
  }
  
  df <- df %>%
    mutate(
      p_plot = if (p_col %in% c("pvalue", "p_value")) p_value else padj,
      p_plot = ifelse(p_plot == 0, .Machine$double.xmin, p_plot),
      neglog10p = -log10(p_plot),
      category = case_when(
        mor > 0 & logfc > 0 ~ "Consistent (Act)",
        mor < 0 & logfc < 0 ~ "Consistent (Rep)",
        TRUE ~ "Inconsistent/Unknown"
      )
    )
  
  write_csv(df %>% arrange(p_plot), out_csv)
  
  lab <- df %>%
    filter(!is.na(gene), gene != "") %>%
    arrange(p_plot) %>%
    head(label_n)
  
  p <- ggplot(df, aes(x = logfc, y = neglog10p)) +
    geom_point(aes(size = abs(mor), color = category), alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_color_manual(values = c(
      "Consistent (Act)" = "firebrick3",
      "Consistent (Rep)" = "dodgerblue3",
      "Inconsistent/Unknown" = "grey60"
    )) +
    ggrepel::geom_label_repel(
      data = lab,
      aes(label = gene),
      size = 5,
      max.overlaps = 50,
      show.legend = FALSE
    ) +
    labs(title = title, x = "log2FC (DESeq2)", y = paste0("-log10(", p_col, ")")) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 11)
    )
  
  save_plot_dual(
    plot = p,
    out_stem = out_stem,
    png_dir = png_dir,
    svg_dir = svg_dir,
    width = 8.5,
    height = 6.5,
    dpi = 300
  )
  
  p
}

## ---------------------------------------------------------------------------
## 5) Run across your comparisons and groups
##    - Save TF activity CSVs
##    - Save TF barplots
##    - Save TF-target plots + CSVs for STAT1 and IRF3
## ---------------------------------------------------------------------------

comparisons <- c("mock_IR_vs_CTL","IAV_IR_vs_CTL","CTL_IAV_vs_mock","IR_IAV_vs_mock")
groups <- c("normal","IPF")

# base dirs
# Short paths prevent Windows path-length errors
tf_base_dir <- file.path(
  revision_results_dir,
  "TF_paper"
)

tf_dir_stat <- file.path(tf_base_dir, "stat")
tf_dir_lfc  <- file.path(tf_base_dir, "log2FC")

tf_plot_dir_stat <- file.path(tf_base_dir, "plots_stat")
tf_plot_dir_lfc  <- file.path(tf_base_dir, "plots_log2FC")

tf_target_dir_lfc <- file.path(tf_base_dir, "targets")

# tabular outputs
dir.create(tf_dir_stat, recursive = TRUE, showWarnings = FALSE)
dir.create(tf_dir_lfc,  recursive = TRUE, showWarnings = FALSE)

# figure outputs with PNG/SVG subfolders
tf_plot_dir_stat_png <- file.path(tf_plot_dir_stat, "PNG")
tf_plot_dir_stat_svg <- file.path(tf_plot_dir_stat, "SVG")
tf_plot_dir_lfc_png  <- file.path(tf_plot_dir_lfc, "PNG")
tf_plot_dir_lfc_svg  <- file.path(tf_plot_dir_lfc, "SVG")

dir.create(tf_plot_dir_stat_png, recursive = TRUE, showWarnings = FALSE)
dir.create(tf_plot_dir_stat_svg, recursive = TRUE, showWarnings = FALSE)
dir.create(tf_plot_dir_lfc_png,  recursive = TRUE, showWarnings = FALSE)
dir.create(tf_plot_dir_lfc_svg,  recursive = TRUE, showWarnings = FALSE)

dir.create(tf_target_dir_lfc, recursive = TRUE, showWarnings = FALSE)

all_tf_stat <- list()
all_tf_lfc  <- list()

tfs_of_interest <- c("STAT1", "IRF3")

for (grp in groups) {
  for (cmp in comparisons) {
    
    res_name <- paste0(grp, "_", cmp)
    res_name_pretty <- pretty_cmp_label(res_name)
    
    ## --- 1) Wald stat on RAW results -> TF activity CSV + barplot ---
    res_raw <- res_list[[res_name]]
    if (!is.null(res_raw) && "stat" %in% colnames(as.data.frame(res_raw))) {
      
      acts_stat <- run_tf_ulm(res_raw, collectri_net, gene_annot, score_col = "stat", min_targets = 5)
      all_tf_stat[[res_name]] <- acts_stat
      
      # CSV
      write_csv(acts_stat, file.path(tf_dir_stat, paste0("TF_ULM_", res_name_pretty, "_stat.csv")))
      
      # Signed barplot
      make_tf_barplot_signed(
        acts_stat,
        title = paste0("Top TFs: ", res_name_pretty, " (Wald stat)"),
        out_stem = paste0("TF_bar_", res_name_pretty, "_stat"),
        png_dir = tf_plot_dir_stat_png,
        svg_dir = tf_plot_dir_stat_svg,
        top_n = 20
      )
      
    }
    
    ## --- 2) log2FC on SHRUNK results -> TF activity CSV + barplot + TF-target plots ---
    res_shr <- res_list_shr[[res_name]]
    if (!is.null(res_shr) && "log2FoldChange" %in% colnames(as.data.frame(res_shr))) {
      
      acts_lfc <- run_tf_ulm(res_shr, collectri_net, gene_annot, score_col = "log2FoldChange", min_targets = 5)
      all_tf_lfc[[res_name]] <- acts_lfc
      
      # CSV
      write_csv(acts_lfc, file.path(tf_dir_lfc, paste0("TF_ULM_", res_name_pretty, "_log2FC.csv")))
      
      # Signed barplot
      make_tf_barplot_signed(
        acts_lfc,
        title = paste0("Top TFs: ", res_name_pretty, " (log2FC shrunk)"),
        out_stem = paste0("TF_bar_", res_name_pretty, "_log2FC"),
        png_dir = tf_plot_dir_lfc_png,
        svg_dir = tf_plot_dir_lfc_svg,
        top_n = 20
      )
      
      # TF-target plots for STAT1 + IRF3
      out_sub <- file.path(tf_target_dir_lfc, res_name_pretty)
      out_sub_png <- file.path(out_sub, "PNG")
      out_sub_svg <- file.path(out_sub, "SVG")
      dir.create(out_sub_png, recursive = TRUE, showWarnings = FALSE)
      dir.create(out_sub_svg, recursive = TRUE, showWarnings = FALSE)
      
      for (tf in tfs_of_interest) {
        make_tf_target_plot(
          tf = tf,
          res_obj = res_shr,
          net = collectri_net,
          gene_annot = gene_annot,
          title = paste0(tf, " targets: ", res_name_pretty),
          out_stem = paste0("targets_", tf),          png_dir = out_sub_png,
          svg_dir = out_sub_svg,
          out_csv = file.path(
            out_sub,
            paste0(tf, ".csv")
          ),
          label_n = 25,
          p_col = "padj"
        )
      }
    }
    
    message("Done: ", res_name_pretty)
  }
}

saveRDS(all_tf_stat, file.path(tf_dir_stat, "TF_ULM_all_stat.rds"))
saveRDS(all_tf_lfc,  file.path(tf_dir_lfc,  "TF_ULM_all_log2FC.rds"))

cat("All TF activity inference complete.\n")


# ============================================================
# Reviewer addition:
# TF activity for baseline IPF vs HLF QUI mock comparison
#
# Uses the same ULM functions, thresholds, network, plots,
# and output folders as the original TF comparisons.
# ============================================================

baseline_tf_name <- "IPF_vs_HLF_QUI_mock"

stopifnot(
  exists("res_baseline_raw"),
  exists("res_baseline_shr"),
  exists("collectri_net"),
  exists("run_tf_ulm"),
  exists("make_tf_barplot_signed"),
  exists("make_tf_target_plot")
)

# ------------------------------------------------------------
# 1) Wald-statistic TF activity
#    Primary computational TF-inference result
# ------------------------------------------------------------

baseline_acts_stat <- run_tf_ulm(
  res_obj = res_baseline_raw,
  net = collectri_net,
  gene_annot = gene_annot,
  score_col = "stat",
  min_targets = 5
)

all_tf_stat[[baseline_tf_name]] <- baseline_acts_stat

readr::write_csv(
  baseline_acts_stat,
  file.path(
    tf_dir_stat,
    paste0("TF_ULM_", baseline_tf_name, "_stat.csv")
  )
)

make_tf_barplot_signed(
  acts = baseline_acts_stat,
  title = "Top TFs: IPF vs HLF, QUI mock (Wald stat)",
  out_stem = paste0("TF_bar_", baseline_tf_name, "_stat"),
  png_dir = tf_plot_dir_stat_png,
  svg_dir = tf_plot_dir_stat_svg,
  top_n = 20
)

# ------------------------------------------------------------
# 2) Shrunken-log2FC TF activity
#    Secondary/sensitivity visualization
# ------------------------------------------------------------

baseline_acts_lfc <- run_tf_ulm(
  res_obj = res_baseline_shr,
  net = collectri_net,
  gene_annot = gene_annot,
  score_col = "log2FoldChange",
  min_targets = 5
)

all_tf_lfc[[baseline_tf_name]] <- baseline_acts_lfc

readr::write_csv(
  baseline_acts_lfc,
  file.path(
    tf_dir_lfc,
    paste0("TF_ULM_", baseline_tf_name, "_log2FC.csv")
  )
)

make_tf_barplot_signed(
  acts = baseline_acts_lfc,
  title = "Top TFs: IPF vs HLF, QUI mock (log2FC shrunk)",
  out_stem = paste0("TF_bar_", baseline_tf_name, "_log2FC"),
  png_dir = tf_plot_dir_lfc_png,
  svg_dir = tf_plot_dir_lfc_svg,
  top_n = 20
)

# ------------------------------------------------------------
# 3) STAT1 and IRF3 target-gene plots
# ------------------------------------------------------------

baseline_tf_target_dir <- file.path(
  tf_target_dir_lfc,
  baseline_tf_name
)

baseline_tf_target_png_dir <- file.path(
  baseline_tf_target_dir,
  "PNG"
)

baseline_tf_target_svg_dir <- file.path(
  baseline_tf_target_dir,
  "SVG"
)

dir.create(
  baseline_tf_target_png_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  baseline_tf_target_svg_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

for (tf in tfs_of_interest) {
  
  make_tf_target_plot(
    tf = tf,
    res_obj = res_baseline_shr,
    net = collectri_net,
    gene_annot = gene_annot,
    title = paste0(
      tf,
      " targets: IPF vs HLF, QUI mock"
    ),
    out_stem = paste0(
      "TFtargets_",
      tf,
      "_",
      baseline_tf_name
    ),
    png_dir = baseline_tf_target_png_dir,
    svg_dir = baseline_tf_target_svg_dir,
    out_csv = file.path(
      baseline_tf_target_dir,
      paste0(
        "TFtargets_",
        tf,
        "_",
        baseline_tf_name,
        ".csv"
      )
    ),
    label_n = 25,
    p_col = "padj"
  )
}

# Re-save combined lists so they contain both the original
# comparisons and the new baseline comparison
saveRDS(
  all_tf_stat,
  file.path(tf_dir_stat, "TF_ULM_all_stat.rds")
)

saveRDS(
  all_tf_lfc,
  file.path(tf_dir_lfc, "TF_ULM_all_log2FC.rds")
)

cat(
  "\nBaseline TF activity analysis completed.\n",
  "Results written under:\n",
  normalizePath(tf_base_dir),
  "\n"
)


# ============================================================
# Reviewer addition:
# Antiviral GSEA setup
# ============================================================

required_gsea_packages <- c(
  "fgsea",
  "msigdbr",
  "data.table",
  "dplyr",
  "readr",
  "ggplot2",
  "svglite"
)

missing_gsea_packages <- required_gsea_packages[
  !vapply(
    required_gsea_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_gsea_packages) > 0) {
  stop(
    "Install these packages before continuing: ",
    paste(missing_gsea_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(data.table)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(svglite)
})

# ------------------------------------------------------------
# Output directories
# ------------------------------------------------------------

gsea_dir <- file.path(
  reviewer_results_dir,
  "GSEA_antiviral"
)

gsea_table_dir <- file.path(gsea_dir, "tables")
gsea_geneset_dir <- file.path(gsea_dir, "gene_sets")
gsea_png_dir <- file.path(gsea_dir, "PNG")
gsea_svg_dir <- file.path(gsea_dir, "SVG")

dir.create(gsea_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_geneset_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_svg_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Record package versions and inspect available MSigDB
# collections before retrieving pathways
# ------------------------------------------------------------

cat("\n================ GSEA PACKAGE VERSIONS ================\n")
cat("fgsea:", as.character(packageVersion("fgsea")), "\n")
cat("msigdbr:", as.character(packageVersion("msigdbr")), "\n")

msig_collections <- msigdbr::msigdbr_collections()

cat("\n================ MSIGDB COLLECTION COLUMNS ================\n")
print(names(msig_collections))

cat("\n================ AVAILABLE MSIGDB COLLECTIONS ================\n")
print(as.data.frame(msig_collections), row.names = FALSE)

cat(
  "\nGSEA output directory:\n",
  normalizePath(gsea_dir),
  "\n"
)


# ============================================================
# End-of-session workspace checkpoint
# ============================================================

# checkpoint_dir <- file.path(
#   revision_dir,
#   "checkpoints"
# )
# 
# dir.create(
#   checkpoint_dir,
#   recursive = TRUE,
#   showWarnings = FALSE
# )
# 
# checkpoint_file <- file.path(
#   checkpoint_dir,
#   "Jena_revision_checkpoint_latest.RData"
# )
# 
# save.image(
#   file = checkpoint_file,
#   compress = "xz"
# )
# 
# writeLines(
#   capture.output(sessionInfo()),
#   file.path(checkpoint_dir, "sessionInfo_latest.txt")
# )
# 
# savehistory(
#   file.path(checkpoint_dir, "R_console_history_latest.Rhistory")
# )
# 
# cat(
#   "\nCheckpoint saved to:\n",
#   normalizePath(checkpoint_file),
#   "\n"
# )
# 

#restarting after checkpoint
# revision_dir <- "C:/Users/jwbre/Documents/Campisi/Jena_collab_exps_all fastqs/Jena_collab_revisions"
# 
# checkpoint_file <- file.path(
#   revision_dir,
#   "checkpoints",
#   "Jena_revision_checkpoint_latest.RData"
# )
# 
# stopifnot(file.exists(checkpoint_file))
# 
# load(checkpoint_file)
# 
# cat(
#   "Checkpoint loaded successfully:\n",
#   normalizePath(checkpoint_file),
#   "\n"
# )


# ============================================================
# GSEA step 1:
# Retrieve and save Hallmark, GO:BP, and Reactome gene sets
# ============================================================

stopifnot(
  exists("gsea_geneset_dir"),
  dir.exists(gsea_geneset_dir)
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
  
  stopifnot(
    "gs_name" %in% colnames(tbl),
    "gene_symbol" %in% colnames(tbl)
  )
  
  tbl %>%
    dplyr::select(
      gs_name,
      gene_symbol,
      dplyr::everything()
    ) %>%
    dplyr::filter(
      !is.na(gs_name),
      gs_name != "",
      !is.na(gene_symbol),
      gene_symbol != ""
    ) %>%
    dplyr::distinct(gs_name, gene_symbol, .keep_all = TRUE)
}

# Hallmark
hallmark_tbl <- get_msigdb_collection(
  collection = "H"
)

# GO Biological Process
gobp_tbl <- get_msigdb_collection(
  collection = "C5",
  subcollection = "GO:BP"
)

# Reactome
reactome_tbl <- get_msigdb_collection(
  collection = "C2",
  subcollection = "CP:REACTOME"
)

# Convert long tables into named pathway lists for fgsea
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

# Remove duplicate genes within each pathway
hallmark_sets <- lapply(hallmark_sets, unique)
gobp_sets <- lapply(gobp_sets, unique)
reactome_sets <- lapply(reactome_sets, unique)

# Save both long tables and fgsea-ready pathway lists
readr::write_csv(
  hallmark_tbl,
  file.path(gsea_geneset_dir, "MSigDB_Hallmark.csv")
)

readr::write_csv(
  gobp_tbl,
  file.path(gsea_geneset_dir, "MSigDB_GO_BP.csv")
)

readr::write_csv(
  reactome_tbl,
  file.path(gsea_geneset_dir, "MSigDB_Reactome.csv")
)

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

cat("\n================ GENE-SET COUNTS ================\n")
cat("Hallmark pathways:", length(hallmark_sets), "\n")
cat("GO:BP pathways:", length(gobp_sets), "\n")
cat("Reactome pathways:", length(reactome_sets), "\n")

cat("\n================ UNIQUE GENES ================\n")
cat("Hallmark genes:", length(unique(hallmark_tbl$gene_symbol)), "\n")
cat("GO:BP genes:", length(unique(gobp_tbl$gene_symbol)), "\n")
cat("Reactome genes:", length(unique(reactome_tbl$gene_symbol)), "\n")


# ============================================================
# GSEA step 2:
# Create complete Wald-statistic-ranked gene lists
# ============================================================

gsea_rank_dir <- file.path(
  gsea_dir,
  "ranked_lists"
)

dir.create(
  gsea_rank_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Convert one DESeq2 result into a named, decreasing rank vector
# ------------------------------------------------------------

make_gsea_rank <- function(res_obj, gene_annot) {
  
  rank_df <- as.data.frame(res_obj) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::mutate(
      gene_id = strip_ens_version(gene_id),
      stat = as.numeric(stat)
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
    # Keep the strongest absolute Wald statistic when multiple
    # Ensembl IDs map to the same gene symbol
    dplyr::slice_max(
      order_by = abs(stat),
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(stat))
  
  ranks <- rank_df$stat
  names(ranks) <- rank_df$gene_name
  
  ranks <- sort(
    ranks,
    decreasing = TRUE
  )
  
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

# ------------------------------------------------------------
# Check that the required original contrasts exist
# ------------------------------------------------------------

required_gsea_results <- c(
  "normal_CTL_IAV_vs_mock",
  "IPF_CTL_IAV_vs_mock",
  "normal_IR_IAV_vs_mock",
  "IPF_IR_IAV_vs_mock"
)

stopifnot(
  all(required_gsea_results %in% names(res_list)),
  exists("res_baseline_raw")
)

# ------------------------------------------------------------
# Create ranked lists
# ------------------------------------------------------------

gsea_result_objects <- list(
  healthy_QUI_IAV_vs_mock =
    res_list[["normal_CTL_IAV_vs_mock"]],
  
  IPF_QUI_IAV_vs_mock =
    res_list[["IPF_CTL_IAV_vs_mock"]],
  
  healthy_SEN_IAV_vs_mock =
    res_list[["normal_IR_IAV_vs_mock"]],
  
  IPF_SEN_IAV_vs_mock =
    res_list[["IPF_IR_IAV_vs_mock"]],
  
  IPF_vs_HLF_QUI_mock =
    res_baseline_raw
)

gsea_ranks <- list()

for (contrast_name in names(gsea_result_objects)) {
  
  rank_result <- make_gsea_rank(
    res_obj = gsea_result_objects[[contrast_name]],
    gene_annot = gene_annot
  )
  
  gsea_ranks[[contrast_name]] <- rank_result$ranks
  
  readr::write_csv(
    rank_result$table,
    file.path(
      gsea_rank_dir,
      paste0(
        "ranked_genes_",
        contrast_name,
        "_Wald_stat.csv"
      )
    )
  )
  
  cat(
    contrast_name,
    ":",
    length(rank_result$ranks),
    "ranked genes\n"
  )
}

saveRDS(
  gsea_ranks,
  file.path(
    gsea_rank_dir,
    "GSEA_ranked_lists_Wald_stat.rds"
  )
)

cat("\n================ RANKED-LIST SUMMARY ================\n")

rank_summary <- tibble::tibble(
  contrast = names(gsea_ranks),
  genes = vapply(
    gsea_ranks,
    length,
    integer(1)
  ),
  maximum_stat = vapply(
    gsea_ranks,
    max,
    numeric(1)
  ),
  minimum_stat = vapply(
    gsea_ranks,
    min,
    numeric(1)
  )
)

print(rank_summary)

# ============================================================
# GSEA step 3:
# Run preranked fgsea across all contrasts and collections
# ============================================================

stopifnot(
  exists("gsea_ranks"),
  exists("hallmark_sets"),
  exists("gobp_sets"),
  exists("reactome_sets"),
  exists("gsea_table_dir")
)

gsea_collections <- list(
  Hallmark = hallmark_sets,
  GO_BP = gobp_sets,
  Reactome = reactome_sets
)

gsea_all_results <- list()
gsea_summary_list <- list()

analysis_number <- 0L

for (contrast_name in names(gsea_ranks)) {
  
  ranks <- gsea_ranks[[contrast_name]]
  
  stopifnot(
    is.numeric(ranks),
    !is.null(names(ranks)),
    all(is.finite(ranks))
  )
  
  for (collection_name in names(gsea_collections)) {
    
    analysis_number <- analysis_number + 1L
    
    cat(
      "\nRunning:",
      contrast_name,
      "|",
      collection_name,
      "\n"
    )
    
    # Different but reproducible seed for each analysis
    set.seed(20260729L + analysis_number)
    
    fgsea_result <- fgsea::fgseaMultilevel(
      pathways = gsea_collections[[collection_name]],
      stats = ranks,
      minSize = 15,
      maxSize = 500,
      eps = 0,
      scoreType = "std"
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
    
    # Preserve the original leading-edge list column in the RDS
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
    
    # Convert leading-edge list to a CSV-compatible string
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
      "Tested:",
      nrow(fgsea_result),
      "| FDR < 0.05:",
      nrow(significant_result),
      "\n"
    )
  }
}

# Save complete results with intact leading-edge gene lists
saveRDS(
  gsea_all_results,
  file.path(
    gsea_table_dir,
    "GSEA_all_collections_all_contrasts.rds"
  )
)

gsea_run_summary <- dplyr::bind_rows(
  gsea_summary_list
)

readr::write_csv(
  gsea_run_summary,
  file.path(
    gsea_table_dir,
    "GSEA_run_summary.csv"
  )
)

cat("\n================ GSEA RUN SUMMARY ================\n")
print(gsea_run_summary)


# ============================================================
# Reviewer addition:
# Formal disease-by-IAV interaction contrasts
# ============================================================

# interaction_dir <- file.path(
#   reviewer_results_dir,
#   "DESeq2",
#   "Disease_by_IAV_interactions"
# )
# 
# dir.create(
#   interaction_dir,
#   recursive = TRUE,
#   showWarnings = FALSE
# )
# 
# stopifnot(
#   exists("dds_3way"),
#   exists("gene_annot")
# )
# 
# print(DESeq2::resultsNames(dds_3way))
# 
# # ------------------------------------------------------------
# # 1) QUI cells
# #
# # Difference in the IAV response:
# # IPF versus healthy under QUI conditions
# # ------------------------------------------------------------
# 
# res_interaction_QUI <- DESeq2::results(
#   dds_3way,
#   name = "diseaseIPF.infectionIAV",
#   alpha = 0.05
# )
# 
# # ------------------------------------------------------------
# # 2) SEN cells
# #
# # Difference in the IAV response:
# # IPF versus healthy under SEN conditions
# #
# # SEN interaction =
# # disease:IAV + disease:IAV:condition
# # ------------------------------------------------------------
# 
# res_interaction_SEN <- DESeq2::results(
#   dds_3way,
#   contrast = list(
#     c(
#       "diseaseIPF.infectionIAV",
#       "diseaseIPF.infectionIAV.conditionIR"
#     )
#   ),
#   alpha = 0.05
# )
# 
# # ------------------------------------------------------------
# # Prepare annotation table
# # ------------------------------------------------------------
# 
# gene_annot_interaction <- gene_annot %>%
#   dplyr::mutate(
#     gene_id = strip_ens_version(gene_id)
#   ) %>%
#   dplyr::distinct(
#     gene_id,
#     .keep_all = TRUE
#   )
# 
# prepare_interaction_table <- function(res_obj) {
#   
#   as.data.frame(res_obj) %>%
#     tibble::rownames_to_column("gene_id_original") %>%
#     dplyr::mutate(
#       gene_id = strip_ens_version(gene_id_original)
#     ) %>%
#     dplyr::left_join(
#       gene_annot_interaction,
#       by = "gene_id"
#     ) %>%
#     dplyr::arrange(padj, pvalue)
# }
# 
# interaction_QUI_table <- prepare_interaction_table(
#   res_interaction_QUI
# )
# 
# interaction_SEN_table <- prepare_interaction_table(
#   res_interaction_SEN
# )
# 
# # ------------------------------------------------------------
# # Save results
# # ------------------------------------------------------------
# 
# readr::write_csv(
#   interaction_QUI_table,
#   file.path(
#     interaction_dir,
#     "DESeq2_disease_by_IAV_interaction_QUI.csv"
#   )
# )
# 
# readr::write_csv(
#   interaction_SEN_table,
#   file.path(
#     interaction_dir,
#     "DESeq2_disease_by_IAV_interaction_SEN.csv"
#   )
# )
# 
# saveRDS(
#   res_interaction_QUI,
#   file.path(
#     interaction_dir,
#     "DESeq2_disease_by_IAV_interaction_QUI.rds"
#   )
# )
# 
# saveRDS(
#   res_interaction_SEN,
#   file.path(
#     interaction_dir,
#     "DESeq2_disease_by_IAV_interaction_SEN.rds"
#   )
# )
# 
# # ------------------------------------------------------------
# # Summary
# # ------------------------------------------------------------
# 
# summarize_interaction <- function(res_obj) {
#   
#   df <- as.data.frame(res_obj)
#   
#   tibble::tibble(
#     genes_with_FDR = sum(!is.na(df$padj)),
#     significant_FDR05 = sum(
#       df$padj < 0.05,
#       na.rm = TRUE
#     ),
#     positive_FDR05 = sum(
#       df$padj < 0.05 &
#         df$log2FoldChange > 0,
#       na.rm = TRUE
#     ),
#     negative_FDR05 = sum(
#       df$padj < 0.05 &
#         df$log2FoldChange < 0,
#       na.rm = TRUE
#     )
#   )
# }
# 
# interaction_summary <- dplyr::bind_rows(
#   QUI = summarize_interaction(res_interaction_QUI),
#   SEN = summarize_interaction(res_interaction_SEN),
#   .id = "condition"
# )
# 
# cat(
#   "\n================ INTERACTION SUMMARY ================\n"
# )
# 
# print(interaction_summary)
# 
# 
# # ============================================================
# # GSEA step 5:
# # Create ranked lists for disease-by-IAV interactions
# #
# # Positive statistic:
# # greater IAV response in IPF than healthy
# #
# # Negative statistic:
# # weaker or more negative IAV response in IPF than healthy
# # ============================================================
# 
# stopifnot(
#   exists("res_interaction_QUI"),
#   exists("res_interaction_SEN"),
#   exists("make_gsea_rank"),
#   exists("gsea_rank_dir")
# )
# 
# interaction_gsea_objects <- list(
#   IPF_vs_healthy_IAV_response_QUI = res_interaction_QUI,
#   IPF_vs_healthy_IAV_response_SEN = res_interaction_SEN
# )
# 
# gsea_interaction_ranks <- list()
# 
# for (contrast_name in names(interaction_gsea_objects)) {
#   
#   rank_result <- make_gsea_rank(
#     res_obj = interaction_gsea_objects[[contrast_name]],
#     gene_annot = gene_annot
#   )
#   
#   gsea_interaction_ranks[[contrast_name]] <- rank_result$ranks
#   
#   readr::write_csv(
#     rank_result$table,
#     file.path(
#       gsea_rank_dir,
#       paste0(
#         "ranked_genes_",
#         contrast_name,
#         "_Wald_stat.csv"
#       )
#     )
#   )
#   
#   cat(
#     contrast_name,
#     ":",
#     length(rank_result$ranks),
#     "ranked genes\n"
#   )
# }
# 
# saveRDS(
#   gsea_interaction_ranks,
#   file.path(
#     gsea_rank_dir,
#     "GSEA_interaction_ranked_lists_Wald_stat.rds"
#   )
# )
# 
# interaction_rank_summary <- tibble::tibble(
#   contrast = names(gsea_interaction_ranks),
#   genes = vapply(
#     gsea_interaction_ranks,
#     length,
#     integer(1)
#   ),
#   maximum_stat = vapply(
#     gsea_interaction_ranks,
#     max,
#     numeric(1)
#   ),
#   minimum_stat = vapply(
#     gsea_interaction_ranks,
#     min,
#     numeric(1)
#   )
# )
# 
# cat(
#   "\n================ INTERACTION RANK SUMMARY ================\n"
# )
# 
# print(
#   interaction_rank_summary,
#   width = Inf
# )


# ============================================================
# GSEA plotting step 1:
# Extract antiviral/innate-immune pathways for the three
# manuscript-focused comparisons
# ============================================================

stopifnot(
  exists("gsea_all_results"),
  exists("gsea_dir")
)

gsea_plot_input_dir <- file.path(
  gsea_dir,
  "plot_inputs"
)

dir.create(
  gsea_plot_input_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Comparisons to show in the manuscript-focused plots
plot_contrasts <- c(
  "healthy_SEN_IAV_vs_mock",
  "IPF_SEN_IAV_vs_mock",
  "IPF_vs_HLF_QUI_mock"
)

# Combine the existing GSEA result objects
gsea_plot_all <- dplyr::bind_rows(
  lapply(
    names(gsea_all_results),
    function(result_name) {
      
      name_parts <- strsplit(
        result_name,
        "__",
        fixed = TRUE
      )[[1]]
      
      if (length(name_parts) != 2) {
        return(NULL)
      }
      
      tibble::as_tibble(
        gsea_all_results[[result_name]]
      ) %>%
        dplyr::mutate(
          contrast = name_parts[1],
          collection = name_parts[2]
        )
    }
  )
) %>%
  dplyr::filter(
    contrast %in% plot_contrasts
  )

# Prespecified antiviral and innate-immune search terms
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
    contrast,
    padj,
    dplyr::desc(abs(NES))
  )

# Save a CSV-compatible version
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
    "GSEA_antiviral_pathway_candidates.csv"
  )
)

cat(
  "\n================ ANTIVIRAL CANDIDATE COUNTS ================\n"
)

candidate_summary <- gsea_antiviral_candidates %>%
  dplyr::group_by(
    contrast,
    collection
  ) %>%
  dplyr::summarise(
    candidate_pathways = dplyr::n(),
    significant_FDR05 = sum(
      padj < 0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(
  candidate_summary,
  width = Inf
)

cat(
  "\n================ TOP SIGNIFICANT CANDIDATES ================\n"
)

top_antiviral_candidates <- gsea_antiviral_candidates %>%
  dplyr::filter(
    !is.na(padj),
    padj < 0.05
  ) %>%
  dplyr::group_by(contrast) %>%
  dplyr::slice_min(
    order_by = padj,
    n = 20,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    contrast,
    collection,
    pathway,
    NES,
    padj,
    size
  )

print(
  top_antiviral_candidates,
  n = Inf,
  width = Inf
)


# ============================================================
# GSEA plotting step 2:
# Antiviral/interferon wheel plots for the three
# manuscript-focused comparisons
# ============================================================

stopifnot(
  exists("gsea_plot_all"),
  exists("gsea_dir")
)

gsea_wheel_dir <- file.path(
  gsea_dir,
  "wheel_plots"
)

gsea_wheel_png_dir <- file.path(
  gsea_wheel_dir,
  "PNG"
)

gsea_wheel_svg_dir <- file.path(
  gsea_wheel_dir,
  "SVG"
)

dir.create(
  gsea_wheel_png_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  gsea_wheel_svg_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Fixed pathway panel used for all three comparisons
# ------------------------------------------------------------

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

wheel_contrasts <- tibble::tribble(
  ~contrast, ~contrast_label, ~file_label,
  
  "healthy_SEN_IAV_vs_mock",
  "Healthy SEN: IAV vs mock",
  "healthy_SEN_IAV_vs_mock",
  
  "IPF_SEN_IAV_vs_mock",
  "IPF SEN: IAV vs mock",
  "IPF_SEN_IAV_vs_mock",
  
  "IPF_vs_HLF_QUI_mock",
  "Baseline IPF vs non-PF: QUI mock",
  "IPF_vs_nonPF_QUI_mock"
)

# ------------------------------------------------------------
# Prepare plotting table
# ------------------------------------------------------------

gsea_wheel_data <- gsea_plot_all %>%
  dplyr::inner_join(
    wheel_pathways,
    by = "pathway"
  ) %>%
  dplyr::inner_join(
    wheel_contrasts,
    by = "contrast"
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
      pmax(
        padj,
        .Machine$double.xmin
      )
    ),
    
    # Prevent exceptionally small FDRs from compressing
    # the color scale for all other pathways
    minus_log10_FDR_plot = pmin(
      minus_log10_FDR,
      50
    )
  ) %>%
  dplyr::arrange(
    contrast,
    pathway_order
  )

# Confirm that all 12 pathways were found for each comparison
wheel_coverage <- gsea_wheel_data %>%
  dplyr::count(
    contrast,
    name = "pathways_found"
  )

cat(
  "\n================ WHEEL-PLOT COVERAGE ================\n"
)

print(
  wheel_coverage,
  width = Inf
)

stopifnot(
  all(wheel_coverage$pathways_found == nrow(wheel_pathways))
)

# Save the exact data used for plotting
gsea_wheel_data_csv <- gsea_wheel_data %>%
  dplyr::mutate(
    leadingEdge = vapply(
      leadingEdge,
      function(x) paste(x, collapse = ";"),
      character(1)
    )
  )

readr::write_csv(
  gsea_wheel_data_csv,
  file.path(
    gsea_wheel_dir,
    "GSEA_antiviral_wheel_plot_data.csv"
  )
)

# Shared scales so all three wheels are directly comparable
wheel_y_max <- ceiling(
  max(
    gsea_wheel_data$abs_NES,
    na.rm = TRUE
  ) + 0.5
)

wheel_fill_max <- max(
  gsea_wheel_data$minus_log10_FDR_plot,
  na.rm = TRUE
)


# ============================================================
# GSEA antiviral/interferon wheel plots
#
# Three manuscript-focused comparisons:
# 1. Healthy SEN: IAV vs mock
# 2. IPF SEN: IAV vs mock
# 3. Baseline IPF vs non-PF: QUI mock
#
# Plot encoding:
# - Bar length = absolute NES
# - + / minus symbol = enrichment direction
# - Fill = -log10(FDR), capped at 50
#
# Formatting:
# - Wrapped pathway labels
# - Enlarged fonts
# - Legend below plot
# - No concentric grey circles
# - No radial-axis numbers
# ============================================================

stopifnot(
  exists("gsea_plot_all"),
  exists("gsea_dir")
)

# ------------------------------------------------------------
# Output directories
# ------------------------------------------------------------

gsea_wheel_dir <- file.path(
  gsea_dir,
  "wheel_plots"
)

gsea_wheel_png_dir <- file.path(
  gsea_wheel_dir,
  "PNG"
)

gsea_wheel_svg_dir <- file.path(
  gsea_wheel_dir,
  "SVG"
)

dir.create(
  gsea_wheel_png_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  gsea_wheel_svg_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Fixed antiviral/interferon pathway panel
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Comparisons and plot titles
# ------------------------------------------------------------

wheel_contrasts <- tibble::tribble(
  ~contrast, ~contrast_label, ~file_label,
  
  "healthy_SEN_IAV_vs_mock",
  "Healthy SEN: IAV vs mock",
  "healthy_SEN_IAV_vs_mock",
  
  "IPF_SEN_IAV_vs_mock",
  "IPF SEN: IAV vs mock",
  "IPF_SEN_IAV_vs_mock",
  
  "IPF_vs_HLF_QUI_mock",
  "Baseline IPF vs non-PF: QUI mock",
  "IPF_vs_nonPF_QUI_mock"
)

# ------------------------------------------------------------
# Prepare plotting data
# ------------------------------------------------------------

gsea_wheel_data <- gsea_plot_all %>%
  dplyr::inner_join(
    wheel_pathways,
    by = "pathway"
  ) %>%
  dplyr::inner_join(
    wheel_contrasts,
    by = "contrast"
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
      pmax(
        padj,
        .Machine$double.xmin
      )
    ),
    
    # Cap exceptionally small FDR values so they do not
    # compress the color scale for all other pathways
    minus_log10_FDR_plot = pmin(
      minus_log10_FDR,
      50
    )
  ) %>%
  dplyr::arrange(
    contrast,
    pathway_order
  )

# ------------------------------------------------------------
# Confirm all 12 pathways are present for each comparison
# ------------------------------------------------------------

wheel_coverage <- gsea_wheel_data %>%
  dplyr::count(
    contrast,
    name = "pathways_found"
  )

cat(
  "\n================ WHEEL-PLOT COVERAGE ================\n"
)

print(
  wheel_coverage,
  width = Inf
)

stopifnot(
  all(
    wheel_coverage$pathways_found ==
      nrow(wheel_pathways)
  )
)

# ------------------------------------------------------------
# Save the exact data used for plotting
# ------------------------------------------------------------

gsea_wheel_data_csv <- gsea_wheel_data %>%
  dplyr::mutate(
    leadingEdge = vapply(
      leadingEdge,
      function(x) paste(x, collapse = ";"),
      character(1)
    )
  )

readr::write_csv(
  gsea_wheel_data_csv,
  file.path(
    gsea_wheel_dir,
    "GSEA_antiviral_wheel_plot_data.csv"
  )
)

# ------------------------------------------------------------
# Shared scales for direct comparison across plots
# ------------------------------------------------------------

wheel_y_max <- ceiling(
  max(
    gsea_wheel_data$abs_NES,
    na.rm = TRUE
  ) + 0.5
)

wheel_fill_max <- max(
  gsea_wheel_data$minus_log10_FDR_plot,
  na.rm = TRUE
)

# ------------------------------------------------------------
# Create one wheel plot per comparison
# Enlarged text throughout
# ------------------------------------------------------------

for (i in seq_len(nrow(wheel_contrasts))) {
  
  current_contrast <- wheel_contrasts$contrast[i]
  current_title <- wheel_contrasts$contrast_label[i]
  current_file <- wheel_contrasts$file_label[i]
  
  plot_data <- gsea_wheel_data %>%
    dplyr::filter(
      contrast == current_contrast
    ) %>%
    dplyr::arrange(
      pathway_order
    )
  
  wheel_plot <- ggplot2::ggplot(
    plot_data,
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
      limits = c(
        0,
        wheel_y_max
      ),
      breaks = seq(
        0,
        wheel_y_max,
        by = 1
      ),
      expand = c(
        0,
        0
      )
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      limits = c(
        0,
        wheel_fill_max
      ),
      name = expression(
        -log[10]("FDR")
      ),
      guide = ggplot2::guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(
          6.5,
          "cm"
        ),
        barheight = grid::unit(
          0.55,
          "cm"
        )
      )
    ) +
    ggplot2::labs(
      title = current_title,
      subtitle = paste0(
        "Bar length = |NES|; +/\u2212 indicates enrichment direction",
        "\nColor scale capped at \u2212log10(FDR) = 50"
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(
      base_size = 18
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 25,
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 8
        )
      ),
      
      plot.subtitle = ggplot2::element_text(
        size = 16,
        hjust = 0.5,
        lineheight = 1.08,
        margin = ggplot2::margin(
          b = 18
        )
      ),
      
      axis.text.x = ggplot2::element_text(
        size = 13,
        face = "bold",
        lineheight = 0.92
      ),
      
      # Remove radial-axis numbers
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      
      # Remove concentric circles and radial grid lines
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      
      legend.position = "bottom",
      legend.direction = "horizontal",
      
      legend.title = ggplot2::element_text(
        face = "bold",
        size = 17
      ),
      
      legend.text = ggplot2::element_text(
        size = 15
      ),
      
      legend.margin = ggplot2::margin(
        t = 14
      ),
      
      plot.margin = ggplot2::margin(
        40,
        85,
        40,
        85
      )
    )
  
  # Save PNG
  ggplot2::ggsave(
    filename = file.path(
      gsea_wheel_png_dir,
      paste0(
        "GSEA_wheel_",
        current_file,
        ".png"
      )
    ),
    plot = wheel_plot,
    width = 11,
    height = 11.5,
    dpi = 300
  )
  
  # Save SVG
  svglite::svglite(
    file = file.path(
      gsea_wheel_svg_dir,
      paste0(
        "GSEA_wheel_",
        current_file,
        ".svg"
      )
    ),
    width = 11,
    height = 11.5
  )
  
  print(wheel_plot)
  grDevices::dev.off()
  
  message(
    "Saved enlarged wheel plot: ",
    current_title
  )
}

cat(
  "\nAll three enlarged GSEA wheel plots completed.\n",
  "Output directory:\n",
  normalizePath(gsea_wheel_dir),
  "\n"
)




# ============================================================
# Reviewer 1, Comment 7
# Metabolic/stress pathway candidate extraction
#
# Uses the existing full GSEA results for:
# 1. Healthy SEN: IAV vs mock
# 2. IPF SEN: IAV vs mock
# 3. Baseline IPF vs non-PF: QUI mock
# ============================================================

stopifnot(
  exists("gsea_plot_all"),
  exists("gsea_dir")
)

metabolic_plot_dir <- file.path(
  gsea_dir,
  "metabolic_stress_wheel_plots"
)

metabolic_plot_input_dir <- file.path(
  metabolic_plot_dir,
  "plot_inputs"
)

dir.create(
  metabolic_plot_input_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

metabolic_contrasts <- c(
  "healthy_SEN_IAV_vs_mock",
  "IPF_SEN_IAV_vs_mock",
  "IPF_vs_HLF_QUI_mock"
)

# Broad prespecified terms covering energy metabolism,
# mitochondrial function, redox stress, and proteostasis
metabolic_terms <- c(
  "OXIDATIVE_PHOSPHORYLATION",
  "GLYCOLYSIS",
  "FATTY_ACID",
  "CHOLESTEROL",
  "PEROXISOME",
  "MTORC1",
  "MTOR_SIGNAL",
  "REACTIVE_OXYGEN",
  "OXIDATIVE_STRESS",
  "HYPOXIA",
  "UNFOLDED_PROTEIN",
  "ENDOPLASMIC_RETICULUM_STRESS",
  "TCA_CYCLE",
  "CITRIC_ACID_CYCLE",
  "RESPIRATORY_ELECTRON_TRANSPORT",
  "ELECTRON_TRANSPORT_CHAIN",
  "MITOCHONDRIAL_TRANSLATION",
  "MITOCHONDRIAL_RESPIRATION",
  "CELLULAR_RESPIRATION",
  "ATP_SYNTHESIS",
  "MITOCHONDRIAL_PROTEIN",
  "MITOCHONDRIAL_GENE_EXPRESSION"
)

metabolic_pattern <- paste(
  metabolic_terms,
  collapse = "|"
)

metabolic_candidates <- gsea_plot_all %>%
  dplyr::filter(
    contrast %in% metabolic_contrasts,
    grepl(
      metabolic_pattern,
      pathway,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::mutate(
    minus_log10_FDR = -log10(
      pmax(
        padj,
        .Machine$double.xmin
      )
    ),
    enrichment_direction = dplyr::case_when(
      NES > 0 ~ "Positive",
      NES < 0 ~ "Negative",
      TRUE ~ "Neutral"
    )
  ) %>%
  dplyr::arrange(
    contrast,
    collection,
    padj,
    dplyr::desc(abs(NES))
  )

# Save all candidate pathways
metabolic_candidates_csv <- metabolic_candidates %>%
  dplyr::mutate(
    leadingEdge = vapply(
      leadingEdge,
      function(x) paste(x, collapse = ";"),
      character(1)
    )
  )

readr::write_csv(
  metabolic_candidates_csv,
  file.path(
    metabolic_plot_input_dir,
    "GSEA_metabolic_stress_pathway_candidates.csv"
  )
)

cat(
  "\n================ METABOLIC CANDIDATE COUNTS ================\n"
)

metabolic_candidate_summary <- metabolic_candidates %>%
  dplyr::group_by(
    contrast,
    collection
  ) %>%
  dplyr::summarise(
    candidate_pathways = dplyr::n(),
    significant_FDR05 = sum(
      padj < 0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(
  metabolic_candidate_summary,
  width = Inf
)

# Print Hallmark pathways and the most relevant Reactome pathways.
# These are the strongest candidates for the final wheel panel.
metabolic_candidates_to_review <- metabolic_candidates %>%
  dplyr::filter(
    collection == "Hallmark" |
      grepl(
        paste(
          c(
            "TCA_CYCLE",
            "CITRIC_ACID_CYCLE",
            "RESPIRATORY_ELECTRON_TRANSPORT",
            "MITOCHONDRIAL_TRANSLATION",
            "MITOCHONDRIAL_RESPIRATION",
            "ATP_SYNTHESIS"
          ),
          collapse = "|"
        ),
        pathway,
        ignore.case = TRUE
      )
  ) %>%
  dplyr::select(
    contrast,
    collection,
    pathway,
    NES,
    padj,
    size
  ) %>%
  dplyr::arrange(
    pathway,
    contrast
  )

cat(
  "\n================ CANDIDATES FOR FINAL WHEEL PANEL ================\n"
)

print(
  metabolic_candidates_to_review,
  n = Inf,
  width = Inf
)


# ============================================================
# Reviewer 1, Comment 7
# Metabolic and cellular-stress GSEA wheel plots
#
# Comparisons:
# 1. Healthy SEN: IAV vs mock
# 2. IPF SEN: IAV vs mock
# 3. Baseline IPF vs non-PF: QUI mock
#
# Plot encoding:
# - Bar length = absolute NES
# - + / minus symbol = enrichment direction
# - Fill = -log10(FDR)
#
# Formatting matches the final antiviral wheel plots:
# - wrapped labels
# - enlarged fonts
# - bottom legend
# - no concentric circles
# - no radial-axis numbers
# ============================================================

stopifnot(
  exists("gsea_plot_all"),
  exists("gsea_dir")
)

# ------------------------------------------------------------
# Output directories
# ------------------------------------------------------------

metabolic_wheel_dir <- file.path(
  gsea_dir,
  "metabolic_stress_wheel_plots"
)

metabolic_wheel_png_dir <- file.path(
  metabolic_wheel_dir,
  "PNG"
)

metabolic_wheel_svg_dir <- file.path(
  metabolic_wheel_dir,
  "SVG"
)

dir.create(
  metabolic_wheel_png_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  metabolic_wheel_svg_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Fixed, nonredundant metabolic/stress pathway panel
# ------------------------------------------------------------

metabolic_wheel_pathways <- tibble::tribble(
  ~pathway, ~pathway_label, ~pathway_order,
  
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "Oxidative phosphorylation",
  1,
  
  "HALLMARK_GLYCOLYSIS",
  "Glycolysis",
  2,
  
  "HALLMARK_FATTY_ACID_METABOLISM",
  "Fatty acid metabolism",
  3,
  
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "Cholesterol homeostasis",
  4,
  
  "HALLMARK_PEROXISOME",
  "Peroxisome",
  5,
  
  "HALLMARK_MTORC1_SIGNALING",
  "mTORC1 signaling",
  6,
  
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "Reactive oxygen species pathway",
  7,
  
  "HALLMARK_HYPOXIA",
  "Hypoxia",
  8,
  
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "Unfolded protein response",
  9,
  
  "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
  "Citric acid cycle / TCA cycle",
  10,
  
  "REACTOME_MITOCHONDRIAL_TRANSLATION",
  "Mitochondrial translation",
  11,
  
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "Respiratory electron transport",
  12
)

# ------------------------------------------------------------
# Comparisons and plot titles
# ------------------------------------------------------------

metabolic_wheel_contrasts <- tibble::tribble(
  ~contrast, ~contrast_label, ~file_label,
  
  "healthy_SEN_IAV_vs_mock",
  "Healthy SEN: IAV vs mock",
  "healthy_SEN_IAV_vs_mock",
  
  "IPF_SEN_IAV_vs_mock",
  "IPF SEN: IAV vs mock",
  "IPF_SEN_IAV_vs_mock",
  
  "IPF_vs_HLF_QUI_mock",
  "Baseline IPF vs non-PF: QUI mock",
  "IPF_vs_nonPF_QUI_mock"
)

# ------------------------------------------------------------
# Prepare plotting data
# ------------------------------------------------------------

metabolic_wheel_data <- gsea_plot_all %>%
  dplyr::inner_join(
    metabolic_wheel_pathways,
    by = "pathway"
  ) %>%
  dplyr::inner_join(
    metabolic_wheel_contrasts,
    by = "contrast"
  ) %>%
  dplyr::mutate(
    pathway_label_wrapped = stringr::str_wrap(
      pathway_label,
      width = 22
    ),
    
    pathway_label_wrapped = factor(
      pathway_label_wrapped,
      levels = stringr::str_wrap(
        metabolic_wheel_pathways$pathway_label,
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
      pmax(
        padj,
        .Machine$double.xmin
      )
    )
  ) %>%
  dplyr::arrange(
    contrast,
    pathway_order
  )

# ------------------------------------------------------------
# Confirm all 12 pathways are present in each comparison
# ------------------------------------------------------------

metabolic_wheel_coverage <- metabolic_wheel_data %>%
  dplyr::count(
    contrast,
    name = "pathways_found"
  )

cat(
  "\n================ METABOLIC WHEEL COVERAGE ================\n"
)

print(
  metabolic_wheel_coverage,
  width = Inf
)

stopifnot(
  all(
    metabolic_wheel_coverage$pathways_found ==
      nrow(metabolic_wheel_pathways)
  )
)

# ------------------------------------------------------------
# Save exact plotting data
# ------------------------------------------------------------

metabolic_wheel_data_csv <- metabolic_wheel_data %>%
  dplyr::mutate(
    leadingEdge = vapply(
      leadingEdge,
      function(x) paste(x, collapse = ";"),
      character(1)
    )
  )

readr::write_csv(
  metabolic_wheel_data_csv,
  file.path(
    metabolic_wheel_dir,
    "GSEA_metabolic_stress_wheel_plot_data.csv"
  )
)

# ------------------------------------------------------------
# Shared scales across all three comparisons
# ------------------------------------------------------------

metabolic_wheel_y_max <- ceiling(
  max(
    metabolic_wheel_data$abs_NES,
    na.rm = TRUE
  ) + 0.5
)

metabolic_wheel_fill_max <- ceiling(
  max(
    metabolic_wheel_data$minus_log10_FDR,
    na.rm = TRUE
  )
)

# Prevent a zero-width color scale
metabolic_wheel_fill_max <- max(
  metabolic_wheel_fill_max,
  1
)

# ------------------------------------------------------------
# Create one wheel plot per comparison
# ------------------------------------------------------------

for (i in seq_len(nrow(metabolic_wheel_contrasts))) {
  
  current_contrast <-
    metabolic_wheel_contrasts$contrast[i]
  
  current_title <-
    metabolic_wheel_contrasts$contrast_label[i]
  
  current_file <-
    metabolic_wheel_contrasts$file_label[i]
  
  plot_data <- metabolic_wheel_data %>%
    dplyr::filter(
      contrast == current_contrast
    ) %>%
    dplyr::arrange(
      pathway_order
    )
  
  metabolic_wheel_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = pathway_label_wrapped,
      y = abs_NES,
      fill = minus_log10_FDR
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
      limits = c(
        0,
        metabolic_wheel_y_max
      ),
      breaks = seq(
        0,
        metabolic_wheel_y_max,
        by = 1
      ),
      expand = c(
        0,
        0
      )
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      limits = c(
        0,
        metabolic_wheel_fill_max
      ),
      name = expression(
        -log[10]("FDR")
      ),
      guide = ggplot2::guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(
          6.5,
          "cm"
        ),
        barheight = grid::unit(
          0.55,
          "cm"
        )
      )
    ) +
    ggplot2::labs(
      title = current_title,
      subtitle = paste0(
        "Bar length = |NES|; +/\u2212 indicates enrichment direction",
        "\nColor indicates \u2212log10(FDR)"
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(
      base_size = 18
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 25,
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 8
        )
      ),
      
      plot.subtitle = ggplot2::element_text(
        size = 16,
        hjust = 0.5,
        lineheight = 1.08,
        margin = ggplot2::margin(
          b = 18
        )
      ),
      
      axis.text.x = ggplot2::element_text(
        size = 13,
        face = "bold",
        lineheight = 0.92
      ),
      
      # Remove radial-axis numbers
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      
      # Remove circles and radial grid lines
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      
      legend.position = "bottom",
      legend.direction = "horizontal",
      
      legend.title = ggplot2::element_text(
        face = "bold",
        size = 17
      ),
      
      legend.text = ggplot2::element_text(
        size = 15
      ),
      
      legend.margin = ggplot2::margin(
        t = 14
      ),
      
      plot.margin = ggplot2::margin(
        40,
        85,
        40,
        85
      )
    )
  
  # Save PNG
  ggplot2::ggsave(
    filename = file.path(
      metabolic_wheel_png_dir,
      paste0(
        "GSEA_metabolic_wheel_",
        current_file,
        ".png"
      )
    ),
    plot = metabolic_wheel_plot,
    width = 11,
    height = 11.5,
    dpi = 300
  )
  
  # Save SVG
  svglite::svglite(
    file = file.path(
      metabolic_wheel_svg_dir,
      paste0(
        "GSEA_metabolic_wheel_",
        current_file,
        ".svg"
      )
    ),
    width = 11,
    height = 11.5
  )
  
  print(metabolic_wheel_plot)
  grDevices::dev.off()
  
  message(
    "Saved metabolic wheel plot: ",
    current_title
  )
}

cat(
  "\nAll three metabolic/stress GSEA wheel plots completed.\n",
  "Output directory:\n",
  normalizePath(metabolic_wheel_dir),
  "\n"
)



# ============================================================
# Reviewer 1, Comment 6
# Step 1: Prepare targeted antiviral-gene plotting table
# ============================================================

stopifnot(
  exists("res_list_shr"),
  exists("res_baseline_shr"),
  exists("gene_annot"),
  exists("reviewer_results_dir")
)

comment6_dir <- file.path(
  reviewer_results_dir,
  "Reviewer1_comment6_antiviral_genes"
)

dir.create(
  comment6_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Prespecified antiviral-gene panel
# ------------------------------------------------------------

antiviral_gene_panel <- tibble::tribble(
  ~gene_symbol, ~gene_group, ~gene_order,
  
  "DDX58", "Viral sensing", 1,
  "IFIH1", "Viral sensing", 2,
  "MAVS",  "Viral sensing", 3,
  "TBK1",  "Viral sensing", 4,
  
  "IRF3",  "Transcriptional regulators", 5,
  "IRF7",  "Transcriptional regulators", 6,
  "STAT1", "Transcriptional regulators", 7,
  "STAT2", "Transcriptional regulators", 8,
  
  "IFIT1", "Antiviral effectors", 9,
  "IFIT2", "Antiviral effectors", 10,
  "IFIT3", "Antiviral effectors", 11,
  "MX1",   "Antiviral effectors", 12,
  "OAS1",  "Antiviral effectors", 13,
  "OAS2",  "Antiviral effectors", 14,
  "OAS3",  "Antiviral effectors", 15,
  "OASL",  "Antiviral effectors", 16,
  "ISG15", "Antiviral effectors", 17,
  "RSAD2", "Antiviral effectors", 18,
  "BST2",  "Antiviral effectors", 19,
  "CXCL10","Antiviral effectors", 20,
  
  "USP18", "Negative regulators", 21,
  "SOCS1", "Negative regulators", 22,
  "SOCS3", "Negative regulators", 23,
  "PIAS1", "Negative regulators", 24,
  "PIAS2", "Negative regulators", 25,
  "PIAS3", "Negative regulators", 26,
  "PIAS4", "Negative regulators", 27
)

# ------------------------------------------------------------
# Helper to convert one DESeq2 result to gene-symbol table
# ------------------------------------------------------------

prepare_antiviral_result <- function(
    res_obj,
    comparison_name,
    comparison_label
) {
  
  as.data.frame(res_obj) %>%
    tibble::rownames_to_column("gene_id_original") %>%
    dplyr::mutate(
      gene_id = strip_ens_version(gene_id_original)
    ) %>%
    dplyr::left_join(
      gene_annot %>%
        dplyr::mutate(
          gene_id = strip_ens_version(gene_id)
        ) %>%
        dplyr::distinct(
          gene_id,
          .keep_all = TRUE
        ),
      by = "gene_id"
    ) %>%
    dplyr::filter(
      !is.na(gene_name),
      gene_name != ""
    ) %>%
    dplyr::group_by(gene_name) %>%
    dplyr::slice_max(
      order_by = abs(log2FoldChange),
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      comparison = comparison_name,
      comparison_label = comparison_label,
      gene_symbol = gene_name,
      log2FoldChange = as.numeric(log2FoldChange),
      pvalue = as.numeric(pvalue),
      padj = as.numeric(padj)
    )
}

# ------------------------------------------------------------
# Extract the three manuscript-focused comparisons
# ------------------------------------------------------------

healthy_sen_antiviral <- prepare_antiviral_result(
  res_obj = res_list_shr[["normal_IR_IAV_vs_mock"]],
  comparison_name = "healthy_SEN_IAV_vs_mock",
  comparison_label = "Healthy SEN:\nIAV vs mock"
)

ipf_sen_antiviral <- prepare_antiviral_result(
  res_obj = res_list_shr[["IPF_IR_IAV_vs_mock"]],
  comparison_name = "IPF_SEN_IAV_vs_mock",
  comparison_label = "IPF SEN:\nIAV vs mock"
)

baseline_antiviral <- prepare_antiviral_result(
  res_obj = res_baseline_shr,
  comparison_name = "IPF_vs_nonPF_QUI_mock",
  comparison_label = "Baseline IPF vs non-PF:\nQUI mock"
)

# ------------------------------------------------------------
# Combine and restrict to selected genes
# ------------------------------------------------------------

comment6_antiviral_table <- dplyr::bind_rows(
  healthy_sen_antiviral,
  ipf_sen_antiviral,
  baseline_antiviral
) %>%
  dplyr::inner_join(
    antiviral_gene_panel,
    by = "gene_symbol"
  ) %>%
  dplyr::mutate(
    significance = dplyr::case_when(
      !is.na(padj) & padj < 0.001 ~ "***",
      !is.na(padj) & padj < 0.01  ~ "**",
      !is.na(padj) & padj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    gene_symbol = factor(
      gene_symbol,
      levels = antiviral_gene_panel$gene_symbol
    ),
    gene_group = factor(
      gene_group,
      levels = c(
        "Viral sensing",
        "Transcriptional regulators",
        "Antiviral effectors",
        "Negative regulators"
      )
    ),
    comparison_label = factor(
      comparison_label,
      levels = c(
        "Healthy SEN:\nIAV vs mock",
        "IPF SEN:\nIAV vs mock",
        "Baseline IPF vs non-PF:\nQUI mock"
      )
    )
  ) %>%
  dplyr::arrange(
    gene_order,
    comparison_label
  )

# ------------------------------------------------------------
# Check for missing genes
# ------------------------------------------------------------

genes_found <- unique(
  as.character(comment6_antiviral_table$gene_symbol)
)

genes_missing <- setdiff(
  antiviral_gene_panel$gene_symbol,
  genes_found
)

cat(
  "\n================ ANTIVIRAL GENE COVERAGE ================\n"
)

cat(
  "Genes requested:",
  nrow(antiviral_gene_panel),
  "\n"
)

cat(
  "Genes found:",
  length(genes_found),
  "\n"
)

cat(
  "Genes missing:",
  length(genes_missing),
  "\n"
)

if (length(genes_missing) > 0) {
  print(genes_missing)
}

# ------------------------------------------------------------
# Save the exact plotting table
# ------------------------------------------------------------

readr::write_csv(
  comment6_antiviral_table,
  file.path(
    comment6_dir,
    "Reviewer1_comment6_antiviral_gene_table.csv"
  )
)

# ------------------------------------------------------------
# Print the core reviewer-requested genes
# ------------------------------------------------------------

core_reviewer_genes <- c(
  "DDX58",
  "IFIH1",
  "IRF3",
  "IRF7",
  "STAT1",
  "IFIT1",
  "MX1",
  "OAS1",
  "ISG15",
  "USP18"
)

cat(
  "\n================ CORE ANTIVIRAL GENES ================\n"
)

print(
  comment6_antiviral_table %>%
    dplyr::filter(
      gene_symbol %in% core_reviewer_genes
    ) %>%
    dplyr::select(
      gene_symbol,
      gene_group,
      comparison_label,
      log2FoldChange,
      padj,
      significance
    ),
  n = Inf,
  width = Inf
)



# ============================================================
# Reviewer 1, Comment 6
# Step 2: Rebuild antiviral table with annotation fallback
# and generate the targeted antiviral-gene heatmap
# ============================================================

stopifnot(
  exists("res_list_shr"),
  exists("res_baseline_shr"),
  exists("gene_annot"),
  exists("antiviral_gene_panel"),
  exists("map_ids_to_symbols"),
  exists("strip_ens_version"),
  exists("comment6_dir")
)

comment6_heatmap_dir <- file.path(
  comment6_dir,
  "heatmap"
)

comment6_heatmap_png_dir <- file.path(
  comment6_heatmap_dir,
  "PNG"
)

comment6_heatmap_svg_dir <- file.path(
  comment6_heatmap_dir,
  "SVG"
)

dir.create(
  comment6_heatmap_png_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  comment6_heatmap_svg_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Robust result-to-symbol mapper
# Uses gene_annot first and org.Hs.eg.db as a fallback through
# the previously defined map_ids_to_symbols() function
# ------------------------------------------------------------

prepare_antiviral_result_fixed <- function(
    res_obj,
    comparison_name,
    comparison_label
) {
  
  result_df <- as.data.frame(res_obj) %>%
    tibble::rownames_to_column("gene_id_original") %>%
    dplyr::mutate(
      gene_id = strip_ens_version(gene_id_original),
      log2FoldChange = as.numeric(log2FoldChange),
      pvalue = as.numeric(pvalue),
      padj = as.numeric(padj)
    ) %>%
    dplyr::filter(
      is.finite(log2FoldChange)
    )
  
  symbol_map <- map_ids_to_symbols(
    result_df$gene_id,
    gene_annot
  ) %>%
    dplyr::distinct(
      gene_id,
      .keep_all = TRUE
    )
  
  result_df %>%
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
      order_by = abs(log2FoldChange),
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      comparison = comparison_name,
      comparison_label = comparison_label,
      gene_symbol = gene_name,
      log2FoldChange,
      pvalue,
      padj
    )
}

# ------------------------------------------------------------
# Rebuild the three comparison tables
# ------------------------------------------------------------

healthy_sen_antiviral <- prepare_antiviral_result_fixed(
  res_obj = res_list_shr[["normal_IR_IAV_vs_mock"]],
  comparison_name = "healthy_SEN_IAV_vs_mock",
  comparison_label = "Healthy SEN:\nIAV vs mock"
)

ipf_sen_antiviral <- prepare_antiviral_result_fixed(
  res_obj = res_list_shr[["IPF_IR_IAV_vs_mock"]],
  comparison_name = "IPF_SEN_IAV_vs_mock",
  comparison_label = "IPF SEN:\nIAV vs mock"
)

baseline_antiviral <- prepare_antiviral_result_fixed(
  res_obj = res_baseline_shr,
  comparison_name = "IPF_vs_nonPF_QUI_mock",
  comparison_label = "Baseline IPF vs non-PF:\nQUI mock"
)

# ------------------------------------------------------------
# Combine and restrict to the prespecified panel
# ------------------------------------------------------------

comparison_levels <- c(
  "Healthy SEN:\nIAV vs mock",
  "IPF SEN:\nIAV vs mock",
  "Baseline IPF vs non-PF:\nQUI mock"
)

gene_group_levels <- c(
  "Viral sensing",
  "Transcriptional regulators",
  "Antiviral effectors",
  "Negative regulators"
)

gene_display_levels <- antiviral_gene_panel %>%
  dplyr::mutate(
    gene_display = dplyr::if_else(
      gene_symbol == "DDX58",
      "DDX58 (RIG-I)",
      gene_symbol
    )
  ) %>%
  dplyr::pull(gene_display)

comment6_antiviral_table <- dplyr::bind_rows(
  healthy_sen_antiviral,
  ipf_sen_antiviral,
  baseline_antiviral
) %>%
  dplyr::inner_join(
    antiviral_gene_panel,
    by = "gene_symbol"
  ) %>%
  dplyr::mutate(
    significance = dplyr::case_when(
      !is.na(padj) & padj < 0.001 ~ "***",
      !is.na(padj) & padj < 0.01  ~ "**",
      !is.na(padj) & padj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    
    gene_display = dplyr::if_else(
      gene_symbol == "DDX58",
      "DDX58 (RIG-I)",
      gene_symbol
    ),
    
    gene_display = factor(
      gene_display,
      levels = rev(gene_display_levels)
    ),
    
    gene_group = factor(
      gene_group,
      levels = gene_group_levels
    ),
    
    comparison_label = factor(
      comparison_label,
      levels = comparison_levels
    )
  ) %>%
  dplyr::arrange(
    gene_group,
    gene_order,
    comparison_label
  )

# ------------------------------------------------------------
# Coverage check
# ------------------------------------------------------------

genes_found <- comment6_antiviral_table %>%
  dplyr::pull(gene_symbol) %>%
  unique()

genes_missing <- setdiff(
  antiviral_gene_panel$gene_symbol,
  genes_found
)

cat(
  "\n================ UPDATED ANTIVIRAL GENE COVERAGE ================\n"
)

cat(
  "Genes requested:",
  nrow(antiviral_gene_panel),
  "\n"
)

cat(
  "Genes found:",
  length(genes_found),
  "\n"
)

cat(
  "Genes missing:",
  length(genes_missing),
  "\n"
)

if (length(genes_missing) > 0) {
  print(genes_missing)
}

cat(
  "\n================ DDX58 / RIG-I RESULTS ================\n"
)

print(
  comment6_antiviral_table %>%
    dplyr::filter(
      gene_symbol == "DDX58"
    ) %>%
    dplyr::select(
      comparison_label,
      log2FoldChange,
      padj,
      significance
    ),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# Save the rebuilt table
# ------------------------------------------------------------

readr::write_csv(
  comment6_antiviral_table,
  file.path(
    comment6_dir,
    "Reviewer1_comment6_antiviral_gene_table_batch_corrected.csv"
  )
)

# ------------------------------------------------------------
# Shared symmetric log2FC color range
# ------------------------------------------------------------

heatmap_limit <- max(
  1,
  ceiling(
    max(
      abs(comment6_antiviral_table$log2FoldChange),
      na.rm = TRUE
    )
  )
)

# ------------------------------------------------------------
# Generate heatmap (label/layout-fixed version)
# ------------------------------------------------------------

comment6_heatmap <- ggplot2::ggplot(
  comment6_antiviral_table,
  ggplot2::aes(
    x = comparison_label,
    y = gene_display,
    fill = log2FoldChange
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = significance),
    size = 5.5,
    fontface = "bold"
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(gene_group),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "dodgerblue3",
    mid = "white",
    high = "firebrick3",
    midpoint = 0,
    limits = c(-heatmap_limit, heatmap_limit),
    oob = scales::squish,
    breaks = c(-heatmap_limit, 0, heatmap_limit),
    name = "Shrunken\nlog2FC"
  ) +
  ggplot2::scale_x_discrete(
    guide = ggplot2::guide_axis(n.dodge = 2)
  ) +
  ggplot2::labs(
    title = "Batch-corrected antiviral gene-expression responses",
    subtitle = "* FDR < 0.05; ** FDR < 0.01; *** FDR < 0.001",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::theme(
    plot.title.position = "plot",
    
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 24,
      hjust = 0.5,
      margin = ggplot2::margin(b = 8)
    ),
    
    plot.subtitle = ggplot2::element_text(
      size = 15,
      hjust = 0.5,
      margin = ggplot2::margin(b = 18)
    ),
    
    axis.text.x = ggplot2::element_text(
      size = 14,
      face = "bold",
      lineheight = 0.95,
      margin = ggplot2::margin(t = 10),
      vjust = 1
    ),
    
    axis.text.y = ggplot2::element_text(
      size = 14,
      face = "bold"
    ),
    
    axis.ticks = ggplot2::element_blank(),
    
    panel.grid = ggplot2::element_blank(),
    
    strip.placement = "outside",
    
    strip.text.y.left = ggplot2::element_text(
      angle = 0,
      face = "bold",
      size = 14
    ),
    
    strip.background = ggplot2::element_rect(
      fill = "grey92",
      color = NA
    ),
    
    legend.title = ggplot2::element_text(
      face = "bold",
      size = 14
    ),
    
    legend.text = ggplot2::element_text(
      size = 12
    ),
    
    plot.margin = ggplot2::margin(
      t = 40,
      r = 35,
      b = 35,
      l = 35
    )
  )

# ------------------------------------------------------------
# Save PNG and SVG
# ------------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    comment6_heatmap_png_dir,
    "Reviewer1_comment6_antiviral_gene_heatmap.png"
  ),
  plot = comment6_heatmap,
  width = 11,
  height = 12.5,
  dpi = 300,
  bg = "white"
)

svglite::svglite(
  file = file.path(
    comment6_heatmap_svg_dir,
    "Reviewer1_comment6_antiviral_gene_heatmap.svg"
  ),
  width = 11,
  height = 12.5,
  bg = "white"
)

print(comment6_heatmap)
grDevices::dev.off()

cat(
  "\nAntiviral heatmap completed.\n",
  "Output directory:\n",
  normalizePath(comment6_heatmap_dir),
  "\n"
)


# ============================================================
# Reviewer 1, Comment 6
# Healthy-versus-IPF antiviral-response scatter plot
#
# Uses batch-corrected, shrunken log2FC values.
# This is a descriptive comparison of response magnitudes,
# not a formal disease-by-infection interaction test.
# ============================================================

stopifnot(
  exists("comment6_antiviral_table"),
  exists("comment6_dir")
)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(svglite)
})

# ------------------------------------------------------------
# Output directories
# ------------------------------------------------------------

comment6_scatter_dir <- file.path(
  comment6_dir,
  "scatter"
)

comment6_scatter_png_dir <- file.path(
  comment6_scatter_dir,
  "PNG"
)

comment6_scatter_svg_dir <- file.path(
  comment6_scatter_dir,
  "SVG"
)

dir.create(
  comment6_scatter_png_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  comment6_scatter_svg_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Prepare one row per antiviral gene
# ------------------------------------------------------------

comment6_scatter_data <- comment6_antiviral_table %>%
  dplyr::filter(
    comparison %in% c(
      "healthy_SEN_IAV_vs_mock",
      "IPF_SEN_IAV_vs_mock"
    )
  ) %>%
  dplyr::transmute(
    gene_symbol = as.character(gene_symbol),
    gene_group = as.character(gene_group),
    comparison,
    log2FoldChange,
    padj
  ) %>%
  tidyr::pivot_wider(
    names_from = comparison,
    values_from = c(
      log2FoldChange,
      padj
    ),
    names_glue = "{.value}_{comparison}"
  ) %>%
  dplyr::rename(
    healthy_log2FC =
      log2FoldChange_healthy_SEN_IAV_vs_mock,
    
    IPF_log2FC =
      log2FoldChange_IPF_SEN_IAV_vs_mock,
    
    healthy_padj =
      padj_healthy_SEN_IAV_vs_mock,
    
    IPF_padj =
      padj_IPF_SEN_IAV_vs_mock
  ) %>%
  dplyr::filter(
    is.finite(healthy_log2FC),
    is.finite(IPF_log2FC)
  ) %>%
  dplyr::mutate(
    healthy_significant =
      !is.na(healthy_padj) &
      healthy_padj < 0.05,
    
    IPF_significant =
      !is.na(IPF_padj) &
      IPF_padj < 0.05,
    
    significance_pattern = dplyr::case_when(
      healthy_significant & IPF_significant ~
        "FDR < 0.05 in both",
      
      healthy_significant & !IPF_significant ~
        "FDR < 0.05 in healthy only",
      
      !healthy_significant & IPF_significant ~
        "FDR < 0.05 in IPF only",
      
      TRUE ~
        "Not significant"
    ),
    
    significance_pattern = factor(
      significance_pattern,
      levels = c(
        "FDR < 0.05 in both",
        "FDR < 0.05 in healthy only",
        "FDR < 0.05 in IPF only",
        "Not significant"
      )
    ),
    
    gene_group = factor(
      gene_group,
      levels = c(
        "Viral sensing",
        "Transcriptional regulators",
        "Antiviral effectors",
        "Negative regulators"
      )
    ),
    
    response_difference =
      IPF_log2FC - healthy_log2FC
  ) %>%
  dplyr::arrange(
    gene_group,
    gene_symbol
  )

cat(
  "\nGenes included in scatter plot:",
  nrow(comment6_scatter_data),
  "\n"
)

# ------------------------------------------------------------
# Save the exact scatter-plot data
# ------------------------------------------------------------

readr::write_csv(
  comment6_scatter_data,
  file.path(
    comment6_scatter_dir,
    "Reviewer1_comment6_antiviral_scatter_data.csv"
  )
)

# ------------------------------------------------------------
# Shared plot limits
# ------------------------------------------------------------

scatter_limit <- max(
  1,
  ceiling(
    max(
      abs(
        c(
          comment6_scatter_data$healthy_log2FC,
          comment6_scatter_data$IPF_log2FC
        )
      ),
      na.rm = TRUE
    ) + 0.5
  )
)

# ------------------------------------------------------------
# Descriptive Spearman correlation
# ------------------------------------------------------------

scatter_spearman <- stats::cor.test(
  comment6_scatter_data$healthy_log2FC,
  comment6_scatter_data$IPF_log2FC,
  method = "spearman",
  exact = FALSE
)

scatter_annotation <- paste0(
  "Spearman rho = ",
  sprintf(
    "%.2f",
    unname(scatter_spearman$estimate)
  ),
  "\nDescriptive comparison; not an interaction test"
)

# ------------------------------------------------------------
# Generate scatter plot
# ------------------------------------------------------------

comment6_scatter <- ggplot2::ggplot(
  comment6_scatter_data,
  ggplot2::aes(
    x = healthy_log2FC,
    y = IPF_log2FC
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    color = "grey75",
    linewidth = 0.7
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    color = "grey75",
    linewidth = 0.7
  ) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      color = gene_group,
      shape = significance_pattern
    ),
    size = 4.8,
    stroke = 1
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(
      label = gene_symbol,
      color = gene_group
    ),
    size = 5,
    fontface = "bold",
    box.padding = 0.45,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.size = 0.6,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::annotate(
    "text",
    x = -scatter_limit + 0.45,
    y = scatter_limit - 0.45,
    label = scatter_annotation,
    hjust = 0,
    vjust = 1,
    size = 5
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Viral sensing" = "#F8766D",
      "Transcriptional regulators" = "#7CAE00",
      "Antiviral effectors" = "#00BFC4",
      "Negative regulators" = "#C77CFF"
    ),
    drop = FALSE
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "FDR < 0.05 in both" = 16,
      "FDR < 0.05 in healthy only" = 18,
      "FDR < 0.05 in IPF only" = 17,
      "Not significant" = 15
    ),
    drop = FALSE
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      -scatter_limit,
      scatter_limit
    ),
    breaks = scales::pretty_breaks(
      n = 7
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      -scatter_limit,
      scatter_limit
    ),
    breaks = scales::pretty_breaks(
      n = 7
    )
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(
    title = paste(
      "Batch-corrected antiviral transcriptional responses",
      "in healthy and IPF SEN fibroblasts",
      sep = "\n"
    ),
    subtitle = paste(
      "Points above the dashed line show larger",
      "shrunken log2FC values in IPF"
    ),
    x = paste(
      "Healthy SEN: IAV vs mock",
      "Shrunken log2FC",
      sep = "\n"
    ),
    y = paste(
      "IPF SEN: IAV vs mock",
      "Shrunken log2FC",
      sep = "\n"
    ),
    color = "Gene group",
    shape = "Significance"
  ) +
  ggplot2::theme_classic(
    base_size = 18
  ) +
  ggplot2::theme(
    plot.title.position = "plot",
    
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 24,
      hjust = 0.5,
      lineheight = 1.05,
      margin = ggplot2::margin(
        b = 10
      )
    ),
    
    plot.subtitle = ggplot2::element_text(
      size = 16,
      hjust = 0.5,
      margin = ggplot2::margin(
        b = 20
      )
    ),
    
    axis.title = ggplot2::element_text(
      face = "bold",
      size = 18
    ),
    
    axis.text = ggplot2::element_text(
      size = 15,
      color = "black"
    ),
    
    legend.title = ggplot2::element_text(
      face = "bold",
      size = 16
    ),
    
    legend.text = ggplot2::element_text(
      size = 14
    ),
    
    legend.position = "right",
    
    plot.margin = ggplot2::margin(
      t = 35,
      r = 35,
      b = 30,
      l = 35
    )
  )

# ------------------------------------------------------------
# Save PNG
# ------------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    comment6_scatter_png_dir,
    "Reviewer1_comment6_antiviral_response_scatter.png"
  ),
  plot = comment6_scatter,
  width = 13.5,
  height = 10.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# Save SVG
# ------------------------------------------------------------

svglite::svglite(
  file = file.path(
    comment6_scatter_svg_dir,
    "Reviewer1_comment6_antiviral_response_scatter.svg"
  ),
  width = 13.5,
  height = 10.5,
  bg = "white"
)

print(comment6_scatter)
grDevices::dev.off()

cat(
  "\nAntiviral-response scatter plot completed.\n",
  "Genes plotted:",
  nrow(comment6_scatter_data),
  "\nSpearman rho:",
  round(
    unname(scatter_spearman$estimate),
    3
  ),
  "\nOutput directory:\n",
  normalizePath(comment6_scatter_dir),
  "\n"
)


# ============================================================
# Break-point checkpoint
# ============================================================

stopifnot(
  exists("revision_dir"),
  dir.exists(revision_dir)
)

checkpoint_dir <- file.path(
  revision_dir,
  "checkpoints"
)

dir.create(
  checkpoint_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

checkpoint_time <- format(
  Sys.time(),
  "%Y%m%d_%H%M"
)

checkpoint_latest <- file.path(
  checkpoint_dir,
  "Jena_revision_checkpoint_latest.RData"
)

checkpoint_dated <- file.path(
  checkpoint_dir,
  paste0(
    "Jena_revision_checkpoint_",
    checkpoint_time,
    ".RData"
  )
)

# Save the complete active R environment
save.image(
  file = checkpoint_latest,
  compress = "xz"
)

# Keep a timestamped backup of the same checkpoint
file.copy(
  from = checkpoint_latest,
  to = checkpoint_dated,
  overwrite = TRUE
)

# Record package and session information
writeLines(
  capture.output(sessionInfo()),
  file.path(
    checkpoint_dir,
    paste0(
      "sessionInfo_",
      checkpoint_time,
      ".txt"
    )
  )
)

# Save Console command history
savehistory(
  file.path(
    checkpoint_dir,
    paste0(
      "R_console_history_",
      checkpoint_time,
      ".Rhistory"
    )
  )
)

cat(
  "\nCheckpoint saved successfully.\n",
  "Latest checkpoint:\n",
  normalizePath(checkpoint_latest),
  "\nTimestamped backup:\n",
  normalizePath(checkpoint_dated),
  "\n"
)
