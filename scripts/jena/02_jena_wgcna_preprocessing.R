# ============================================================
# Jena batch-adjusted WGCNA preprocessing and post-WGCNA analysis
# ============================================================
# Primary RNA-seq dataset: GEO GSE334185
#
# Purpose:
#   1) verify aligned tximport counts and metadata
#   2) fit the batch-adjusted DESeq2 model
#   3) generate the batch-corrected, variance-filtered matrix used by PyWGCNA
#   4) perform DEG/module overlap and module enrichment analyses after PyWGCNA
#
# Required local inputs include:
#   Jena_collab_exps_results/txi_rawcounts.rds
#   Jena_collab_exps_results/meta_aligned.rds
#   DEG result files under results_batch/DESeq2/DEGs_3way_clear_labels
#   PyWGCNA module outputs under WGCNA/
#
# Statistical models and filtering thresholds are unchanged.
# ============================================================

# ----------------------------
# Setup
# ----------------------------

# Set this to the local Jena project directory.
project_dir <- "PATH/TO/JENA_PROJECT_DIRECTORY"

if (!dir.exists(project_dir)) {
  stop(
    "Project directory not found. Update 'project_dir' at the top of this script."
  )
}

base_dir <- file.path(project_dir, "WGCNA")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
setwd(base_dir)

jena_results_dir <- file.path(project_dir, "Jena_collab_exps_results")


library(tximport)
library(readr)
library(DESeq2)
library(dplyr)
library(svglite)
library(limma)

txi  <- readRDS(file.path(jena_results_dir, "txi_rawcounts.rds"))
meta <- readRDS(file.path(jena_results_dir, "meta_aligned.rds"))

# Batch-adjusted version: make sure metadata is current and aligned
stopifnot("batch" %in% colnames(meta))
stopifnot(identical(rownames(meta), colnames(txi$counts)))
print(table(meta$batch, useNA = "ifany"))

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
resultsNames(dds_3way)


# normalized counts for WGCNA -------------------------------------------------------

library(matrixStats)

#perform vst normalization with deseq2 output
# Batch-adjusted version:
#   - DESeq2 model includes batch for differential expression
#   - For WGCNA module construction, remove the batch effect from the VST matrix
#     while preserving the disease * infection * condition biology.
Jena.y.vsd <- varianceStabilizingTransformation(dds_3way, fitType = 'local') #line of code might not be needed; but the fitType being local is usually used for smaller datasets
Jena.wpn.vsd <- getVarianceStabilizedData(dds_3way) #gets matrix of transformed counts directly from deseq object

design_wgcna_bc <- model.matrix(
  ~ disease * infection * condition,
  data = meta
)

Jena.wpn.vsd.batch_corrected <- limma::removeBatchEffect(
  Jena.wpn.vsd,
  batch = meta$batch,
  design = design_wgcna_bc
)

Jena.rv.wpn <- rowVars(Jena.wpn.vsd.batch_corrected)
summary(Jena.rv.wpn)

# keep the top 25% most variable genes (above the 75th-percentile variance cutoff)
q75_wpn <- quantile(rowVars(Jena.wpn.vsd.batch_corrected), .75)
expr_normalized <- Jena.wpn.vsd.batch_corrected[Jena.rv.wpn > q75_wpn, ]
dim(expr_normalized)


#optimized for larger datasets  (n>30)
# Jena.y.vsd.2<- vst(y.deseq2, blind=FALSE) #blind = TRUE would ignore the experimental design and assume all the samples are from the same condition
# Jena.vst.counts.2 <- assay(Jena.y.vsd.2)
# 
# # get row-wise variance for the VST-transformed counts
# Jena.rv.vst.2 <- rowVars(Jena.vst.counts.2)
# Jena.rv.vst.q75.2 <- quantile(Jena.rv.vst.2, 0.75)
# expr.vst.q75.2 <- Jena.vst.counts.2[Jena.rv.vst.2 > Jena.rv.vst.q75.2, ]
# dim(expr.vst.q75.2)
#need new name for csv file

# write.csv(expr.vst.q75.2, file= file.path(base_dir,"Jena_RNAseq_DEseq2norm_diffNames_vstNew_q75.csv"), row.names=TRUE) #JO080_count_DEseq2norm.csv


#add sample column as the first column with each variable specifically from deseq
meta_v2 <- meta %>%
  # mutate(sample = paste0(tube_folder_label, "_", disease, "_", condition, "_", infection)) %>%
  mutate(sample = tube_folder_label) %>%
  relocate(sample, .before = 1) #%>% #moves the sample column to the first column, necessary for WGCNA
  #relocate(age, .before = sex) #moves continuous variable from last position (as other columns will be removed in the WGCNA) which is necessary for the WGCNA (can't have a continuous variable last)

# Rename the column headers in expr_normalized to sample in Jena_samples_v2
# only if they match to the folder_label title (shouldnt use folder label title and columns should be sample name now)
colnames(expr_normalized) <- meta_v2$sample[
  match(colnames(expr_normalized), meta_v2$tube_folder_label)
]

#write out new csv with "sample" column necessary for WGCNA
write.csv(meta_v2, file = file.path(base_dir, "Jena_metadata_modified2_batch.csv"), row.names = FALSE)

#rename column headers as "sample" as necessary for WGCNA
# Ensure expr and meta are aligned
stopifnot(ncol(expr_normalized) == nrow(meta_v2))
stopifnot(identical(colnames(expr_normalized), meta_v2$tube_folder_label))

# Set column names to the sample names (already done above, but keep if you want it explicit)
colnames(expr_normalized) <- meta_v2$sample

write.csv(expr_normalized,
          file = file.path(base_dir, "Jena_RNAseq_DEseq2norm_batchCorrected_diffNames_vstTyler_q75.csv"),
          row.names = TRUE)





# post WGCNA gene list and heatmap comparisons ----------------------------

### add in GO enrichment from lists here, not from WGCNA package
### make so the bubble plot format is matching from previous code


#whole loop for csvs with overlaps

# ============================================================
# Step 1 (DEG vs Module overlaps) — UPDATED OUTPUT COLUMNS
# Writes per (DEG file × module):
#   ensembl, symbol, connectivity, module_pval, log2FoldChange, padj
# Notes:
# - DEGs: uses significant genes only (padj<p_cut and |LFC|>=lfc_cut depending on direction)
# - Module genes: no significance cutoff, but we pull connectivity + module_pval from module CSV when present
# - No DEG pvalue column included
# - Robust to AnnotationDbi::select() collisions by using dplyr::select explicitly
# - Robust symbol mapping (strips ENS version, filters to ENSG only, filters to valid keys)
# ============================================================

# Continue using project_dir/base_dir configured at the top of this script.

library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)

# ---------------------------
# Paths
# ---------------------------
deg_dir <- file.path(project_dir, "results_batch", "DESeq2", "DEGs_3way_clear_labels")

# IMPORTANT: point to WGCNA output root
wgcna_root <- file.path(base_dir, "Jena_count_DEseq2_batchCorrected_nonzero_50mod_0.4thresh_robustZ_diffNames")
# Update this path if your new batch-corrected WGCNA run has a different folder name.
stopifnot(dir.exists(wgcna_root))

output_base <- file.path(base_dir, "results_batch")
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

out_overlap_dir <- file.path(output_base, "DEG_vs_Module_overlaps")
dir.create(out_overlap_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Helpers
# ---------------------------
strip_ensembl_version <- function(x) sub("\\.\\d+$", "", x)

# Robust ENS->SYMBOL mapping (won't crash if keys invalid)
ens2sym_tbl <- function(keys) {
  keys <- strip_ensembl_version(keys)
  keys <- unique(keys)
  keys <- keys[!is.na(keys) & keys != ""]
  keys <- keys[grepl("^ENSG\\d+$", keys)]
  
  if (length(keys) == 0) {
    return(tibble(ENSEMBL = character(), SYMBOL = character()))
  }
  
  valid <- keys %in% keys(org.Hs.eg.db, keytype = "ENSEMBL")
  keys_valid <- keys[valid]
  
  if (length(keys_valid) == 0) {
    return(tibble(ENSEMBL = character(), SYMBOL = character()))
  }
  
  suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys     = keys_valid,
      keytype  = "ENSEMBL",
      columns  = c("ENSEMBL", "SYMBOL")
    )
  ) %>%
    as_tibble() %>%
    dplyr::distinct(ENSEMBL, .keep_all = TRUE)
}

# Read DEG csv -> significant genes + stats (ensembl, log2FoldChange, padj)
read_deg_sig_df <- function(csv_path, p_cut = 0.05, lfc_cut = 0.5,
                            direction = c("up","down","either")) {
  direction <- match.arg(direction)
  
  df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  gene_col <- c("gene_id","gene","ensembl","ENSEMBL","X") |> keep(~ .x %in% names(df))
  if (length(gene_col) == 0) stop("Couldn't find a gene id column in: ", basename(csv_path))
  gene_col <- gene_col[1]
  
  if (!("padj" %in% names(df))) stop("padj column not found in: ", basename(csv_path))
  if (!("log2FoldChange" %in% names(df))) stop("log2FoldChange column not found in: ", basename(csv_path))
  
  df <- df %>%
    mutate(ensembl = strip_ensembl_version(.data[[gene_col]])) %>%
    filter(!is.na(padj), padj < p_cut, !is.na(log2FoldChange))
  
  if (direction == "up") {
    df <- df %>% filter(log2FoldChange >= lfc_cut)
  } else if (direction == "down") {
    df <- df %>% filter(log2FoldChange <= -lfc_cut)
  } else {
    df <- df %>% filter(abs(log2FoldChange) >= lfc_cut)
  }
  
  df %>%
    dplyr::select(ensembl, log2FoldChange, padj) %>%
    dplyr::distinct(ensembl, .keep_all = TRUE)
}

# Read module CSV -> per-gene (ensembl, connectivity, module_pval)
# read_module_df <- function(module_csv_path) {
#   
#   mdf <- read.csv(module_csv_path, stringsAsFactors = FALSE, check.names = FALSE)
#   
#   # Fix bad/blank/NA column names
#   cn <- colnames(mdf)
#   cn[is.na(cn) | cn == ""] <- paste0("X", which(is.na(cn) | cn == ""))
#   colnames(mdf) <- make.unique(cn)
#   
#   # Pick gene id column
#   gene_col <- c("ensembl","ENSEMBL","gene_id","X","var_names","id") |> keep(~ .x %in% colnames(mdf))
#   if (length(gene_col) == 0) gene_col <- colnames(mdf)[1] else gene_col <- gene_col[1]
#   
#   mdf <- mdf %>%
#     dplyr::mutate(ensembl = strip_ensembl_version(.data[[gene_col]]))
#   
#   # Connectivity column (auto-detect)
#   conn_col <- c("kTotal","kWithin","connectivity","kME","intramodularConnectivity") |> keep(~ .x %in% colnames(mdf))
#   if (length(conn_col) == 0) {
#     mdf$connectivity <- NA_real_
#   } else {
#     mdf$connectivity <- suppressWarnings(as.numeric(mdf[[conn_col[1]]]))
#   }
#   
#   # Module p-value column (auto-detect)
#   pval_candidates <- c("p-val","pval","p_val","p.value","P.Value","PVAL","pvalue","p_value","PValue")
#   pval_col <- pval_candidates[pval_candidates %in% colnames(mdf)]
#   if (length(pval_col) == 0) {
#     mdf$module_pval <- NA_real_
#   } else {
#     mdf$module_pval <- suppressWarnings(as.numeric(mdf[[pval_col[1]]]))
#   }
#   
#   mdf %>%
#     dplyr::select(ensembl, connectivity, module_pval) %>%
#     dplyr::distinct(ensembl, .keep_all = TRUE)
# }

read_module_df <- function(module_csv_path) {
  
  mdf <- read.csv(module_csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Fix bad/blank/NA column names
  cn <- colnames(mdf)
  cn[is.na(cn) | cn == ""] <- paste0("X", which(is.na(cn) | cn == ""))
  colnames(mdf) <- make.unique(cn)
  
  # Pick gene id column
  gene_candidates <- c("ensembl","ENSEMBL","gene_id","var_names","id","Unnamed: 0","X","X1")
  gene_col <- gene_candidates[gene_candidates %in% colnames(mdf)]
  if (length(gene_col) == 0) gene_col <- colnames(mdf)[1] else gene_col <- gene_col[1]
  
  mdf <- mdf %>%
    dplyr::mutate(ensembl = strip_ensembl_version(as.character(.data[[gene_col]])))
  
  # Connectivity column (auto-detect)
  conn_candidates <- c("kTotal","kWithin","connectivity","kME","intramodularConnectivity")
  conn_col <- conn_candidates[conn_candidates %in% colnames(mdf)]
  if (length(conn_col) == 0) {
    mdf$connectivity <- NA_real_
  } else {
    mdf$connectivity <- suppressWarnings(as.numeric(mdf[[conn_col[1]]]))
  }
  
  # NEW: connectivity z-score (your file has connectivity_zscore)
  z_candidates <- c("connectivity_zscore","connectivity_z_score","connectivity.zscore")
  z_col <- z_candidates[z_candidates %in% colnames(mdf)]
  if (length(z_col) == 0) {
    mdf$connectivity_zscore <- NA_real_
  } else {
    mdf$connectivity_zscore <- suppressWarnings(as.numeric(mdf[[z_col[1]]]))
  }
  
  # Module p-value column (auto-detect)
  pval_candidates <- c("p-val","pval","p_val","p.value","P.Value","PVAL","pvalue","p_value","PValue")
  pval_col <- pval_candidates[pval_candidates %in% colnames(mdf)]
  if (length(pval_col) == 0) {
    mdf$module_pval <- NA_real_
  } else {
    mdf$module_pval <- suppressWarnings(as.numeric(mdf[[pval_col[1]]]))
  }
  
  mdf %>%
    dplyr::select(ensembl, connectivity, connectivity_zscore, module_pval) %>%
    dplyr::distinct(ensembl, .keep_all = TRUE)
}

# Name modules from filenames (adjust colors list if needed)
module_name_from_file <- function(x) {
  nm <- basename(x)
  nm <- sub("\\.csv$", "", nm)
  nm <- sub(".*?(lightcoral|dimgrey|black|snow|silver|firebrick).*", "\\1", nm, ignore.case = TRUE)
  nm
}

# ---------------------------
# Plot saving helper: PNG + SVG
# ---------------------------
save_plot_png_svg <- function(plot_obj, png_file, svg_file,
                              width = 6, height = 4.5, dpi = 300,
                              bg = "white") {
  
  ggsave(
    filename = png_file,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    bg = bg
  )
  
  ggsave(
    filename = svg_file,
    plot = plot_obj,
    device = svglite::svglite,
    width = width,
    height = height,
    bg = bg
  )
}

# ---------------------------
# Locate files
# ---------------------------
deg_files <- list.files(deg_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(deg_files) > 0)

module_files <- list.files(
  wgcna_root,
  pattern = "module.*genes.*\\.csv$|module_genes.*\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(module_files) == 0) {
  module_files <- list.files(wgcna_root, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  message("No obvious module gene CSVs matched the pattern; found ", length(module_files),
          " CSVs total under wgcna_root. You may need to tighten the pattern.")
}
stopifnot(length(module_files) > 0)

modules <- tibble(
  module_file = module_files,
  module = map_chr(module_files, module_name_from_file)
) %>%
  dplyr::distinct(module, module_file)

# ---------------------------
# Preload all module per-gene tables (faster + consistent)
# ---------------------------
module_df_list <- setNames(
  lapply(modules$module_file, read_module_df),
  modules$module
)

# Build ENS->SYMBOL mapping from all module genes + all DEG genes
all_module_ens <- unlist(lapply(module_df_list, function(x) x$ensembl), use.names = FALSE)

all_deg_ens <- deg_files %>%
  purrr::map(function(f) {
    df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    gene_col <- c("gene_id","gene","ensembl","ENSEMBL","X") |> keep(~ .x %in% names(df))
    if (length(gene_col) == 0) return(character())
    strip_ensembl_version(df[[gene_col[1]]])
  }) %>%
  unlist(use.names = FALSE)

ens2sym <- ens2sym_tbl(c(all_module_ens, all_deg_ens))

# ---------------------------
# Main: overlap DEGs with modules + annotate
# ---------------------------
direction <- "either"
p_cut <- 0.05
lfc_cut <- 0.5

overlap_summary <- list()

for (deg_path in deg_files) {
  
  deg_name <- sub("\\.csv$", "", basename(deg_path))
  deg_sig  <- read_deg_sig_df(deg_path, p_cut = p_cut, lfc_cut = lfc_cut, direction = direction)
  
  # if no sig genes, still write summary rows (optional); here we keep going but overlaps will be empty
  for (i in seq_len(nrow(modules))) {
    
    mod <- modules$module[i]
    moddf <- module_df_list[[mod]]
    
    # overlap by ensembl IDs
    inter_ens <- base::intersect(deg_sig$ensembl, moddf$ensembl)
    
    # assemble output table with requested columns
    out_tbl <- tibble(ensembl = inter_ens) %>%
      left_join(ens2sym, by = c("ensembl" = "ENSEMBL")) %>%
      dplyr::rename(symbol = SYMBOL) %>%
      left_join(moddf, by = "ensembl") %>%
      left_join(deg_sig, by = "ensembl") %>%
      dplyr::select(ensembl, symbol, connectivity, connectivity_zscore, module_pval, log2FoldChange, padj) %>%
      dplyr::arrange(dplyr::desc(connectivity_zscore))
    
    out_csv <- file.path(out_overlap_dir, paste0(deg_name, "__", mod, "__overlap.csv"))
    write.csv(out_tbl, out_csv, row.names = FALSE)
    
    overlap_summary[[length(overlap_summary) + 1]] <- tibble(
      deg = deg_name,
      module = mod,
      n_deg_sig = nrow(deg_sig),
      n_module = nrow(moddf),
      n_overlap = length(inter_ens),
      overlap_file = out_csv
    )
  }
}

overlap_summary_df <- bind_rows(overlap_summary) %>%
  arrange(desc(n_overlap))

write.csv(overlap_summary_df,
          file.path(out_overlap_dir, "overlap_summary.csv"),
          row.names = FALSE)

overlap_summary_df %>% head(20)



#venn diagram overlap

library(tidyverse)
library(ggvenn)

# output folder for venns
out_venn_dir <- file.path(output_base, "DEG_vs_Module_venns")
out_venn_png_dir <- file.path(out_venn_dir, "PNG")
out_venn_svg_dir <- file.path(out_venn_dir, "SVG")

dir.create(out_venn_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_venn_png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_venn_svg_dir, recursive = TRUE, showWarnings = FALSE)

# If you didn't already define it above, keep your helper here
make_ggvenn2 <- function(A, B,
                         A_name = "DEG (sig)", B_name = "Module",
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
    set_name_size = 5,
    text_size = 4,
    show_percentage = TRUE,
    show_elements = FALSE
  ) +
    ggtitle(title) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# Choose direction for DEG set: "up", "down", or "either"
direction <- "either"

# Loop and plot
for (deg_path in deg_files) {
  
  deg_name <- sub("\\.csv$", "", basename(deg_path))
  deg_set <- read_deg_sig_df(
    deg_path,
    p_cut = 0.05,
    lfc_cut = 0.5,
    direction = direction
  )$ensembl
  
  # If a DEG comparison yields 0 sig genes, skip (prevents empty plots)
  if (length(deg_set) == 0) {
    message("Skipping Venns for ", deg_name, " (0 significant genes).")
    next
  }
  
  for (i in seq_len(nrow(modules))) {
    
    mod <- modules$module[i]
    mod_set <- module_df_list[[mod]]$ensembl
    
    p <- make_ggvenn2(
      A = deg_set,
      B = mod_set,
      A_name = "DEG (sig)",
      B_name = paste0("Module: ", mod),
      title = paste0(deg_name, " vs ", mod),
      fill = c("#b2182b", "#2166ac")  # red/blue; change if you want
    )
    
    out_png <- file.path(out_venn_png_dir, paste0(deg_name, "__", mod, "__venn.png"))
    out_svg <- file.path(out_venn_svg_dir, paste0(deg_name, "__", mod, "__venn.svg"))
    
    save_plot_png_svg(
      plot_obj = p,
      png_file = out_png,
      svg_file = out_svg,
      width = 6.0,
      height = 4.5,
      dpi = 300
    )
  }
}

message("Saved venn diagrams to: ", out_venn_dir,
        " | PNG: ", out_venn_png_dir,
        " | SVG: ", out_venn_svg_dir)




#venn diagrams for further comparisons
# ============================================================
# Module-overlap vs Module-overlap (Normal vs IPF) for ALL pairs
# - For each DEG pair (normal vs IPF) with same suffix:
#     overlap_normal = sigDEG_normal ∩ moduleGenes
#     overlap_ipf    = sigDEG_ipf    ∩ moduleGenes
#   -> Venn(overlap_normal, overlap_ipf)
#   -> CSV of shared / normal_only / ipf_only with SYMBOL, connectivity, module_pval, log2FC+padj
# ============================================================

if (requireNamespace("conflicted", quietly = TRUE)) {
  library(conflicted)
  conflict_prefer("select", "dplyr")
  conflict_prefer("filter", "dplyr")
  conflict_prefer("rename", "dplyr")
  conflict_prefer("mutate", "dplyr")
}

conflicted::conflicts_prefer(base::intersect)
conflicted::conflicts_prefer(base::setdiff)



strip_ensembl_version <- function(x) sub("\\.\\d+$", "", x)

# ---- 1) Read DEG CSV and keep significant genes + stats ----
read_deg_sig_df <- function(csv_path, p_cut = 0.05, lfc_cut = 0.5,
                            direction = c("up","down","either")) {
  direction <- match.arg(direction)
  
  df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  gene_col <- c("gene_id") |> keep(~ .x %in% names(df))
  if (length(gene_col) == 0) stop("No gene id column found in: ", basename(csv_path))
  gene_col <- gene_col[1]
  
  if (!("padj" %in% names(df))) stop("padj not found in: ", basename(csv_path))
  if (!("log2FoldChange" %in% names(df))) stop("log2FoldChange not found in: ", basename(csv_path))
  
  df <- df %>%
    mutate(ensembl = strip_ensembl_version(.data[[gene_col]])) %>%
    filter(!is.na(padj), padj < p_cut, !is.na(log2FoldChange))
  
  if (direction == "up") {
    df <- df %>% filter(log2FoldChange >= lfc_cut)
  } else if (direction == "down") {
    df <- df %>% filter(log2FoldChange <= -lfc_cut)
  } else {
    df <- df %>% filter(abs(log2FoldChange) >= lfc_cut)
  }
  
  keep_cols <- intersect(c("ensembl","log2FoldChange","padj","stat","lfcSE","baseMean"), names(df))
  
  df %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    dplyr::distinct(ensembl, .keep_all = TRUE)
}

# ---- 2) Ensembl -> SYMBOL ----
ens2sym_tbl <- function(keys) {
  keys <- strip_ensembl_version(keys)
  keys <- unique(keys)
  keys <- keys[!is.na(keys) & keys != ""]
  
  # keep only ENSG IDs (prevents garbage keys from crashing AnnotationDbi)
  keys <- keys[grepl("^ENSG\\d+$", keys)]
  
  if (length(keys) == 0) {
    return(tibble(ENSEMBL=character(), SYMBOL=character()))
  }
  
  # Keep only keys that org.Hs.eg.db recognizes (prevents hard error)
  valid <- keys %in% keys(org.Hs.eg.db, keytype = "ENSEMBL")
  keys_valid <- keys[valid]
  
  if (length(keys_valid) == 0) {
    return(tibble(ENSEMBL=character(), SYMBOL=character()))
  }
  
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys     = keys_valid,
    keytype  = "ENSEMBL",
    columns  = c("ENSEMBL","SYMBOL")
  ) %>%
    as_tibble() %>%
    dplyr::distinct(ENSEMBL, .keep_all = TRUE)
}


# ---- 3) Load module genes + connectivity + module_pval ----
read_module_df <- function(module_csv_path) {
  
  mdf <- read.csv(module_csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Fix bad/blank/NA column names
  cn <- colnames(mdf)
  cn[is.na(cn) | cn == ""] <- paste0("X", which(is.na(cn) | cn == ""))
  colnames(mdf) <- make.unique(cn)
  
  # Pick gene id column
  gene_col <- c("ensembl","ENSEMBL","gene_id","X","var_names","id") |> keep(~ .x %in% colnames(mdf))
  if (length(gene_col) == 0) gene_col <- colnames(mdf)[1] else gene_col <- gene_col[1]
  
  mdf <- mdf %>%
    dplyr::mutate(ensembl = strip_ensembl_version(.data[[gene_col]]))
  
  # Connectivity (if present)
  conn_col <- c("kTotal","kWithin","connectivity","kME","intramodularConnectivity") |> keep(~ .x %in% colnames(mdf))
  if (length(conn_col) == 0) {
    mdf$connectivity <- NA_real_
  } else {
    mdf$connectivity <- suppressWarnings(as.numeric(mdf[[conn_col[1]]]))
  }
  
  # Module p-value (auto-detect)
  pval_candidates <- c("p-val","pval","p_val","p.value","P.Value","PVAL","pvalue","p_value","PValue")
  pval_col <- pval_candidates[pval_candidates %in% colnames(mdf)]
  if (length(pval_col) == 0) {
    mdf$module_pval <- NA_real_
  } else {
    mdf$module_pval <- suppressWarnings(as.numeric(mdf[[pval_col[1]]]))
  }
  
  mdf %>%
    dplyr::select(ensembl, connectivity, module_pval) %>%
    dplyr::distinct(ensembl, .keep_all = TRUE)
}

# ---- 4) Make 2-set Venn ----
venn2 <- function(A, B, A_name, B_name, title="") {
  ggvenn(setNames(list(unique(A), unique(B)), c(A_name, B_name)),
         fill_color = c("#b2182b","#2166ac"),
         stroke_size = 0.6,
         set_name_size = 5,
         text_size = 4,
         show_percentage = TRUE,
         show_elements = FALSE) +
    ggtitle(title) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face="bold", size=13, hjust=0.5))
}

# ============================================================
# INPUTS YOU ALREADY HAVE FROM STEP 1/2
#   deg_dir, wgcna_root, output_base, modules, deg_files
# If deg_files/modules not in memory, uncomment the locate blocks
# ============================================================

# If you need: locate deg_files
# deg_files <- list.files(deg_dir, pattern="\\.csv$", full.names=TRUE)

# If you need: locate module files and modules tibble
# module_files <- list.files(wgcna_root, pattern="\\.csv$", full.names=TRUE, recursive=TRUE)
# module_name_from_file <- function(x) {
#   nm <- basename(x)
#   nm <- sub("\\.csv$", "", nm)
#   nm <- sub(".*?(lightcoral|dimgrey|black|snow|silver|firebrick).*", "\\1", nm, ignore.case = TRUE)
#   nm
# }
# modules <- tibble(module_file = module_files, module = map_chr(module_files, module_name_from_file)) %>%
#   distinct(module, module_file)

# thresholds
p_cut <- 0.05
lfc_cut <- 0.5
direction <- "either"

# output folder
out_pair_dir <- file.path(output_base, "ModuleOverlapPairs_ALL")
out_pair_png_dir <- file.path(out_pair_dir, "PNG")
out_pair_svg_dir <- file.path(out_pair_dir, "SVG")

dir.create(out_pair_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_pair_png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_pair_svg_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Build DEG pairs: match DEG_normal_* with DEG_IPF_* by suffix
# e.g. DEG_normal_IR_IAV_vs_normal_IR_mock  <->  DEG_IPF_IR_IAV_vs_IPF_IR_mock
# ============================================================

deg_tbl <- tibble(path = deg_files, file = basename(deg_files)) %>%
  mutate(stem = sub("\\.csv$", "", file)) %>%
  mutate(group = case_when(
    grepl("^DEG_normal_", stem, ignore.case = TRUE) ~ "normal",
    grepl("^DEG_IPF_", stem, ignore.case = TRUE)    ~ "ipf",
    TRUE ~ NA_character_
  )) %>%
  dplyr::filter(!is.na(group)) %>%
  # remove leading DEG_normal_ / DEG_IPF_
  mutate(suffix_raw = sub("^DEG_(normal|IPF)_", "", stem, ignore.case = TRUE)) %>%
  # normalize internal normal/IPF tokens so they match across groups
  mutate(pair_key = suffix_raw %>%
           tolower() %>%
           gsub("ipf", "group", ., fixed = TRUE) %>%
           gsub("normal", "group", ., fixed = TRUE)
  )

pairs_tbl <- deg_tbl %>%
  dplyr::select(group, pair_key, path, stem) %>%
  tidyr::pivot_wider(names_from = group, values_from = c(path, stem)) %>%
  dplyr::filter(!is.na(path_normal) & !is.na(path_ipf))

# sanity check:
nrow(pairs_tbl)
pairs_tbl %>% dplyr::select(pair_key, stem_normal, stem_ipf)


# ============================================================
# Cache DEG sig tables so we don't re-read in inner loops
# ============================================================

deg_cache <- new.env(parent = emptyenv())

get_deg_cached <- function(path) {
  key <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (exists(key, envir = deg_cache, inherits = FALSE)) {
    return(get(key, envir = deg_cache, inherits = FALSE))
  }
  df <- read_deg_sig_df(path, p_cut = p_cut, lfc_cut = lfc_cut, direction = direction)
  assign(key, df, envir = deg_cache)
  df
}

# ============================================================
# Loop: all pairs x all modules
# ============================================================

summary_list <- list()

for (pair_i in seq_len(nrow(pairs_tbl))) {
  
  normal_path <- pairs_tbl$path_normal[pair_i]
  ipf_path    <- pairs_tbl$path_ipf[pair_i]
  
  normal_name <- pairs_tbl$stem_normal[pair_i]
  ipf_name    <- pairs_tbl$stem_ipf[pair_i]
  suffix <- pairs_tbl$pair_key[pair_i]
  
  
  degN <- get_deg_cached(normal_path)
  degI <- get_deg_cached(ipf_path)
  
  # skip empty DEG sets
  if (nrow(degN) == 0 && nrow(degI) == 0) next
  
  for (m in seq_len(nrow(modules))) {
    
    module_color <- modules$module[m]
    module_csv   <- modules$module_file[m]
    
    moddf <- read_module_df(module_csv)
    module_set <- moddf$ensembl
    
    normal_overlap <- intersect(degN$ensembl, module_set)
    ipf_overlap    <- intersect(degI$ensembl, module_set)
    
    # skip if both empty
    if (length(normal_overlap) == 0 && length(ipf_overlap) == 0) next
    
    # venn
    p <- venn2(
      normal_overlap, ipf_overlap,
      A_name = "Normal (sig ∩ module)",
      B_name = "IPF (sig ∩ module)",
      title  = paste0("Module ", module_color, ": ", suffix)
    )
    
    venn_path_png <- file.path(out_pair_png_dir, paste0("module_", module_color, "__", suffix, "__overlapVenn.png"))
    venn_path_svg <- file.path(out_pair_svg_dir, paste0("module_", module_color, "__", suffix, "__overlapVenn.svg"))
    
    save_plot_png_svg(
      plot_obj = p,
      png_file = venn_path_png,
      svg_file = venn_path_svg,
      width = 6.5,
      height = 4.8,
      dpi = 300
    )
    
    # sets
    shared   <- intersect(normal_overlap, ipf_overlap)
    normal_u <- setdiff(normal_overlap, ipf_overlap)
    ipf_u    <- setdiff(ipf_overlap, normal_overlap)
    
    # annotate
    ens_all <- unique(c(shared, normal_u, ipf_u))
    map_tbl <- ens2sym_tbl(ens_all)
    
    annotate_set <- function(ens_vec, label) {
      tibble(ensembl = ens_vec) %>%
        left_join(map_tbl, by = c("ensembl" = "ENSEMBL")) %>%
        dplyr::rename(symbol = SYMBOL) %>%
        left_join(moddf, by = "ensembl") %>%
        left_join(degN, by = "ensembl", suffix = c("", "_normal")) %>%
        left_join(degI, by = "ensembl", suffix = c("", "_ipf")) %>%
        dplyr::mutate(set = label) %>%
        dplyr::select(
          set, ensembl, symbol,
          connectivity, module_pval,
          log2FoldChange, padj,
          log2FoldChange_ipf, padj_ipf
        ) %>%
        dplyr::arrange(dplyr::desc(connectivity))
    }
    
    df_out <- bind_rows(
      annotate_set(shared,   "shared"),
      annotate_set(normal_u, "normal_only"),
      annotate_set(ipf_u,    "ipf_only")
    )
    
    out_csv <- file.path(out_pair_dir, paste0("module_", module_color, "__", suffix, "__overlap_sets.csv"))
    write.csv(df_out, out_csv, row.names = FALSE)
    
    summary_list[[length(summary_list) + 1]] <- tibble(
      module = module_color,
      suffix = suffix,
      normal_deg = normal_name,
      ipf_deg = ipf_name,
      n_normal_overlap = length(normal_overlap),
      n_ipf_overlap = length(ipf_overlap),
      n_shared = length(shared),
      venn_png = venn_path_png,
      overlap_csv = out_csv
    )
  }
}

summary_df <- bind_rows(summary_list) %>%
  arrange(desc(n_shared), desc(n_normal_overlap + n_ipf_overlap))

write.csv(summary_df, file.path(out_pair_dir, "ModuleOverlapPairs_summary.csv"), row.names = FALSE)

cat("Done. Outputs in:\n", out_pair_dir, "\n",
    "Summary:\n", file.path(out_pair_dir, "ModuleOverlapPairs_summary.csv"), "\n")





#GO enrichment for modules
#helper
geneRatio_to_numeric <- function(x) {
  sapply(strsplit(as.character(x), "/"),
         function(z) as.numeric(z[1]) / as.numeric(z[2]))
}

palette_by_direction <- function(direction = c("UP", "DOWN")) {
  direction <- match.arg(direction)
  if (direction == "UP") {
    list(low = "#FEE0D2", high = "#99000D")
  } else {
    list(low = "#DEEBF7", high = "#08519C")
  }
}

bubble_enrich <- function(enrich_obj, title, direction,
                          showCategory = 20, wrap_width = 40) {
  df <- as.data.frame(enrich_obj)
  if (nrow(df) == 0) return(NULL)
  
  df <- df[order(df$p.adjust), ]
  df <- head(df, showCategory)
  df <- df[order(df$Count, decreasing = TRUE), ]
  
  df$GeneRatio_num <- geneRatio_to_numeric(df$GeneRatio)
  df$Description_wrapped <- str_wrap(df$Description, wrap_width)
  
  df$Description_wrapped <- factor(
    df$Description_wrapped,
    levels = rev(df$Description_wrapped)
  )
  
  pal <- palette_by_direction(direction)
  
  ggplot(df, aes(GeneRatio_num, Description_wrapped)) +
    geom_point(aes(size = Count, color = p.adjust), alpha = 0.9) +
    scale_color_gradient(
      name = "p.adjust",
      trans = "log10",
      low = pal$low,
      high = pal$high
    ) +
    scale_size_continuous(name = "Count") +
    labs(x = "GeneRatio", y = NULL, title = title) +
    theme_bw() +
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

sym_to_entrez <- function(symbols) {
  df <- bitr(symbols,
             fromType = "SYMBOL",
             toType   = "ENTREZID",
             OrgDb    = org.Hs.eg.db)
  unique(df$ENTREZID)
}

library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(stringr)
library(ggplot2)

# ---------------------------
# Paths
# ---------------------------
# WGCNA output root for the new batch-corrected WGCNA run
wgcna_root <- "Jena_count_DEseq2_batchCorrected_nonzero_50mod_0.4thresh_robustZ_diffNames"
# Update this path if your new batch-corrected WGCNA run has a different folder name.
stopifnot(dir.exists(wgcna_root))

go_out_dir <- file.path(base_dir, "results_batch", "GO_by_module")
go_png_dir <- file.path(go_out_dir, "PNG")
go_svg_dir <- file.path(go_out_dir, "SVG")
go_png_dir_top10 <- file.path(go_out_dir, "PNG_top10")
go_svg_dir_top10 <- file.path(go_out_dir, "SVG_top10")

dir.create(go_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(go_png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(go_svg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(go_png_dir_top10, recursive = TRUE, showWarnings = FALSE)
dir.create(go_svg_dir_top10, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Helpers
# ---------------------------
strip_ensembl_version <- function(x) sub("\\.\\d+$", "", x)

# robustly read module genes from CSV
get_module_genes <- function(csv_path) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  cand <- c("ensembl","ENSEMBL","gene_id","X","var_names","id") |> keep(~ .x %in% names(df))
  genes <- if (length(cand) > 0) df[[cand[1]]] else df[[1]]
  genes <- strip_ensembl_version(as.character(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  unique(genes)
}

# map Ensembl -> Entrez (drop NAs)
ens_to_entrez <- function(ens) {
  ens <- unique(strip_ensembl_version(ens))
  ens <- ens[grepl("^ENSG\\d+$", ens)]
  if (length(ens) == 0) return(character())
  suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = ens,
      keytype = "ENSEMBL",
      columns = c("ENSEMBL","ENTREZID")
    )
  ) %>%
    as_tibble() %>%
    filter(!is.na(ENTREZID) & ENTREZID != "") %>%
    distinct(ENTREZID) %>%
    pull(ENTREZID)
}

# identify module name from file name
module_name_from_file <- function(x) {
  nm <- basename(x)
  nm <- sub("\\.csv$", "", nm)
  nm <- sub(".*?(lightcoral|dimgrey|black|snow|silver|firebrick).*", "\\1", nm, ignore.case = TRUE)
  tolower(nm)
}

# ---------------------------
# Find module gene list CSVs
# ---------------------------
module_files <- list.files(
  wgcna_root,
  pattern = "module.*genes.*\\.csv$|module_genes.*\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(module_files) == 0) {
  module_files <- list.files(wgcna_root, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  message("No obvious module gene CSVs matched the pattern; found ", length(module_files),
          " CSVs total under wgcna_root. You may need to tighten the pattern.")
}
stopifnot(length(module_files) > 0)

modules_tbl <- tibble(
  module_file = module_files,
  module = map_chr(module_files, module_name_from_file)
) %>%
  distinct(module, module_file)

# ---------------------------
# Run GO per module
# ---------------------------
ontologies <- c("BP","CC","MF")

summary_rows <- list()

for (i in seq_len(nrow(modules_tbl))) {
  
  mod <- modules_tbl$module[i]
  f   <- modules_tbl$module_file[i]
  
  ens <- get_module_genes(f)
  entrez <- ens_to_entrez(ens)
  
  if (length(entrez) < 10) {
    message("Skipping module ", mod, ": too few Entrez IDs (", length(entrez), ")")
    next
  }
  
  for (ont in ontologies) {
    
    ego <- enrichGO(
      gene          = entrez,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
    
    # save table
    res_df <- as.data.frame(ego)
    out_csv <- file.path(go_out_dir, paste0("GO_", ont, "_module_", mod, ".csv"))
    write.csv(res_df, out_csv, row.names = FALSE)
    
    # bubble plot with your helper
    # modules don't have inherent UP/DOWN; choose one palette (UP=red) unless you want to define per module.
    p <- bubble_enrich(
      enrich_obj = ego,
      title      = paste0("GO-", ont, " enrichment: module ", mod),
      direction  = "UP",
      showCategory = 20
    )
    
    if (!is.null(p)) {
      out_png <- file.path(go_png_dir, paste0("GO_", ont, "_module_", mod, "_bubble.png"))
      out_svg <- file.path(go_svg_dir, paste0("GO_", ont, "_module_", mod, "_bubble.svg"))
      
      n_show <- 20
      height_in <- max(6, 0.35 * n_show + 2)
      
      save_plot_png_svg(
        plot_obj = p,
        png_file = out_png,
        svg_file = out_svg,
        width = 9,
        height = height_in,
        dpi = 300
      )
    }
    
    p_top10 <- bubble_enrich(
      enrich_obj = ego,
      title      = paste0("GO-", ont, " enrichment: module ", mod, " (top 10)"),
      direction  = "UP",
      showCategory = 10
    )
    
    if (!is.null(p_top10)) {
      out_png_top10 <- file.path(go_png_dir_top10, paste0("GO_", ont, "_module_", mod, "_bubble_top10.png"))
      out_svg_top10 <- file.path(go_svg_dir_top10, paste0("GO_", ont, "_module_", mod, "_bubble_top10.svg"))
      
      n_show_top10 <- 10
      height_in_top10 <- max(6, 0.35 * n_show_top10 + 2)
      
      save_plot_png_svg(
        plot_obj = p_top10,
        png_file = out_png_top10,
        svg_file = out_svg_top10,
        width = 9,
        height = height_in_top10,
        dpi = 300
      )
    }
    
    summary_rows[[length(summary_rows) + 1]] <- tibble(
      module = mod,
      ontology = ont,
      n_input_entrez = length(entrez),
      n_terms = nrow(res_df),
      results_csv = out_csv
    )
  }
}

go_summary <- bind_rows(summary_rows) %>% arrange(module, ontology)
write.csv(go_summary, file.path(go_out_dir, "GO_module_summary.csv"), row.names = FALSE)

go_summary

