#!/usr/bin/env Rscript

# IBD_multiome_analysis_GSEA_RNA_Status.R
#
# Gene set enrichment analysis (GSEA) of pseudobulk RNA differential-expression
# statistics comparing endoscopic healing status (I vs. NI), resolved by
# cell type, disease and tissue contexts.
#
# Associated analysis: Yang, Wayman, Denson, et al. (2026), Figure 3A
#
# NOTE:
#   All paths below are project-relative and should be updated to match the
#   local directory structure before running this script.

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(42)

suppressPackageStartupMessages({
  library(ggplot2)
  library(magrittr)
  library(sva)
  library(fgsea)
  library(reshape2)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(scales)
  library(matrixStats)
  library(dplyr)
  library(patchwork)
  library(tibble)
})

# ==============================================================================
# Configuration
# ==============================================================================

# Project directories ----------------------------------------------------------
project_dir <- "."
script_dir <- file.path(project_dir, "scripts")
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")

# Output directory
base_dir_out <- file.path(results_dir, "GSEA_RNA_Status")

dir.create(base_dir_out, showWarnings = FALSE, recursive = TRUE)

# DESeq2 results ---------------------------------------------------------------
name_res <- "1Factor_CellType_Downsample"
file_data <- file.path(
  results_dir, "res_DESeq2_IBD_multiome_RNA_Disease_Tisssue_Status_raw_unfiltered.txt"
)

# Ligand-receptor pair information --------------------------------------------
# Immune Dictionary
file_lr <- file.path(
  data_dir, "gene_sets", "ImmuneDictionary", "ImmDict_LR_pairs_by_cellchatDB.txt"
)
df_lr <- read.table(file_lr, header = TRUE, sep = "\t")

# CytoSig alternative:
# file_lr <- file.path(
#   data_dir, "gene_sets", "CytoSig", "CytoSig_LR_pairs_by_cellchatDB.txt"
# )
# df_lr <- read.table(file_lr, header = TRUE, sep = "\t")

# Gene sets --------------------------------------------------------------------
# Immune Dictionary: ImmDict_max200
name_gene_set <- "ImmDict_max200"
file_gene_set <- file.path(
  data_dir, "gene_sets", "ImmuneDictionary", "gene_set_combined_ImmuneDict_max200.rds"
)

# PROGENy + CD64 alternative:
# name_gene_set <- "PROGENy_CD64_max200"
# file_gene_set <- c(
#   file.path(data_dir, "gene_sets", "PROGENy_CD64", "gene_set_combined_PROGENy_max200.rds"),
#   file.path(data_dir, "gene_sets", "PROGENy_CD64", "gene_set_combined_Gao2020CD64_max200.rds")
# )

# CytoSig alternative:
# name_gene_set <- "CytoSig"
# file_gene_set <- file.path(
#   data_dir, "gene_sets", "CytoSig", "CytoSig_top200_gene.rds"
# )

# Background genes -------------------------------------------------------------
# Set to NULL when no restricted background is required (e.g., PROGENy, CD64,
# or CytoSig).
file_bkgrd_keep <- file.path(
  data_dir, "gene_sets", "ImmuneDictionary", "gene_set_ImmuneDict_Mm2Hs_bkgrd.txt"
)
# file_bkgrd_keep <- NULL

# Analysis parameters ----------------------------------------------------------
file_save <- "IBD_multiome"
select_disease_tissue <- "CD_TI"  # Options: CD_R, UC_R, CD_TI

order_celltype <- c(
  "CD4_T_cell", "CD4_T_Eff", "CD8_T_cell", "Naive_T_cell", "gdT_cell",
  "NK_cell", "B_cell", "Plasma_cell", "Macrophage", "DC", "Cycling_Immune",
  "Fibroblast", "BEST4_Colonocyte", "Colonocyte", "Immature_Colonocyte",
  "Enterocyte", "Immature_Enterocyte", "Enteroendocrine", "Endothelial",
  "Goblet", "Immature_Goblet", "TA", "Cycling_TA", "Stem", "Tuft",
  "Paneth", "M_cell"
)

# Remove mitochondrial and ribosomal genes for Immune Dictionary analyses.
# Use NULL for PROGENy, CD64, or CytoSig analyses.
del_gene_pattern <- "^MT-|^RPS|^RPL"
# del_gene_pattern <- NULL

n_cores <- 12
min_gene <- 150 # minimum gene set size for analysis

# Enrichment method: "fgsea" or "RNAEnrich"
type_method <- "fgsea"

# Ranking metric for fgsea:
#   "LFC"      = log2 fold-change
#   "SignPadj" = sign(LFC) * -log10(adjusted p-value)
#   "SignPraw" = sign(LFC) * -log10(raw p-value)
#   "statWald" = Wald statistic (LFC / SE[LFC])
type_ranking <- "statWald"

fdr <- 0.1

# ==============================================================================
# Set output paths and load optional methods
# ==============================================================================

if (type_method == "fgsea") {
  dir_out <- file.path(
    base_dir_out,
    paste(name_res, type_method, type_ranking, sep = "_"),
    select_disease_tissue,
    name_gene_set
  )
} else if (type_method == "RNAEnrich") {
  dir_out <- file.path(
    base_dir_out,
    paste(name_res, type_method, sep = "_"),
    select_disease_tissue,
    name_gene_set
  )
}

dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

if (type_method == "RNAEnrich") {
  source(file.path(script_dir, "RNA-Enrich.r"))
}

# ==============================================================================
# Prepare background genes
# ==============================================================================

if (!is.null(file_bkgrd_keep)) {
  gene_bkgrd_keep <- readLines(file_bkgrd_keep)
}

if (!is.null(del_gene_pattern)) {
  gene_bkgrd_keep <- gene_bkgrd_keep[!grepl(del_gene_pattern, gene_bkgrd_keep)]
}

# ==============================================================================
# Load differential-expression statistics
# ==============================================================================

message("Loading differential-expression statistics...")

results <- read.delim(file_data, sep = "\t", header = TRUE)
results$Context <- paste(results$CellType, results$DiseaseTissue, sep = "_")

# Orient fold changes consistently as I vs. NI.
flip_idx <- results$Ident1 == "NI" & results$Ident2 == "I"
results[flip_idx, "Log2FC"] <- -results[flip_idx, "Log2FC"]

# Define context order from the specified cell-type and disease/tissue order.
order_clustersubtype <- NULL
for (ix in order_celltype) {
  for (jx in select_disease_tissue) {
    order_clustersubtype <- c(order_clustersubtype, paste(ix, jx, sep = "_"))
  }
}
order_clustersubtype <- intersect(order_clustersubtype, unique(results$Context))

# ==============================================================================
# Load gene sets
# ==============================================================================

message("Loading gene sets...")

gene_sets <- list()
for (ix in file_gene_set) {
  gene_sets <- c(gene_sets, readRDS(ix))
}

gene_sets <- gene_sets[lengths(gene_sets) >= min_gene]

file_out <- file.path(dir_out, paste0(name_gene_set, "_unique_cytokine_names.txt"))
writeLines(names(gene_sets), file_out)

# ==============================================================================
# Run enrichment analysis by cell type / disease / tissue context
# ==============================================================================

message("Running enrichment analysis...")

gsea_all <- list()
sig_sets <- NULL

for (ix in order_clustersubtype) {
  message("Context: ", ix)

  # Define the nominally expressed background for this context.
  curr_nominally_expressed_gene <- as.character(unique(subset(results, Context == ix)$Gene))

  if (exists("gene_bkgrd_keep")) {
    curr_nominally_expressed_gene <- intersect(curr_nominally_expressed_gene, gene_bkgrd_keep)
  }

  # Restrict DE statistics to the context-specific background.
  curr_data <- subset(results, Context == ix)
  curr_data <- subset(curr_data, Gene %in% curr_nominally_expressed_gene)
  curr_data <- unique(curr_data[, c(
    "Context", "Gene", "baseMean", "Log2FC", "padj", "pval", "stat"
  )])

  if (type_method == "fgsea") {
    # Rank genes using the selected DE statistic.
    if (type_ranking == "LFC") {
      gene_rank <- setNames(curr_data$Log2FC, curr_data$Gene)
    } else if (type_ranking == "SignPadj") {
      gene_rank <- setNames(sign(curr_data$Log2FC) * -log10(curr_data$padj), curr_data$Gene)
    } else if (type_ranking == "SignPraw") {
      gene_rank <- setNames(sign(curr_data$Log2FC) * -log10(curr_data$pval), curr_data$Gene)
    } else if (type_ranking == "statWald") {
      gene_rank <- setNames(curr_data$stat, curr_data$Gene)
    }

    gene_rank <- sort(gene_rank, decreasing = TRUE)

    # Randomly resolve ties while retaining the original workflow.
    order_gene_rank <- names(rank(gene_rank, ties.method = "random"))
    gene_rank <- gene_rank[order_gene_rank]

    curr_gsea <- fgsea(
      pathways = gene_sets,
      stats = gene_rank,
      scoreType = "std",
      minSize = 1,
      nproc = n_cores
    )
  } else if (type_method == "RNAEnrich") {
    # Prepare differential-expression results for RNA-Enrich.
    results_df <- curr_data %>%
      as.data.frame() %>%
      filter(Gene %in% curr_data$Gene, !is.na(pval), !is.na(baseMean))

    # Restrict custom gene sets to the current background.
    gene_sets_bg <- lapply(gene_sets, function(x) unique(x[x %in% curr_data$Gene]))
    gene_sets_bg <- gene_sets_bg[lengths(gene_sets_bg) >= 10]

    # Convert the list of gene sets to RNA-Enrich long format.
    conceptList <- rep(names(gene_sets_bg), lengths(gene_sets_bg))
    nullsetList <- unlist(gene_sets_bg, use.names = FALSE)

    curr_gsea <- rna_enrich(
      sigvals = results_df$pval,
      geneids = results_df$Gene,
      avg_readcount = results_df$baseMean,
      species = "sce",
      database = "custom",
      conceptList = conceptList,
      nullsetList = nullsetList,
      direction = sign(results_df$Log2FC),
      min.g = 5,
      max.g = NA,
      sig.cutoff = fdr,
      read_lim = 0,
      plot_file = file.path(dir_out, "RNA-Enrich_plot.jpg"),
      results_file = file.path(dir_out, "RNA-Enrich_custom_results.txt")
    )

    curr_gsea <- curr_gsea[, -1]
    curr_gsea <- curr_gsea[, -2]
    colnames(curr_gsea)[colnames(curr_gsea) == "Concept.name"] <- "pathway"
    colnames(curr_gsea)[colnames(curr_gsea) == "p.value"] <- "pval"
    colnames(curr_gsea)[colnames(curr_gsea) == "FDR"] <- "padj"
  }

  curr_sig <- curr_gsea$pathway[curr_gsea$padj <= fdr]
  gsea_all[[ix]] <- curr_gsea
  sig_sets <- c(sig_sets, curr_sig)
}

sig_sets <- unique(sig_sets)

# ==============================================================================
# Save results
# ==============================================================================

message("Saving GSEA results...")

# Save all GSEA results.
file_out <- file.path(
  dir_out, paste0("gsea_", file_save, "_", name_gene_set, "_I_vs_NI_Disease_Tissue_All.rds")
)
saveRDS(gsea_all, file_out)

# Save the union of significant gene sets across contexts.
file_out <- file.path(dir_out, paste0("sig_sets", file_save, "_", name_gene_set, "_All.txt"))
writeLines(sig_sets, file_out)

message("Done.")
