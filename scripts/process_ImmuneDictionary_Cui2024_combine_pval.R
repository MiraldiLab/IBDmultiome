#!/usr/bin/env Rscript

# process_ImmuneDictionary_Cui2024_combine_pval.R
#
# Process cytokine-response signatures from the Immune Dictionary
# (Cui et al., Nature, 2024).
#
# Differential-expression results are combined across immune cell populations
# using an unweighted Stouffer Z-score meta-analysis. Positive-response gene
# sets are then defined for each cytokine and optionally mapped from mouse to
# human orthologs.
#
# NOTE:
#   All paths are project-relative and should be updated to match the local
#   directory structure before running this script.

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(42)

suppressPackageStartupMessages({
  library(ggplot2)
})

# ==============================================================================
# Configuration
# ==============================================================================

project_dir <- "."
data_dir <- file.path(project_dir, "data")
dir_out <- file.path(project_dir, "processed", "wt_pval")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Mouse-human ortholog files
file_orth_path <- file.path(data_dir, "meta_cytokine_alternative.txt")
file_orth_gene <- file.path(data_dir, "orthologs_hg38_2_mm10_07032024.txt")

# Analysis parameters
file_save <- "ImmuneDict"
species <- "Hs"

celltypes <- c(
  "B_cell", "cDC1", "cDC2", "eTAC", "ILC", "Langerhans", "Macrophage",
  "Mast_cell", "MigDC", "Monocyte", "Neutrophil", "NK_cell", "pDC",
  "T_cell_CD4", "T_cell_CD8", "T_cell_gd", "Treg"
)

min_npop_up <- 1   # Minimum number of populations in which a gene is upregulated
min_set_size <- 5  # Minimum number of genes required to retain a gene set
max_set_size <- 400

# ==============================================================================
# Load and combine differential-expression results
# ==============================================================================

res_deg <- NULL

for (ix in celltypes) {
  curr_file <- file.path(
    data_dir, "gene_sets", "ImmuneDictionary",
    paste0("tbl_ImmuneDict_cytokine_signature_", ix, ".txt")
  )

  curr_res <- read.delim(curr_file, header = TRUE, check.names = FALSE)
  curr_res <- curr_res[, c("Celltype_Str", "Cytokine_Str", "Gene", "Avg_log2FC", "FDR")]
  colnames(curr_res) <- c("CellType", "Pathway", "Gene", "Avg_log2FC", "Padj")

  res_deg <- rbind(res_deg, curr_res)
}

# ==============================================================================
# Map mouse genes to human orthologs
# ==============================================================================

if (species == "Hs") {
  
  # Optional pathway-name ortholog mapping retained from the original workflow.
  # orth_path <- read.delim(file_orth_path, header = TRUE, sep = "\t", check.names = FALSE)
  # orth_path <- orth_path[, c("Cytokine_Str", "Human gene symbol")]
  # colnames(orth_path) <- c("Pathway", "Pathway_Hs")
  # res_deg_hs <- merge(res_deg, orth_path, by = "Pathway")

  orth_gene <- read.delim(file_orth_gene, header = TRUE, sep = "\t", check.names = FALSE)
  colnames(orth_gene) <- c("Gene", "Gene_Hs")

  res_deg_hs <- merge(res_deg, orth_gene, by = "Gene")
  res_deg_hs$Gene <- as.character(res_deg_hs$Gene_Hs)
  res_deg_hs <- res_deg_hs[, -ncol(res_deg_hs)]

  res_deg <- res_deg_hs
  species_name <- "Mm2Hs"

  # Background genes available after mouse-to-human ortholog mapping
  gene_bkgrd <- sort(unique(orth_gene$Gene_Hs))
  
} else {
  species_name <- "Mm"
}

dir_out <- file.path(dir_out, species_name)
dir.create(dir_out, showWarnings = FALSE)

# ==============================================================================
# Build cytokine-response gene sets
# ==============================================================================

pathways <- sort(unique(res_deg$Pathway))
gene_sets <- list()
df_gene_sets <- NULL

for (ix in pathways) {
  curr_res <- subset(res_deg, Pathway == ix)

  # Convert adjusted p-values to signed Z-scores.
  curr_res$Padj <- pmin(pmax(curr_res$Padj, 1e-15), 1 - 1e-16)
  curr_res$Z <- sign(curr_res$Avg_log2FC) * qnorm(1 - curr_res$Padj / 2)

  # Split observations by gene.
  gene_split <- split(seq_len(nrow(curr_res)), curr_res$Gene)

  # Number of contributing cell types per gene.
  k_vec <- sapply(gene_split, length)

  # Unweighted Stouffer Z-score meta-analysis.
  Z_meta_vec <- sapply(gene_split, function(idx) {
    z <- curr_res$Z[idx]
    sum(z) / sqrt(length(z))
  })
  genes <- names(Z_meta_vec)

  # Two-sided meta-analysis p-values.
  p_meta_vec <- 2 * pnorm(-abs(Z_meta_vec))

  # Number and fraction of cell types with positive log2 fold-change.
  n_up_vec <- sapply(gene_split, function(idx) sum(curr_res$Avg_log2FC[idx] > 0))
  frac_up_vec <- sapply(gene_split, function(idx) mean(curr_res$Avg_log2FC[idx] > 0))

  meta_unweighted <- data.frame(
    Pathway = ix,
    gene = genes,
    k = k_vec[genes],
    Z_meta = signif(Z_meta_vec[genes], digits = 4),
    p_meta = signif(p_meta_vec[genes], digits = 4),
    n_up = n_up_vec[genes],
    frac_up = frac_up_vec[genes]
  )
  meta_unweighted <- meta_unweighted[order(meta_unweighted$p_meta), ]

  # Retain positively regulated genes observed in at least min_npop_up populations.
  signature_unweighted <- meta_unweighted[
    meta_unweighted$Z_meta > 0 & meta_unweighted$n_up >= min_npop_up,
  ]

  if (nrow(signature_unweighted) >= min_set_size) {
    # if (nrow(signature_unweighted) >= 500) stop()
    if (max(signature_unweighted$Z_meta) > 100) stop()

    pval_cutoff <- signature_unweighted[
      min(nrow(signature_unweighted), max_set_size), "p_meta"
    ]
    signature_filtered <- subset(signature_unweighted, p_meta <= pval_cutoff)
    curr_gene_set <- as.character(signature_filtered$gene)

    gene_sets[[ix]] <- curr_gene_set

    df_gene_set <- data.frame(
      Path1 = paste(file_save, "Up", ix, sep = "_"),
      Path2 = paste(file_save, "Up", ix, sep = "_"),
      GeneSet = paste(curr_gene_set, collapse = "|")
    )
    df_gene_sets <- rbind(df_gene_sets, df_gene_set)
  }
}

# ==============================================================================
# Save gene sets
# ==============================================================================

file_out <- file.path(
  dir_out,
  paste0(
    "gene_set_", file_save, "_minPop", min_npop_up, "_maxSize", max_set_size,
    "_", species_name, "_3col.txt"
  )
)
write.table(df_gene_sets, file_out, quote = FALSE, sep = " ",
            col.names = FALSE, row.names = FALSE)

file_out <- file.path(
  dir_out,
  paste0(
    "gene_set_", file_save, "_minPop", min_npop_up, "_maxSize", max_set_size,
    "_", species_name, ".rds"
  )
)
saveRDS(gene_sets, file_out)

# Save background genes
file_out <- file.path(
  dir_out, paste0("gene_set_", file_save, "_", species_name, "_bkgrd.txt")
)
writeLines(gene_bkgrd, file_out)
