# IBD_multiome_process_integrate_multiome.R
# module load R/4.2.3 bedtools/2.30.0 hdf5/1.8.13 igvtools/2.4.0
# Integrate multiome cells that pass both RNA and ATAC quality control.
# Perform Harmony integration for each modality, followed by weighted
# nearest-neighbor (WNN) analysis.
#
# Before running, replace the /path/to placeholders below.
# Sample metadata must include Sample and Tissue. Input objects must provide
# Sample and barcode metadata, plus the regression, Harmony, and plotting
# covariates configured below. The merge workflow expects multiple samples.
# The reference ATAC object must contain matching sample-prefixed cell names;
# select a reference object appropriate for tissue_subset.
rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(42)
suppressPackageStartupMessages({
    library(Seurat)
    library(Signac)
    library(harmony)
    library(GenomicRanges)
    library(ggplot2)
    library(magrittr)
    library(sva)
    library(Matrix)
    library(rtracklayer)
    library(scales)
})

# ---- Input paths and analysis parameters --------------------------------------

# Output directory
dir_out <- '/path/to/output/integrate_multiome'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Sample metadata
file_meta <- '/path/to/sample_metadata/meta_multiomeIBD_sample.txt'

# Directory of filtered RNA+ATAC objects
dir_filtered <- '/path/to/qc_filter'

# ATAC merged object
file_scatac <- '/path/to/integrate_atac/TI/obj_multiomeIBD_ATAC_merge_TI.rds'

# Analysis parameters
file_save <- 'multiomeIBD' # file save base
file_ref_prefix <- 'obj_multiomeIBD'
tissue_subset <- 'TI' # tissue subset ('All' for all cells)
min_cells_scrna <- 50 # RNA gene filter - min cells
min_cells_scatac <- 100 # ATAC peak filter - min cells
vars_regress <- c('nCount_RNA', 'gex_pct_mito') # RNA variable to regress off during scaling
k_lsi <- 30 # ATAC dimension reduction
k_pca <- 30 # PCA dimension reduction
cutoff_cor_lsi <- 0.75 # cutoff correlation with count depth to exclude LSI components
vars_harmony <- c('Patient') # Harmony covariates
iter_harmony <- 30 # maximum Harmony iterations
n_var_feat <- 2000 # Number of variable features for integration/dimension reduction

# Metadata to plot
plot_meta <- c('Disease', 'Tissue', 'Status', 'Sex', 'Patient', 'AILD', 'Sample')

# ---- Load and merge data ------------------------------------------------------

# Create the tissue-specific output directory.
dir_out <- file.path(dir_out, tissue_subset)
dir.create(dir_out, showWarnings = FALSE)

# Load sample metadata
meta_sample <- read.delim(file_meta, header = TRUE, sep = '\t')
if (tissue_subset != 'All') {
    meta_sample <- subset(meta_sample, Tissue == tissue_subset)
}
name_sample <- meta_sample$Sample

# Merge RNA+ATAC objects
print('Object merge')
sc_list <- list()
for (ix in 1:length(name_sample)) {
    file_data <- file.path(dir_filtered, name_sample[ix], paste0(file_ref_prefix, '_', name_sample[ix], '_Multi.rds'))
    curr_scrna <- readRDS(file_data)
    DefaultAssay(curr_scrna) <- 'RNA'
    curr_scrna <- RenameCells(curr_scrna, new.names = paste(as.character(curr_scrna@meta.data$Sample),
                                                            as.character(curr_scrna@meta.data$barcode), sep = '_'))
    curr_scrna <- DietSeurat(curr_scrna, assays = 'RNA', layers = 'counts')
    sc_list[[name_sample[ix]]] <- curr_scrna
}
rm(curr_scrna)
sc_merge <- merge(x = sc_list[[1]], y = sc_list[2:length(sc_list)], merge.data = FALSE)
sc_merge <- JoinLayers(sc_merge)

# Load ATAC object
print('Load ATAC object')
scatac <- readRDS(file_scatac)
scatac <- DietSeurat(scatac, assays = 'ATAC', layers = 'counts')

# Add the reference ATAC assay, retaining cells present in the RNA merge.
sc_merge[['ATAC']] <- subset(scatac, cells = colnames(sc_merge))[['ATAC']]
rm(scatac)

# ---- Feature filtering and selection ------------------------------------------

# Retain genes and peaks detected in the configured minimum number of cells.
keep_feat_scrna <- rownames(CreateSeuratObject(LayerData(sc_merge, assay = 'RNA', layer = 'counts'),
                                min.cells = min_cells_scrna))
keep_feat_scatac <- rownames(CreateSeuratObject(LayerData(sc_merge, assay = 'ATAC', layer = 'counts'),
                                min.cells = min_cells_scatac))
sc_merge <- subset(sc_merge, features = c(keep_feat_scrna, keep_feat_scatac))
DefaultAssay(sc_merge) <- 'RNA'

# Select shared integration features for RNA dimensionality reduction.
feat_int <- SelectIntegrationFeatures(sc_list, nfeatures = 10000)
feat_int <- intersect(feat_int, rownames(sc_merge))[1:n_var_feat]
VariableFeatures(sc_merge) <- feat_int
rm(sc_list)

# ---- RNA normalization and Harmony integration --------------------------------
print('Harmony RNA')
# sc_merge <- FindVariableFeatures(sc_merge, verbose=FALSE)
sc_merge <- NormalizeData(sc_merge, assay = 'RNA', verbose = FALSE)
sc_merge <- ScaleData(sc_merge, vars.to.regress = vars_regress, verbose = FALSE)
sc_merge <- RunPCA(sc_merge, verbose = FALSE)
sc_merge <- RunHarmony(object = sc_merge, group.by.vars = vars_harmony, reduction.use = 'pca',
                            assay.use = 'RNA', reduction.save = 'harmony_RNA', project.dim = FALSE,
                            max_iter = iter_harmony, lambda = rep(1, length(vars_harmony)))
sc_merge <- RunUMAP(sc_merge, dims = 1:k_pca, reduction = 'harmony_RNA', reduction.name = 'umap_RNA',
                        verbose = FALSE, return.model = TRUE)
for (ix in plot_meta) {
    DimPlot(sc_merge, group.by = ix, shuffle = TRUE) + labs(x = 'UMAP-1', y = 'UMAP-2')
    file_out <- file.path(dir_out, paste0('umap_', file_save, '_Multi_harmony_RNA_', ix, '.png'))
    ggsave(file_out, width = 7)
}

# ---- ATAC normalization and Harmony integration -------------------------------
print('Harmony ATAC')
DefaultAssay(sc_merge) <- 'ATAC'
sc_merge <- RunTFIDF(sc_merge, verbose = FALSE)
sc_merge <- FindTopFeatures(sc_merge, min.cutoff = 'q0', verbose = FALSE)
sc_merge <- RunSVD(sc_merge, n = k_lsi, reduction.name = 'lsi', reduction.key = 'LSI_', verbose = FALSE)
sc_merge <- RunHarmony(object = sc_merge, group.by.vars = vars_harmony, reduction.use = 'lsi',
                            assay.use = 'ATAC', reduction.save = 'harmony_ATAC', project.dim = FALSE,
                            max_iter = iter_harmony, lambda = rep(1, length(vars_harmony)))
cor_depth <- DepthCor(sc_merge, assay = 'ATAC', reduction = 'harmony_ATAC', n = k_lsi)
keep_dim_lsi <- which(abs(cor_depth$data$nCount_ATAC) < cutoff_cor_lsi)
sc_merge <- RunUMAP(sc_merge, dims = keep_dim_lsi, reduction = 'harmony_ATAC', reduction.name = 'umap_ATAC',
                        verbose = FALSE, return.model = TRUE)
for (ix in plot_meta) {
    DimPlot(sc_merge, group.by = ix, shuffle = TRUE) + labs(x = 'UMAP-1', y = 'UMAP-2')
    file_out <- file.path(dir_out, paste0('umap_', file_save, '_Multi_harmony_ATAC_', ix, '.png'))
    ggsave(file_out, width = 7)
}

# ---- Joint weighted nearest-neighbor analysis ---------------------------------
print('WNN')
sc_merge <- FindMultiModalNeighbors(sc_merge, reduction.list = list('harmony_RNA', 'harmony_ATAC'),
                                    dims.list = list(1:k_pca, keep_dim_lsi))
sc_merge <- RunUMAP(sc_merge, nn.name = 'weighted.nn', reduction.name = 'umap_WNN', reduction.key = 'wnnUMAP_',
                        verbose = FALSE, return.model = TRUE)
for (ix in plot_meta) {
    DimPlot(sc_merge, group.by = ix, reduction = 'umap_WNN', shuffle = TRUE) + labs(x = 'UMAP-1', y = 'UMAP-2')
    file_out <- file.path(dir_out, paste0('umap_', file_save, '_Multi_harmony_WNN_', ix, '.png'))
    ggsave(file_out, width = 7)
}
DefaultAssay(sc_merge) <- 'RNA'

# ---- Save the integrated object -----------------------------------------------
file_out <- file.path(dir_out, paste0('obj_', file_save, '_Multi_harmony_', tissue_subset, '.rds'))
saveRDS(sc_merge, file_out)

print('DONE')
