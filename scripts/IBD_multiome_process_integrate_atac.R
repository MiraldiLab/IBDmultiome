# IBD_multiome_process_integrate_atac.R
#
# Purpose: Merge sample ATAC objects quantified against a common peak reference,
# integrate LSI embeddings with Harmony, and export objects and diagnostic UMAPs.
#
# Setup: Replace all /path/to/... placeholders with local paths.
# Sample metadata must contain Sample and Tissue columns.
# Input objects: <dir_ref>/<Sample>/<file_ref_prefix>_<Sample>_ATAC.rds.
# Objects must contain an ATAC counts layer and the metadata columns specified
# in vars_harmony and plot_meta. The original workflow expects multiple samples.
# Optional file_rm_cells contains one merged cell ID per line, including the
# sample prefix added by merge(); leave NULL to skip additional cell removals.
#
# Original environment:
# module load R/4.2.3 bedtools/2.30.0 hdf5/1.8.13 igvtools/2.4.0

# Session setup and packages --------------------------------------------------
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

# Inputs and integration parameters ------------------------------------------

# Output directory
dir_out <- '/path/to/output/integrate_atac'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Sample metadata
file_meta <- '/path/to/sample_metadata/meta_multiomeIBD_sample.txt'

# Directory containing ATAC objects quantified against the common peak reference
dir_ref <- '/path/to/atac_peaks_ref/sizeFrag200_FDR0p1_lambda_llocal10000'

# Analysis parameters
file_save <- 'multiomeIBD' # file save base
file_ref_prefix <- 'obj_multiomeIBD'
tissue_subset <- 'TI' # tissue subset ('All' for all cells)
min_cells_scatac <- 100 # Minimum cells in which a peak must be detected to retain it
k_lsi <- 30 # Number of LSI components for the merged object
vars_harmony <- c('Patient') # Harmony covariates
iter_harmony <- 30 # maximum Harmony iterations
cutoff_cor_lsi <- 0.75 # cutoff correlation with count depth to exclude LSI components
file_rm_cells <- NULL # cells to remove (low quality)
plot_meta <- c('Disease', 'Tissue', 'Status', 'Sex', 'Patient', 'AILD', 'Sample') # metadata to plot

keep_layers <- list()
keep_layers[['ATAC']] <- c('counts')

# Configure outputs and select samples ---------------------------------------

# Output directory
dir_out <- file.path(dir_out, tissue_subset)
dir.create(dir_out, showWarnings = FALSE)

# Load sample metadata
meta_sample <- read.delim(file_meta, header = TRUE, sep = '\t')
if (tissue_subset != 'All') {
    meta_sample <- subset(meta_sample, Tissue == tissue_subset)
}
name_sample <- meta_sample$Sample

# Load sample objects and retain the requested assay layers -------------------
print('Object list')
sc_list <- list()
for (ix in 1:length(name_sample)) {
    file_data <- file.path(dir_ref, name_sample[ix], paste0(file_ref_prefix, '_', name_sample[ix], '_ATAC.rds'))
    curr_scatac <- readRDS(file_data)
    DefaultAssay(curr_scatac) <- 'ATAC'
    curr_scatac <- DietSeurat(curr_scatac, assays = names(keep_layers), layers = keep_layers)
    sc_list[[name_sample[ix]]] <- curr_scatac
}
rm(curr_scatac)

# Merge objects with sample-prefixed cell identifiers -------------------------
print('Merge objects')
sc_merge <- merge(x = sc_list[[1]], y = sc_list[2:length(sc_list)], merge.data = FALSE, add.cell.ids = names(sc_list))
rm(sc_list)

# Remove cells
if (!is.null(file_rm_cells)) {
    cells_rm <- readLines(file_rm_cells)
    sc_merge <- subset(sc_merge, cells = setdiff(colnames(sc_merge), cells_rm))
}

# Filter features
print('Filter features')
keep_feat_scatac <- rownames(CreateSeuratObject(LayerData(sc_merge, assay = 'ATAC', layer = 'counts'),
                                min.cells = min_cells_scatac))
sc_merge <- subset(sc_merge, features = keep_feat_scatac)

# Compute LSI and integrate ATAC embeddings with Harmony ----------------------
print('Harmony ATAC')
sc_merge <- RunTFIDF(sc_merge, verbose = FALSE)
sc_merge <- FindTopFeatures(sc_merge, min.cutoff = 'q0', verbose = FALSE)
sc_merge <- RunSVD(sc_merge, n = k_lsi, reduction.name = 'lsi', reduction.key = 'LSI_', verbose = FALSE)
sc_integrated <- RunHarmony(object = sc_merge, group.by.vars = vars_harmony, reduction.use = 'lsi',
                            assay.use = 'ATAC', reduction.save = 'harmony_ATAC', project.dim = FALSE,
                            max_iter = iter_harmony, lambda = rep(1, length(vars_harmony)))
cor_depth <- DepthCor(sc_integrated, assay = 'ATAC', reduction = 'harmony_ATAC', n = k_lsi)
keep_dim_lsi <- which(abs(cor_depth$data$nCount_ATAC) < cutoff_cor_lsi)
sc_integrated <- RunUMAP(sc_integrated, dims = keep_dim_lsi, reduction = 'harmony_ATAC', reduction.name = 'umap_ATAC',
                            verbose = FALSE, return.model = TRUE)
for (ix in plot_meta) {
    DimPlot(sc_integrated, group.by = ix, shuffle = TRUE) + labs(x = 'UMAP-1', y = 'UMAP-2')
    file_out <- file.path(dir_out, paste0('umap_', file_save, '_ATAC_harmony_', ix, '.png'))
    ggsave(file_out, height = 6, width = 6.5)
}

# Export integrated and merged objects ---------------------------------------
print('Save objects')
file_out <- file.path(dir_out, paste0('obj_', file_save, '_ATAC_harmony_', tissue_subset, '.rds'))
saveRDS(sc_integrated, file_out)
file_out <- file.path(dir_out, paste0('obj_', file_save, '_ATAC_merge_', tissue_subset, '.rds'))
saveRDS(sc_merge, file_out)

# Export cell metadata as tab-delimited text.
# The original .rds filename extension is retained; this is not an RDS object.
file_out <- file.path(dir_out, paste0('meta_', file_save, '_ATAC_harmony_', tissue_subset, '.rds'))
write.table(sc_integrated@meta.data, file_out, quote = FALSE, sep = '\t', col.names = NA, row.names = TRUE)

print('DONE')
