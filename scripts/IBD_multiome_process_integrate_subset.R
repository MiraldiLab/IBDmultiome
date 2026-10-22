# IBD_multiome_process_integrate_subset.R
# module load R/4.2.3 bedtools/2.30.0 hdf5/1.8.13 igvtools/2.4.0
# Reintegrate selected cell-type groups for finer downstream annotation.
# Run RNA and ATAC Harmony integration and, when both are enabled, joint
# weighted nearest-neighbor (WNN) analysis.
#
# Replace the /path/to placeholders before running. The input object must
# contain RNA and ATAC count layers and the configured metadata covariates.
# The annotation table must be tab-delimited, with cell names in its first
# column and an annotation column matching name_annot.
# tissue_subset labels output directories and files; it does not filter cells.
# The RNA-enabled metadata plotting block expects umap_WNN to be available;
# the default configuration enables both modalities to generate it.
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
dir_out <- '/path/to/output/integrate_subset'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Integrated multiome object
file_obj <- '/path/to/integrate_multiome/All/obj_multiomeIBD_Multi_harmony_All.rds'

# Cell annotation
file_annot <- '/path/to/cell_annotations/annot_cell_All_wsnn_res.1_CellTypeL0.txt'
name_annot <- 'CellTypeL0'

# Analysis parameters
flag_scrna <- TRUE # Perform scRNA integration
flag_scatac <- TRUE # Perform scATAC integration
file_save <- 'multiomeIBD' # file save base
tissue_subset <- 'All' # name subset of data
scrna_min_cells <- 20 # Minimum number of cells to keep gene
scatac_min_cells <- 100 # minimum number of cells for feature to be considered in dimension reduction
n_var_feat <- 2000 # Number of variable features for integration/dimension reduction
vars_regress <- c('nCount_RNA', 'gex_pct_mito') # RNA variable to regress off during scaling
k_pca <- 30 # PCA dimension reduction
k_lsi <- 30 # ATAC dimension reduction
cutoff_cor_lsi <- 0.75 # cutoff correlation with count depth to exclude LSI components
vars_harmony <- c('Patient') # Harmony covariates
iter_harmony <- 30 # max iterations for Harmony

# Optional text file containing one cell name per line; use NULL for none.
# Example: file_rm_cells <- '/path/to/cell_annotations/cells_to_remove.txt'
file_rm_cells <- NULL

# Named groups to reintegrate; values must match labels in name_annot.
subcluster_group <- list()
subcluster_group[['B_cell']] <- c('B_cell')
subcluster_group[['T_cell']] <- c('T_cell')
subcluster_group[['Epithelial']] <- c('Epithelial')
# subcluster_group[['Goblet']] <- c('Goblet')
# subcluster_group[['Fibroblast']] <- c('Fibroblast')
# subcluster_group[['Myeloid']] <- c('Myeloid')
# subcluster_group[['Plasma_cell']] <- c('Plasma_cell')
# subcluster_group[['Enteroendocrine']] <- c('Enteroendocrine')
# subcluster_group[['Endothelial']] <- c('Endothelial')
# subcluster_group[['Tuft']] <- c('Tuft')
# subcluster_group[['Cycling']] <- c('Cycling')

# Metadata to plot
plot_meta_cell <- c('nCount_RNA', 'nFeature_RNA', 'gex_pct_mito', 'gex_pct_ribo',
                    'atac_fragments', 'nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment',
                    'atac_pct_mito', 'atac_pct_blacklist', 'nucleosome_signal', 'atac_pct_reads_in_peaks')
plot_meta_sample <- c('Disease', 'Tissue', 'Status', 'Sex', 'Patient', 'AILD', 'Sample')

# ---- Load data and annotations ------------------------------------------------

# Create the output directory for this data label.
dir_out <- file.path(dir_out, tissue_subset)
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Load object
print('Load object')
obj_all <- readRDS(file_obj)

# Cell-type metadata
meta_celltype <- read.delim(file_annot, header = TRUE, row.names = 1, sep = '\t')
obj_all <- AddMetaData(obj_all, meta_celltype)
Idents(obj_all) <- name_annot

# Cells to remove
if (!is.null(file_rm_cells)) {
    rm_cells <- readLines(file_rm_cells)
} else {
    rm_cells <- NULL
}

# ---- Reintegrate each selected group -------------------------------------------
print('Integrate groups')
for (ix in names(subcluster_group)) {

    # Create the group directory and select cells, excluding removal-list entries.
    dir_out_sub <- file.path(dir_out, ix)
    dir.create(dir_out_sub, showWarnings = FALSE)
    obj_sub <- subset(obj_all, cells = setdiff(WhichCells(obj_all, idents = subcluster_group[[ix]]), rm_cells))

    # Save parameters
    vparam <- c(paste0('Tissue: ', tissue_subset), paste0('Genes min cells: ', scrna_min_cells),
                paste0('Peaks min cells: ', scatac_min_cells), paste0('N Variable Genes: ', n_var_feat),
                paste0('k_pca = ', k_pca), paste0('k_lsi = ', k_lsi))
    file_out <- file.path(dir_out_sub, paste0('param_', file_save, '_sub_', ix, '.txt'))
    writeLines(vparam, file_out)

    # Retain RNA and ATAC count layers for reintegration.
    obj_sub <- DietSeurat(obj_sub, assays = c('RNA', 'ATAC'), layers = c('counts', 'counts'))

    # Filter RNA genes by detection count and retain ATAC features with min.cutoff = 1.
    if ('RNA' %in% names(obj_sub@assays)) {
        keep_feat_scrna <- rownames(CreateSeuratObject(LayerData(obj_sub, assay = 'RNA', layer = 'counts'),
                                    min.cells = scrna_min_cells))
    } else {
        keep_feat_scrna <- NULL
    }
    DefaultAssay(obj_sub) <- 'ATAC'
    keep_feat_scatac <- VariableFeatures(FindTopFeatures(obj_sub, min.cutoff = 1, verbose = FALSE))
    keep_feat <- c(keep_feat_scrna, keep_feat_scatac)
    obj_sub <- subset(obj_sub, features = keep_feat)

    if (flag_scrna) {
        # Normalize RNA, scale covariates, and integrate PCA embeddings with Harmony.
        print('- Harmony RNA')
        DefaultAssay(obj_sub) <- 'RNA'
        obj_sub <- FindVariableFeatures(obj_sub, nfeatures = n_var_feat, verbose = TRUE)
        VariableFeatures(obj_sub) <- VariableFeatures(obj_sub)[1:n_var_feat]
        obj_sub <- NormalizeData(obj_sub, assay = 'RNA', verbose = FALSE)
        obj_sub <- ScaleData(obj_sub, vars.to.regress = vars_regress, verbose = FALSE)
        obj_sub <- RunPCA(obj_sub, verbose = FALSE)
        obj_sub <- RunHarmony(object = obj_sub, group.by.vars = vars_harmony, reduction.use = 'pca',
                                    assay.use = 'RNA', reduction.save = 'harmony_RNA', project.dim = FALSE,
                                    max_iter = iter_harmony, lambda = rep(1, length(vars_harmony)))
        obj_sub <- RunUMAP(obj_sub, dims = 1:k_pca, reduction = 'harmony_RNA', reduction.name = 'umap_RNA',
                            verbose = FALSE, return.model = TRUE)
        DimPlot(obj_sub, label = TRUE, reduction = 'umap_RNA', shuffle = TRUE) + theme_void() + NoLegend()
        file_out <- file.path(dir_out_sub, paste0('umap_', file_save, '_sub_', ix, '_harmony_RNA_', name_annot, '.png'))
        ggsave(file_out)
    }

    if (flag_scatac) {
        # Normalize ATAC and integrate LSI embeddings with Harmony.
        print('- Harmony ATAC')
        DefaultAssay(obj_sub) <- 'ATAC'
        obj_sub <- RunTFIDF(obj_sub, verbose = FALSE)
        obj_sub <- FindTopFeatures(obj_sub, min.cutoff = scatac_min_cells, verbose = FALSE)
        obj_sub <- RunSVD(obj_sub, n = k_lsi, reduction.name = 'lsi', reduction.key = 'LSI_', verbose = FALSE)
        obj_sub <- RunHarmony(object = obj_sub, group.by.vars = vars_harmony, reduction.use = 'lsi',
                                    assay.use = 'ATAC', reduction.save = 'harmony_ATAC', project.dim = FALSE,
                                    max_iter = iter_harmony, lambda = rep(1, length(vars_harmony)))
        cor_depth <- DepthCor(obj_sub, assay = 'ATAC', reduction = 'harmony_ATAC', n = k_lsi)
        keep_dim_lsi <- which(abs(cor_depth$data$nCount_ATAC) < cutoff_cor_lsi)
        obj_sub <- RunUMAP(obj_sub, dims = keep_dim_lsi, reduction = 'harmony_ATAC',
                            reduction.name = 'umap_ATAC', verbose = FALSE, return.model = TRUE)
        DimPlot(obj_sub, label = TRUE, reduction = 'umap_ATAC', shuffle = TRUE) + theme_void() + NoLegend()
        file_out <- file.path(dir_out_sub, paste0('umap_', file_save, '_sub_', ix, '_harmony_ATAC_', name_annot, '.png'))
        ggsave(file_out)
    }

    if (flag_scrna & flag_scatac) {
        # Combine RNA and ATAC embeddings using weighted nearest neighbors.
        print('- WNN')
        obj_sub <- FindMultiModalNeighbors(obj_sub, reduction.list = list('harmony_RNA', 'harmony_ATAC'),
                                            dims.list = list(1:k_pca, keep_dim_lsi), verbose = FALSE)
        obj_sub <- RunUMAP(obj_sub, nn.name = 'weighted.nn', reduction.name = 'umap_WNN',
                                reduction.key = 'wnnUMAP_', verbose = FALSE, return.model = TRUE)
        DimPlot(obj_sub, label = TRUE, reduction = 'umap_WNN', shuffle = TRUE) + theme_void() + NoLegend()
        file_out <- file.path(dir_out_sub, paste0('umap_', file_save, '_sub_', ix, '_harmony_WNN_', name_annot, '.png'))
        ggsave(file_out)
    }

    if (flag_scrna) {
        # Cell cycle scores
        DefaultAssay(obj_sub) <- 'RNA'
        obj_sub <- CellCycleScoring(obj_sub, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes,
                                    set.ident = FALSE)
        # Plot continuous and categorical metadata on the joint WNN embedding.
        for (jx in c(plot_meta_cell, 'G2M.Score', 'S.Score')) {
            FeaturePlot(obj_sub, features = jx, reduction = 'umap_WNN', order = TRUE, min.cutoff = 'q5', max.cutoff = 'q95') +
                            labs(x = 'UMAP-1', y = 'UMAP-2')
            file_out <- file.path(dir_out_sub, paste0('umap_', file_save, '_sub_', ix, '_harmony_WNN_feat_', jx, '.png'))
            ggsave(file_out, width = 7)
        }
        for (jx in intersect(plot_meta_sample, colnames(obj_sub@meta.data))) {
            DimPlot(obj_sub, group.by = jx, reduction = 'umap_WNN', shuffle = TRUE) + labs(x = 'UMAP-1', y = 'UMAP-2')
            file_out <- file.path(dir_out_sub, paste0('umap_', file_save, '_sub_', ix, '_harmony_WNN_meta_', jx, '.png'))
            ggsave(file_out, width = 7)

            DimPlot(obj_sub, group.by = jx, reduction = 'umap_WNN', shuffle = TRUE, label = TRUE, label.size = 4, repel = TRUE) +
                    theme_void() + NoLegend()
            file_out <- file.path(dir_out_sub,
                                    paste0('umap_', file_save, '_sub_', ix, '_harmony_WNN_meta_', jx, '_label.png'))
            ggsave(file_out)
        }
    }

    # Save integrated object
    print('- Save integrated object')
    file_out <- file.path(dir_out_sub, paste0('obj_', file_save, '_sub_', ix, '_harmony_', tissue_subset, '.rds'))
    saveRDS(obj_sub, file_out)
}

print('DONE')
