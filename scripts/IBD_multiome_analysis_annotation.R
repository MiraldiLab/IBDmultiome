# IBD_multiome_analysis_annotation.R
# module load R/4.2.3 bedtools/2.30.0 hdf5/1.8.13 igvtools/2.4.0
# Map reviewed cluster annotations to cells and visualize them on the WNN UMAP.
# Export a cell-level annotation table and cluster/annotation plots.
#
# Replace the /path/to placeholders before running.
# file_clust: two tab-delimited columns without a header (cell name, cluster).
# file_annot: two tab-delimited columns without a header (cluster, annotation).
# Cell names must match the multiome object, which must contain an RNA assay
# and the umap_WNN reduction. Cluster identifiers must match between tables.
# tissue_subset labels outputs; it does not subset the input object.
# The annotated object is used for plotting but is not saved as an RDS file.
rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(42)
suppressPackageStartupMessages({
    library(Seurat)
    library(Signac)
    library(GenomicRanges)
    library(ggplot2)
})

# ---- Input paths and annotation parameters -------------------------------------

# Output directory
dir_out <- '/path/to/output/annotation'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Multiome object
file_data <- '/path/to/integrate_multiome/All/obj_multiomeIBD_Multi_harmony_All.rds'

# Cell clustering
name_clust <- 'wsnn_res.1'
file_clust <- '/path/to/clustering/meta_multiomeIBD_Multi_All_clusterRes1.txt'

# Cluster annotation
name_annot <- 'CellTypeL0'
file_annot <- '/path/to/cluster_annotations/annot_clust_All_wsnn_res.1_CellTypeL0.txt'

# Output labels
file_save <- 'multiomeIBD' # file save base
tissue_subset <- 'All' # name subset of data

# ---- Load clustering and reviewed annotations ----------------------------------

# Output directory
dir_out <- file.path(dir_out, tissue_subset)
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Load cell clustering
cell_clust <- read.delim(file_clust, sep = '\t', header = FALSE)
colnames(cell_clust) <- c('Cell', name_clust)

# Load cluster annotation
annot_clust <- read.delim(file_annot, sep = '\t', header = FALSE)
colnames(annot_clust) <- c(name_clust, name_annot)

# Map clusters to annotations using an inner join on the cluster identifier.
# Only matching entries are retained in the exported annotation table.
meta_annot <- merge(cell_clust, annot_clust, by = name_clust)
rownames(meta_annot) <- meta_annot$Cell
meta_annot <- meta_annot[, c(name_clust, name_annot)]

# ---- Export cell-level annotations ---------------------------------------------
file_out <- file.path(dir_out, paste0('annot_cell_', tissue_subset, '_', name_clust, '_', name_annot, '.txt'))
write.table(meta_annot, file_out, quote = FALSE, sep = '\t', col.names = NA, row.names = TRUE)

# ---- Attach annotations to the multiome object ---------------------------------
print('Load multiome data')
obj <- readRDS(file_data)
DefaultAssay(obj) <- 'RNA'
obj <- AddMetaData(obj, meta_annot)

# ---- Plot clusters and annotations ---------------------------------------------

# Order cluster identities by decreasing cell count.
Idents(obj) <- name_clust
Idents(obj) <- factor(Idents(obj), levels = names(sort(table(Idents(obj)), decreasing = TRUE)))
DimPlot(obj, reduction = 'umap_WNN', label = TRUE, label.size = 6, shuffle = TRUE) + theme_void() + NoLegend()
file_out <- file.path(dir_out, paste0('umap_', tissue_subset, '_WNN_clust_', name_clust, '.png'))
ggsave(file_out, height = 7, width = 7)

# Order annotation identities by decreasing cell count.
Idents(obj) <- name_annot
Idents(obj) <- factor(Idents(obj), levels = names(sort(table(Idents(obj)), decreasing = TRUE)))
DimPlot(obj, reduction = 'umap_WNN', label = TRUE, label.size = 6, shuffle = TRUE) + theme_void() + NoLegend()
file_out <- file.path(dir_out, paste0('umap_', tissue_subset, '_WNN_clust_', name_clust, '_annot_', name_annot, '.png'))
ggsave(file_out, height = 7, width = 7)

print('DONE')
