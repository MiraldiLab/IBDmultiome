# IBD_multiome_analysis_cluster.R
# module load R/4.2.3 bedtools/2.30.0 hdf5/1.8.13 igvtools/2.4.0
# Visualize transferred labels and cluster the integrated multiome WNN graph.
# Export cluster assignments and label-composition heatmaps.
#
# Replace the /path/to placeholders before running, including the Python
# executable used for Leiden clustering.
# The ATAC object must include Sample metadata and a UMAP embedding; the
# multiome object must include the wsnn graph and umap_WNN reduction.
# Transfer tables must contain the configured predicted-label columns and
# unprefixed cell barcodes as row names; sample prefixes are added below.
# tissue_subset labels outputs; it does not subset either input object.
#
# Existing helper dependency: net_sparse_2_full() and cluster_hierarchical()
# must be available when the heatmap block runs. This script does not define
# or explicitly source these functions.
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
    library(reticulate)
})

# ---- Input paths and clustering parameters -------------------------------------

# Output directory
dir_out <- '/path/to/output/cluster'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# ATAC object
file_scatac <- '/path/to/integrate_atac/All/obj_multiomeIBD_ATAC_harmony_All.rds'

# Integrated RNA+ATAC object
file_scrna <- '/path/to/integrate_multiome/All/obj_multiomeIBD_Multi_harmony_All.rds'

# Label transfer names
name_transfer <- c(
'Elmentaite2021',
'Smillie2019'
)
# Label transfer metadata directory
dir_meta_transfer <- '/path/to/label_transfer_metadata'

# Clustering parameters
file_save <- 'multiomeIBD' # file save base
tissue_subset <- 'All' # name subset of data
res_clust <- 1 # clustering resolution scRNA+scATAC
algo_clust <- 4 # 1 = Louvain, 3 = SLM algorithm, 4 = Leiden

# ---- Configure Python and load integrated objects ------------------------------

# Select the Python executable for Leiden clustering.
if (algo_clust == 4) {
    use_python('/path/to/python/environment/bin/python')
}

# Output directory
dir_out <- file.path(dir_out, tissue_subset)
dir.create(dir_out, showWarnings = FALSE)

# Load ATAC object
print('Load ATAC object')
scatac <- readRDS(file_scatac)

# Load RNA+ATAC object
print('Load RNA+ATAC object')
scrna <- readRDS(file_scrna)

# ---- Attach transferred labels ------------------------------------------------

# Read per-sample labels and match sample-prefixed cell names.
meta_transfer <- NULL
for (ix in unique(scatac@meta.data$Sample)) {
    curr_file <- file.path(dir_meta_transfer, paste0('meta_transfer_', file_save, '_', ix, '_Multi_RNA.txt'))
    curr_meta <- read.delim(curr_file, header = TRUE, row.names = 1)
    col_transfer <- paste0('Transfer_', name_transfer)
    col_transfer <- paste0(col_transfer, '_predicted.id_filtered')
    curr_meta <- curr_meta[, col_transfer]
    rownames(curr_meta) <- paste(ix, rownames(curr_meta), sep = '_')
    meta_transfer <- rbind(meta_transfer, curr_meta)
}
scatac <- AddMetaData(scatac, meta_transfer)
scrna <- AddMetaData(scrna, meta_transfer)
for (ix in name_transfer) {
    curr_name <- paste0('Transfer_', ix, '_predicted.id_filtered')
    scatac@meta.data[which(is.na(scatac@meta.data[, curr_name])), curr_name] <- ' '
    scatac@meta.data[which(scatac@meta.data[, curr_name] == 'Unclassified'), curr_name] <- ' '
    scrna@meta.data[which(is.na(scrna@meta.data[, curr_name])), curr_name] <- ' '
    scrna@meta.data[which(scrna@meta.data[, curr_name] == 'Unclassified'), curr_name] <- ' '
}
for (ix in name_transfer) {
    curr_name <- paste0('Transfer_', ix, '_predicted.id_filtered')
    curr_levels <- names(sort(table(scatac@meta.data[, curr_name]), decreasing = TRUE))
    curr_levels <- c(setdiff(curr_levels, ' '), ' ')
    scatac@meta.data[, curr_name] <- factor(scatac@meta.data[, curr_name], levels = curr_levels)
    curr_levels <- names(sort(table(scrna@meta.data[, curr_name]), decreasing = TRUE))
    curr_levels <- c(setdiff(curr_levels, ' '), ' ')
    scrna@meta.data[, curr_name] <- factor(scrna@meta.data[, curr_name], levels = curr_levels)
}

# ---- Plot transferred labels --------------------------------------------------

# ATAC label overview and individual label highlights.
for (ix in name_transfer) {
    curr_name <- paste0('Transfer_', ix, '_predicted.id_filtered')
    curr_levels <- levels(scatac@meta.data[, curr_name])
    curr_order <- NULL
    Idents(scatac) <- curr_name
    curr_order <- WhichCells(scatac, idents = setdiff(levels(scatac), ' '))
    curr_order <- c(curr_order, WhichCells(scatac, idents = ' '))
    DimPlot(scatac, group.by = curr_name, label = TRUE, repel = TRUE, pt.size = 0.5, cells = rev(curr_order),
            cols = c(hue_pal()(length(levels(scatac@meta.data[, curr_name]))-1), 'grey90'), shuffle = TRUE) +
            NoLegend() + labs(x = 'UMAP-1', y = 'UMAP-2')
    file_out <- file.path(dir_out, paste0('umap_', file_save, '_ATAC_transfer_RNA_', ix, '.png'))
    ggsave(file_out)
    dir_out_transfer <- file.path(dir_out, paste0('Transfer_ATAC_', ix))
    dir.create(dir_out_transfer, showWarnings = FALSE)
    for (jx in levels(scatac)) {
        DimPlot(scatac, pt.size = 0.5, sizes.highlight = 0.5, cells.highlight = WhichCells(scatac, idents = jx),
                    cols.highlight = 'grey20') +
                    labs(x = 'UMAP-1', y = 'UMAP-2') + theme_void() + NoLegend()
        file_out <- file.path(dir_out_transfer, paste0('umap_', file_save, '_ATAC_transfer_RNA_', ix, '_', jx, '.png'))
        ggsave(file_out)
    }
}
# RNA+ATAC label overview and individual label highlights.
for (ix in name_transfer) {
    curr_name <- paste0('Transfer_', ix, '_predicted.id_filtered')
    curr_levels <- levels(scrna@meta.data[, curr_name])
    curr_order <- NULL
    Idents(scrna) <- curr_name
    curr_order <- WhichCells(scrna, idents = setdiff(levels(scrna), ' '))
    if (' ' %in% levels(scrna)) {
        curr_order <- c(curr_order, WhichCells(scrna, idents = ' '))
    }
    DimPlot(scrna, reduction = 'umap_WNN', group.by = curr_name, label = TRUE, repel = TRUE, pt.size = 0.5, cells = rev(curr_order),
            cols = c(hue_pal()(length(levels(scrna@meta.data[, curr_name]))-1), 'grey90'), shuffle = TRUE) +
            NoLegend() + labs(x = 'UMAP-1', y = 'UMAP-2')
    file_out <- file.path(dir_out, paste0('umap_', file_save, '_Multi_transfer_RNA_', ix, '.png'))
    ggsave(file_out)
    dir_out_transfer <- file.path(dir_out, paste0('Transfer_Multi_', ix))
    dir.create(dir_out_transfer, showWarnings = FALSE)
    for (jx in levels(scrna)) {
        DimPlot(scrna, reduction = 'umap_WNN', pt.size = 0.5, sizes.highlight = 0.5,
                    cells.highlight = WhichCells(scrna, idents = jx), cols.highlight = 'grey20') +
                    labs(x = 'UMAP-1', y = 'UMAP-2') + theme_void() + NoLegend()
        file_out <- file.path(dir_out_transfer, paste0('umap_', file_save, '_Multi_transfer_RNA_', ix, '_', jx, '.png'))
        ggsave(file_out)
    }
}

# ---- Cluster the multiome WNN graph and summarize labels ------------------------
print('Cluster RNA+ATAC')
scrna <- FindClusters(scrna, graph.name = 'wsnn', resolution = res_clust, algorithm = algo_clust,
                        verbose = FALSE, method = ifelse(algo_clust == 4, 'igraph', 'matrix'))
DimPlot(scrna, reduction = 'umap_WNN', label = TRUE, label.size = 6, shuffle = TRUE) + labs(x = 'UMAP-1', y = 'UMAP-2') +
            theme_void() + NoLegend()
file_out <- file.path(dir_out, paste0('umap_', file_save, '_Multi_clusterRes', res_clust, '.png'))
ggsave(file_out)
for (ix in name_transfer) {
    df_type <- NULL
    curr_name <- paste0('Transfer_', ix, '_predicted.id_filtered')
    for (jx in levels(scrna)) {
        curr_type <- as.data.frame(table(scrna@meta.data[WhichCells(scrna, idents = jx), curr_name]))
        colnames(curr_type) <- c('CellType', 'Count')
        curr_type <- subset(curr_type, Count != 0)
        curr_type <- subset(curr_type, CellType != ' ')
        curr_type$Cluster <- jx
        curr_type$Count <- curr_type$Count/sum(curr_type$Count)
        df_type <- rbind(df_type, curr_type)
    }
    df_type$Count <- 100*df_type$Count
    df_type <- subset(df_type, Count >= 5)
    df_type <- df_type[, c('Cluster', 'CellType', 'Count')]
    mat_type <- net_sparse_2_full(df_type)
    order_type <- rownames(mat_type)[cluster_hierarchical(t(mat_type))]
    df_type$Cluster <- factor(df_type$Cluster, levels = as.numeric(levels(scrna)))
    df_type$CellType <- factor(df_type$CellType, levels = order_type)
    ggplot(df_type, aes(Cluster, CellType, fill = Count)) + geom_tile() +
            scale_x_discrete(position = 'top') +
            theme(axis.text.x.top = element_text(vjust = 0.5)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 1, size = 14, face = 'bold')) +
            theme(axis.text.y = element_text(size = 12, face = 'bold')) +
            labs(x = '', y = '', fill = '% Label') +
            scale_fill_gradient2(low = 'white', high = 'red', oob = squish) +
            theme(panel.grid.major = element_blank(), panel.background = element_blank()) +
            geom_hline(yintercept = seq(.5, length(unique(df_type$CellType))+1, 1), lwd = 0.5, colour = 'gray30') +
            geom_vline(xintercept = seq(.5, length(unique(df_type$Cluster))+1, 1), lwd = 0.5, colour = 'gray30') +
            theme(axis.ticks = element_blank())
    file_out <- file.path(dir_out, paste0('heat_', file_save, '_Multi_clusterRes', res_clust, '_pctLabel_', ix, '.pdf'))
    ggsave(file_out, height = 7, width = 9)
}

# ---- Export cluster assignments -----------------------------------------------

# Write a two-column, tab-delimited table (Cell, Cluster) without a header.
scrna[['Cluster']] <- Idents(scrna)
meta_cluster <- scrna@meta.data[, 'Cluster', drop = FALSE]
meta_cluster$Cell <- rownames(meta_cluster)
meta_cluster <- meta_cluster[c('Cell', 'Cluster')]
file_out <- file.path(dir_out, paste0('meta_', file_save, '_Multi_', tissue_subset, '_clusterRes', res_clust, '.txt'))
write.table(meta_cluster, file_out, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

print('DONE')
