# IBD_multiome_analysis_cluster_subset.R
# module load R/4.2.3 bedtools/2.30.0 hdf5/1.8.13 igvtools/2.4.0
# Cluster an integrated multiome subset using its weighted nearest-neighbor
# graph, then visualize cluster QC metrics and cell-cycle scores.
# Optionally summarize transferred cell-type labels.
#
# Replace the /path/to placeholders before running. The input object must
# contain the wsnn graph, umap_WNN reduction, and configured QC metadata.
# CellCycleScoring uses the input object's active assay.
# tissue_subset and name_subset label outputs; they do not filter cells.
#
# Existing dependencies: use_python() must be available (from reticulate).
# Optional label heatmaps require net_sparse_2_full() and, when applicable,
# cluster_hierarchical(). These helpers are not explicitly loaded here.
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

# ---- Input paths and clustering parameters -------------------------------------

# Output directory
dir_out <- '/path/to/output/cluster_subset'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Multiome object
file_data <- '/path/to/integrate_subset/All/Stromal/obj_multiomeIBD_sub_Stromal_harmony_All.rds'
name_subset <- 'Stromal'

# Label transfer names
name_transfer <- c(
# 'Elmentaite2021',
# 'Smillie2019',
# 'Hickey2023_Epithelial_colon',
# 'Hickey2023_Epithelial_ileum',
# 'Hickey2023_Immune'
    'Hickey2023_Stromal'
# 'Elmentaite2021_Tcell'
)
# Optional label-transfer metadata directory. The NULL assignment below
# disables transfer metadata loading and plots in the default workflow.
dir_meta_transfer <- '/path/to/label_transfer_subset_metadata'
dir_meta_transfer <- NULL

# Clustering parameters
file_save <- 'multiomeIBD' # file save base
tissue_subset <- 'All' # name subset of data
res_clust <- 0.8 # clustering resolution scRNA+scATAC
algo_clust <- 4 # 1 = Louvain, 3 = SLM algorithm, 4 = Leiden

meta_plot <- c('nCount_RNA', 'nFeature_RNA', 'gex_pct_mito',
                    'atac_fragments', 'TSS.enrichment', 'nucleosome_signal', 'atac_pct_reads_in_peaks',
                    'atac_pct_mito', 'atac_pct_blacklist')

# ---- Configure Python and load the integrated subset ---------------------------

# Select the Python executable for Leiden clustering.
if (algo_clust == 4) {
    use_python('/path/to/python/environment/bin/python')
}

# Output directory
dir_out <- file.path(dir_out, tissue_subset, name_subset)
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Load multiome object
print('Load multiome object')
scrna <- readRDS(file_data)

# ---- Optional transferred labels ----------------------------------------------

# Transfer table row names must match the input object cell names.
if (!is.null(dir_meta_transfer)) {
    file_meta <- file.path(dir_meta_transfer, tissue_subset, name_subset,
                                    paste0('meta_transfer_', file_save, '_', name_subset, '_RNA.txt'))
    meta_transfer <- read.delim(file_meta, header = TRUE, row.names = 1)
    col_transfer <- paste0('Transfer_', name_transfer)
    col_transfer <- paste0(col_transfer, '_predicted.id_filtered')
    meta_transfer <- meta_transfer[, col_transfer, drop = FALSE]
    scrna <- AddMetaData(scrna, meta_transfer)
    for (ix in name_transfer) {
        curr_name <- paste0('Transfer_', ix, '_predicted.id_filtered')
        scrna@meta.data[which(is.na(scrna@meta.data[, curr_name])), curr_name] <- ' '
        scrna@meta.data[which(scrna@meta.data[, curr_name] == 'Unclassified'), curr_name] <- ' '
    }
    for (ix in name_transfer) {
        curr_name <- paste0('Transfer_', ix, '_predicted.id_filtered')
        curr_levels <- names(sort(table(scrna@meta.data[, curr_name]), decreasing = TRUE))
        curr_levels <- c(setdiff(curr_levels, ' '), ' ')
        scrna@meta.data[, curr_name] <- factor(scrna@meta.data[, curr_name], levels = curr_levels)
    }

    # Plot transferred labels and highlight individual populations.
    for (ix in name_transfer) {
        curr_name <- paste0('Transfer_', ix, '_predicted.id_filtered')
        curr_levels <- levels(scrna@meta.data[, curr_name])
        curr_order <- NULL
        Idents(scrna) <- curr_name
        curr_order <- WhichCells(scrna, idents = setdiff(levels(scrna), ' '))
        if (' ' %in% names(table(Idents(scrna)))) {
            curr_order <- c(curr_order, WhichCells(scrna, idents = ' '))
        }
        DimPlot(scrna, reduction = 'umap_WNN', group.by = curr_name, label = TRUE, repel = TRUE, pt.size = 0.5, cells = rev(curr_order),
                cols = c(hue_pal()(length(levels(scrna@meta.data[, curr_name]))-1), 'grey90')) +
                NoLegend() + labs(x = 'UMAP-1', y = 'UMAP-2')
        file_out <- file.path(dir_out, paste0('umap_', file_save, '_transfer_RNA_', ix, '.png'))
        ggsave(file_out)
        dir_out_transfer <- file.path(dir_out, paste0('Transfer_', ix))
        dir.create(dir_out_transfer, showWarnings = FALSE)
        for (jx in levels(scrna)) {
            DimPlot(scrna, reduction = 'umap_WNN', pt.size = 0.5, sizes.highlight = 0.5,
                        cells.highlight = WhichCells(scrna, idents = jx), cols.highlight = 'grey20') +
                        labs(x = 'UMAP-1', y = 'UMAP-2') + theme_void() + NoLegend()
            file_out <- file.path(dir_out_transfer, paste0('umap_', file_save, '_transfer_RNA_', ix, '_', jx, '.png'))
            ggsave(file_out)
        }
    }
}

# ---- Cluster the multiome WNN graph --------------------------------------------
print('Cluster RNA+ATAC')
scrna <- FindClusters(scrna, graph.name = 'wsnn', resolution = res_clust, algorithm = algo_clust,
                        verbose = FALSE, method = ifelse(algo_clust == 4, 'igraph', 'matrix'))
DimPlot(scrna, reduction = 'umap_WNN', label = TRUE, label.size = 6) + labs(x = 'UMAP-1', y = 'UMAP-2') + theme_void() + NoLegend()
file_out <- file.path(dir_out, paste0('umap_', file_save, '_clusterRes', res_clust, '.png'))
ggsave(file_out)

# ---- Optional label-composition heatmaps ---------------------------------------
if (!is.null(dir_meta_transfer)) {
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
        if (nrow(mat_type) > 2) {
            order_type <- rownames(mat_type)[cluster_hierarchical(t(mat_type))]
            df_type$Cluster <- factor(df_type$Cluster, levels = as.numeric(levels(scrna)))
            df_type$CellType <- factor(df_type$CellType, levels = order_type)
        }
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
        file_out <- file.path(dir_out, paste0('heat_', file_save, '_clusterRes', res_clust, '_pctLabel_', ix, '.pdf'))
        ggsave(file_out, height = 7, width = 9)
    }
}

# ---- Plot QC metrics by cluster ------------------------------------------------
Idents(scrna) <- factor(Idents(scrna), levels = rev(levels(scrna)))
VlnPlot(scrna, features = rev(meta_plot), stack = TRUE)
file_out <- file.path(dir_out, paste0('vln_', file_save, '_clusterRes', res_clust, '_qc_metrics.pdf'))
ggsave(file_out, width = 10, height = 8)
Idents(scrna) <- factor(Idents(scrna), levels = rev(levels(scrna)))

# ---- Export cluster assignments -----------------------------------------------

# Write a two-column, tab-delimited table (Cell, Cluster) without a header.
scrna[['Cluster']] <- Idents(scrna)
meta_cluster <- scrna@meta.data[, 'Cluster', drop = FALSE]
meta_cluster$Cell <- rownames(meta_cluster)
meta_cluster <- meta_cluster[c('Cell', 'Cluster')]
file_out <- file.path(dir_out, paste0('meta_', file_save, '_clusterRes', res_clust, '.txt'))
write.table(meta_cluster, file_out, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

# ---- Score and visualize cell-cycle expression --------------------------------
scrna <- CellCycleScoring(scrna, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
FeaturePlot(scrna, features = c('G2M.Score', 'S.Score'), min.cutoff = 'q5', max.cutoff = 'q95', reduction = 'umap_WNN',
            cols = c('lightgrey', 'red'), order = TRUE) + labs(x = 'UMAP-1', y = 'UMAP-2') + NoLegend()
file_out <- file.path(dir_out, paste0('umap_', file_save, '_cell_cycle.png'))
ggsave(file_out, width = 13)

print('DONE')
