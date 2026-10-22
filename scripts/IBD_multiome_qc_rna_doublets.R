# IBD_multiome_qc_rna_doublets.R
#
# Purpose: Identify putative RNA doublets per sample using DoubletFinder.
# Export per-sample and combined doublet lists, UMAPs, and QC plots.
#
# Setup: Replace each /path/to/... placeholder with the appropriate local path.
# The tab-delimited sample metadata must contain:
#   Sample             Sample identifier used in output directories and filenames.
#   dir_cellranger_arc Cell Ranger ARC sample directory containing
#                      outs/filtered_feature_bc_matrix.h5.
#   dir_cellbender     CellBender directory containing an HDF5 file named after
#                      the directory: <directory_basename>.h5.
# The filtered gene list must contain one gene name per line.
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
    library(DoubletFinder)
    library(ggplot2)
})

# Inputs and analysis parameters ---------------------------------------------

# Output directory
dir_out <- '/path/to/output/rna_doublets'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Sample metadata
file_meta <- '/path/to/sample_metadata/meta_multiomeIBD_sample.txt'

# Filtered genes
file_gene <- '/path/to/gene_filter/gene_set_refdata-cellranger-arc-GRCh38-2020-A-2.0.0_filtered.txt'

# Analysis parameters
file_save <- 'multiomeIBD' # output file name
pattern_pct_mito <- '^MT-' # gene pattern for percent mitochondrial expression
min_nCount_RNA <- 1000 # minimum RNA UMI counts
max_nCount_RNA <- 20000 # maximum RNA UMI counts
min_nFeature_RNA <- 200 # minimum RNA genes
max_nFeature_RNA <- 8000 # maximum RNA genes
max_gex_pct_mito <- 25 # maximum RNA percent mitochondrial gene
min_cells <- 20 # Min number of cells to consider gene for scRNA dimension reduction
vars_regress <- c('nCount_RNA', 'gex_pct_mito') # Variables to regress during scaling
k_pca <- 30 # Number of PCs for scRNA dimension reduction
clust_res <- 0.5 # Clustering resolution

# Load metadata and gene list ------------------------------------------------

# Load sample metadata
meta_sample <- read.delim(file_meta, header = TRUE, sep = '\t')
name_sample <- meta_sample$Sample

# Load filtered genes
keep_gene <- readLines(file_gene)

# Process each sample --------------------------------------------------------
cells_all_doublets <- NULL
for (ix in 1:length(name_sample)) {

    # Sample directory
    print(paste0('Doublet removal: ', name_sample[ix]))
    dir_out_sample <- file.path(dir_out, name_sample[ix])
    dir.create(dir_out_sample, showWarnings = FALSE)

    # Define eligible barcodes from the filtered Cell Ranger ARC matrix.
    file_data_scrna <- paste0(meta_sample$dir_cellranger_arc[ix], '/outs/filtered_feature_bc_matrix.h5')
    data_scrna_tmp <- Read10X_h5(file_data_scrna)
    valid_barcodes <- colnames(data_scrna_tmp[['Gene Expression']])

    # Estimate the doublet rate using the original analysis assumption:
    # 0.8% per 1,000 barcodes in the filtered Cell Ranger ARC matrix.
    rate_doublet <- 0.008 * length(valid_barcodes)/1000

    # Load scRNA from CellBender output
    file_data_scrna <- paste0(meta_sample$dir_cellbender[ix], '/', basename(meta_sample$dir_cellbender[ix]), '.h5')
    data_scrna <- Read10X_h5(file_data_scrna)
    scrna <- CreateSeuratObject(data_scrna, min.cells = 1)

    # Filter cells
    scrna <- subset(scrna, cells = valid_barcodes)
    scrna[['gex_pct_mito']] <- PercentageFeatureSet(scrna, pattern = pattern_pct_mito)
    scrna <- subset(scrna, subset = nCount_RNA >= min_nCount_RNA & nCount_RNA <= max_nCount_RNA &
                                    nFeature_RNA >= min_nFeature_RNA & nFeature_RNA <= max_nFeature_RNA &
                                    gex_pct_mito <= max_gex_pct_mito)

    # Filter genes
    keep_feat_scrna <- rownames(CreateSeuratObject(LayerData(scrna, assay = 'RNA', layer = 'counts'),
                                    min.cells = min_cells))
    keep_feat_scrna <- intersect(keep_feat_scrna, keep_gene)
    scrna <- subset(scrna, features = keep_feat_scrna)

    # Normalize, scale, and cluster using the original Seurat workflow order.
    scrna <- FindVariableFeatures(scrna, verbose = FALSE)
    scrna <- NormalizeData(scrna, assay = 'RNA', verbose = FALSE)
    scrna <- ScaleData(scrna, vars.to.regress = vars_regress, verbose = FALSE)
    scrna <- RunPCA(scrna, verbose = FALSE)
    scrna <- RunUMAP(scrna, dims = 1:k_pca, reduction.name = 'umap_RNA', verbose = FALSE)
    scrna <- FindNeighbors(scrna, dims = 1:k_pca, reduction = 'pca', assay = 'RNA', verbose = FALSE)
    scrna <- FindClusters(scrna, graph.name = 'RNA_snn', resolution = clust_res, verbose = FALSE)
    DimPlot(scrna, reduction = 'umap_RNA', label = TRUE, label.size = 6) + NoLegend() + theme_void()
    file_out <- file.path(dir_out_sample, paste0('umap_', file_save, '_', name_sample[ix], '_clust', clust_res, '.png'))
    ggsave(file_out)

    # Select pK and classify doublets with DoubletFinder.
    sweep_res <- paramSweep(scrna, PCs = 1:k_pca, sct = FALSE, num.cores = 8)
    sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
    bcmvn <- find.pK(sweep_stats)
    homotypic_prop <- modelHomotypic(Idents(scrna))
    nExp_poi <- round(rate_doublet * length(colnames(scrna)))
    # Retain the homotypic-adjusted estimate; classification below uses nExp_poi.
    nExp_poi_adj <- round(nExp_poi * (1-homotypic_prop))
    bcmvn$pK_num <- as.numeric(levels(bcmvn$pK))
    opt_pK <- as.numeric(bcmvn$pK_num[which(bcmvn$BCmetric == max(bcmvn$BCmetric))])
    ggplot(bcmvn, aes(x = pK_num, y = BCmetric)) + geom_point() +
            geom_vline(xintercept = opt_pK, linetype = 'dashed', color = 'red') +
            labs(x = 'pK', y = 'Bimodality Coefficient (BCmvn)')
    file_out <- file.path(dir_out_sample, paste0('DoubletFinder_', name_sample[ix], '_pK_vs_BCmetric.pdf'))
    ggsave(file_out)
    scrna <- doubletFinder(scrna, PCs = 1:k_pca, pN = 0.25, pK = opt_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    metaname_doublet <- colnames(scrna@meta.data)[which(grepl('DF.classifications', colnames(scrna@meta.data)))]
    doublet_id_DoubletFinder <- scrna@meta.data[, metaname_doublet]
    names(doublet_id_DoubletFinder) <- rownames(scrna@meta.data)
    scrna <- AddMetaData(scrna, doublet_id_DoubletFinder, col.name = 'doublet_id_DoubletFinder')

    # Doublet UMAP
    DimPlot(scrna, group.by = 'doublet_id_DoubletFinder') + labs(x = 'UMAP-1', y = 'UMAP-2')
    file_out <- file.path(dir_out_sample, paste0('umap_DoubletFinder_pK', gsub('[:.:]', 'p', opt_pK), '.pdf'))
    ggsave(file_out)

    # Save doublets
    curr_doublets <- data.frame(Sample = name_sample[ix], Cell = names(which(doublet_id_DoubletFinder == 'Doublet')))
    cells_all_doublets <- rbind(cells_all_doublets, curr_doublets)
    file_out <- file.path(dir_out_sample, paste0('cells_', file_save, '_qc_scrna_', name_sample[ix], '_DoubletFinder.txt'))
    write.table(curr_doublets, file_out, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

    # Plot quality metrics
    curr_df <- scrna@meta.data[, c('nCount_RNA', 'nFeature_RNA', 'gex_pct_mito'), drop = FALSE]
    curr_df$Quality <- 'Singlet'
    curr_df$Sample <- name_sample[ix]
    curr_df[curr_doublets$Cell, 'Quality'] <- 'Doublet'
    curr_df$Quality <- factor(curr_df$Quality, levels = c('Singlet', 'Doublet'))
    p1 <- ggplot(curr_df, aes(x = Sample, y = nCount_RNA, fill = Quality)) +
            geom_violin(width = 0.6) + geom_boxplot(width = 0.6)  +
            scale_fill_manual(values = c('grey40', 'grey70')) +
            labs(x = '', y = 'Total UMIs') +
            theme(axis.title = element_text(size = 15, face = 'plain')) +
            theme(axis.text = element_text(size = 12, face = 'plain')) +
            theme(legend.position = 'none')
    p2 <- ggplot(curr_df, aes(x = Sample, y = nFeature_RNA, fill = Quality)) +
            geom_violin(width = 0.6) + geom_boxplot(width = 0.6)  +
            scale_fill_manual(values = c('grey40', 'grey70')) +
            labs(x = '', y = 'Genes Detected') +
            theme(axis.title = element_text(size = 15, face = 'plain')) +
            theme(axis.text = element_text(size = 12, face = 'plain')) +
            theme(legend.position = 'none')
    p3 <- ggplot(curr_df, aes(x = Sample, y = gex_pct_mito, fill = Quality)) +
            geom_violin(width = 0.6) + geom_boxplot(width = 0.6)  +
            scale_fill_manual(values = c('grey40', 'grey70')) +
            labs(x = '', y = '% Mitochondrial') +
            theme(axis.title = element_text(size = 15, face = 'plain')) +
            theme(axis.text = element_text(size = 12, face = 'plain'))
    p <- p1 + p2 + p3
    file_out <- file.path(dir_out_sample, paste0('vln_DoubletFinder_pK', gsub('[:.:]', 'p', opt_pK), '_qc.pdf'))
    ggsave(file_out, height = 5, width = 7)

}

# Export combined doublet list -----------------------------------------------
file_out <- file.path(dir_out, paste0('cells_', file_save, '_qc_scrna_DoubletFinder.txt'))
write.table(cells_all_doublets, file_out, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

print('DONE')
