# IBD_multiome_qc_filter.R
#
# Purpose: Build per-sample Multi and ATAC Seurat objects from precomputed
# QC-passing cell lists, generate diagnostic plots, and summarize cell recovery.
# RNA-passing cells are counted; this workflow exports Multi and ATAC objects.
#
# Setup: Replace all /path/to/... placeholders with local paths.
# Required sample metadata columns:
#   Sample, dir_cellranger_arc, dir_cellbender, DiseaseTissueStatus, nNuclei.
# DiseaseTissueStatus controls sample ordering (e.g., CD_R_NI); nNuclei records
# the number of nuclei loaded for the recovery summaries.
# Cell Ranger ARC directories must contain outs/filtered_feature_bc_matrix.h5
# and outs/atac_fragments.tsv.gz with its index. Each CellBender directory must
# contain <directory_basename>.h5.
# QC tables are outputs of IBD_multiome_qc_metrics.R and include Sample, barcode,
# and the per-cell QC metrics. The gene list contains one gene name per line.
# Optional removal files contain one sample-prefixed cell ID per line
# (<Sample>_<barcode>); leave file_cells_rm as NULL to skip additional removals.
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
    library(ggplot2)
    library(patchwork)
    library(EnsDb.Hsapiens.v86)
    library(dplyr)
    library(Matrix)
    library(scales)
    library(rtracklayer)
})
source('/path/to/helper_scripts/plot_utils.R')

# Inputs and analysis parameters ---------------------------------------------

# Output directory
dir_out <- '/path/to/output/qc_filter'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Cells passing QC
file_cells_both <- '/path/to/qc_metrics/meta_multiomeIBD_filtered_both.txt'
file_cells_scatac <- '/path/to/qc_metrics/meta_multiomeIBD_filtered_scatac.txt'
file_cells_scrna <- '/path/to/qc_metrics/meta_multiomeIBD_filtered_scrna.txt'

# Sample metadata
file_meta <- '/path/to/sample_metadata/meta_multiomeIBD_sample.txt'

# Gene annotation
file_annot <- '/path/to/cellranger_arc_reference/genes/genes.gtf'

# Filtered genes
file_gene <- '/path/to/gene_filter/gene_set_refdata-cellranger-arc-GRCh38-2020-A-2.0.0_filtered.txt'

# Optional additional cell removals (NULL to retain the supplied QC cell sets)
file_cells_rm <- NULL

# Analysis and plotting parameters
file_save <- 'multiomeIBD'
order_disease <- c('CD', 'UC')
order_tissue <- c('R', 'TI')
order_status <- c('NI', 'I')
genome <- 'hg38'
keep_chr <- paste0('chr', c(1:22, 'X', 'M')) # chromosomes to keep
min_atac_frag <- 2000 # minimum ATAC fragments
max_atac_frag <- 1E5 # maximum ATAC fragments
min_nCount_RNA <- 1000 # minimum RNA UMI counts
max_nCount_RNA <- 20000 # maximum RNA UMI counts
scrna_min_cells <- 20 # Min number of cells to consider gene for scRNA dimension reduction
scatac_min_cells <- 50 # Min number of cells to consider peak for scATAC dimension reduction
vars_regress <- c('nCount_RNA', 'gex_pct_mito') # Variables to regress during scaling
k_pca <- 30 # Number of PCs for scRNA dimension reduction
k_lsi <- 30 # Number of LSI components for scATAC dimension reduction
gene_upstream <- 2000 # bases upstream TSS for activity estimates
gene_downstream <- 500 # bases downstream TSS for activity estimates

# Metadata to plot
plot_meta <- c('nCount_RNA', 'nFeature_RNA', 'gex_pct_mito', 'gex_pct_ribo',
                    'atac_fragments', 'nCount_ATAC', 'nFeature_ATAC', 'atac_pct_reads_in_peaks',
                    'atac_pct_mito', 'atac_pct_blacklist', 'TSS.enrichment', 'nucleosome_signal')

# Load metadata, cell lists, and reference annotation ------------------------

# Load sample metadata
meta_sample <- read.delim(file_meta, header = TRUE, sep = '\t')
name_sample <- meta_sample$Sample

# Load cells passing QC
cells_both <- read.delim(file_cells_both, header = TRUE, sep = '\t')
cells_scatac <- read.delim(file_cells_scatac, header = TRUE, sep = '\t')
cells_scrna <- read.delim(file_cells_scrna, header = TRUE, sep = '\t')

# Load filtered genes
keep_gene <- readLines(file_gene)

# Load gene annotation
annotations <- rtracklayer::import(file_annot)
annotations$gene_biotype <- annotations$gene_type
annotations <- subset(annotations, gene_name %in% keep_gene)
annotations <- keepSeqlevels(annotations, keep_chr, pruning.mode = 'coarse')
seqlevelsStyle(annotations) <- 'UCSC'

# Load cells to remove
if (!is.null(file_cells_rm[1])) {
    cells_rm <- NULL
    for (ix in file_cells_rm) {
        cells_rm <- c(cells_rm, readLines(ix))
    }
}

# Order samples by metadata
order_sample <- NULL
for (ix in order_disease) {
    for (jx in order_tissue) {
        for (kx in order_status) {
            curr_group <- paste(ix, jx, kx, sep = '_')
            curr_sample <- sort(subset(meta_sample, DiseaseTissueStatus %in% curr_group)$Sample)
            order_sample <- c(order_sample, curr_sample)
        }
    }
}

# Process each sample --------------------------------------------------------
print('*** Process samples ***')
df_ncell <- NULL
for (ix in 1:length(name_sample)) {

    # Output directory
    dir_out_sample <- file.path(dir_out, name_sample[ix])
    dir.create(dir_out_sample, showWarnings = FALSE)

    # Cells passing QC
    curr_cells_both <- subset(cells_both, Sample == name_sample[ix])
    rownames(curr_cells_both) <- curr_cells_both$barcode
    curr_cells_scatac <- subset(cells_scatac, Sample == name_sample[ix])
    rownames(curr_cells_scatac) <- curr_cells_scatac$barcode
    curr_cells_scrna <- subset(cells_scrna, Sample == name_sample[ix])
    rownames(curr_cells_scrna) <- curr_cells_scrna$barcode

    # Remove cells
    if (!is.null(file_cells_rm[1])) {
        curr_cells_rm <- gsub(paste0(name_sample[ix], '_'), '', cells_rm[grepl(name_sample[ix], cells_rm)])
        curr_cells_both <- curr_cells_both[setdiff(rownames(curr_cells_both), curr_cells_rm), ]
        curr_cells_scatac <- curr_cells_scatac[setdiff(rownames(curr_cells_scatac), curr_cells_rm), ]
        curr_cells_scrna <- curr_cells_scrna[setdiff(rownames(curr_cells_scrna), curr_cells_rm), ]
    }

    # Load scATAC from Cell Ranger multiome output
    file_data_scatac <- paste0(meta_sample$dir_cellranger_arc[ix], '/outs/filtered_feature_bc_matrix.h5')
    sc_data_scatac <- Read10X_h5(file_data_scatac)
    atac_counts <- sc_data_scatac[['Peaks']]
    grange_counts <- StringToGRanges(rownames(atac_counts), sep = c(':', '-'))
    grange_use <- seqnames(grange_counts) %in% setdiff(standardChromosomes(grange_counts), 'chrY')
    atac_counts <- atac_counts[as.vector(grange_use), ]
    file_frag <- paste0(meta_sample$dir_cellranger_arc[ix], '/outs/atac_fragments.tsv.gz')
    chrom_assay <- CreateChromatinAssay(counts = atac_counts, sep = c(':', '-'), genome = genome,
                                            fragments = file_frag, min.cells = 5, annotation = annotations)
    keep_cells <- colnames(chrom_assay)

    # Load scRNA from CellBender output
    file_data_scrna <- paste0(meta_sample$dir_cellbender[ix], '/', basename(meta_sample$dir_cellbender[ix]), '.h5')
    sc_data_scrna <- Read10X_h5(file_data_scrna)
    sc_data_scrna <- sc_data_scrna[, keep_cells]

    # Initialize the object with ATAC-passing cells and their QC metadata.
    sc <- CreateSeuratObject(sc_data_scrna, min.cells = 1, project = name_sample[ix])
    sc[['ATAC']] <- chrom_assay
    sc <- subset(sc, cells = curr_cells_scatac$barcode)
    sc <- AddMetaData(sc, curr_cells_scatac)
    df_meta <- curr_cells_scatac

    # Build the joint RNA/ATAC subset and retain features detected in >=1 cell.
    print('- Filter cells: RNA & ATAC')
    sc_filtered_both <- subset(sc, cells = curr_cells_both$barcode)
    keep_feat_scrna <- rownames(CreateSeuratObject(LayerData(sc_filtered_both, assay = 'RNA', layer = 'counts'),
                                    min.cells = 1))
    keep_feat_scatac <- rownames(CreateSeuratObject(LayerData(sc_filtered_both, assay = 'ATAC', layer = 'counts'),
                                    min.cells = 1))
    sc_filtered_both <- subset(sc_filtered_both, features = c(keep_feat_scrna, keep_feat_scatac))
    n_cell_both <- nrow(curr_cells_both)
    n_cell_scatac <- nrow(curr_cells_scatac)
    n_cell_scrna <- nrow(curr_cells_scrna)
    print(paste0(' - quality ATAC cells: ', n_cell_scatac))
    print(paste0(' - quality RNA cells: ', n_cell_scrna))
    print(paste0(' - quality RNA & ATAC cells: ', n_cell_both))
    curr_ncell <- data.frame(Sample = rep(name_sample[ix], 3), Ncell = c(n_cell_both, n_cell_scrna, n_cell_scatac),
                                        Type = c('ATAC & RNA', 'RNA', 'ATAC'))
    df_ncell <- rbind(df_ncell, curr_ncell)

    # Save number of quality cells
    file_out <- file.path(dir_out_sample, paste0('ncell_', file_save, '_', name_sample[ix], '.txt'))
    write.table(curr_ncell, file_out, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

    # Scatter plot - QC metrics
    df_meta$Quality <- 'Low'
    df_meta[colnames(sc_filtered_both), 'Quality'] <- 'High'
    ggplot(df_meta) + geom_point(aes(x = atac_fragments, y = nCount_RNA, color = Quality), size = 0.5) +
        theme_minimal() +
        labs(x = 'ATAC: Total Fragments', y = 'RNA: Total UMIs', color = 'Quality') +
        theme(axis.title = element_text(size = 14, face = 'plain')) +
        theme(axis.text = element_text(size = 11, face = 'plain')) +
        theme(legend.text = element_text(size = 12)) +
        scale_y_continuous(trans = 'log10') + annotation_logticks(sides = 'l') +
        scale_x_continuous(trans = 'log10') + annotation_logticks(sides = 'b') +
        scale_color_manual(values = c('black', 'grey70')) +
        geom_hline(yintercept = c(min_nCount_RNA, max_nCount_RNA), linetype = 'dashed', color = 'red') +
        geom_vline(xintercept = c(min_atac_frag, max_atac_frag), linetype = 'dashed', color = 'red')
    file_out <- file.path(dir_out_sample, paste0('qc_', file_save, '_raw_scatter_counts_quality_', name_sample[ix], '.pdf'))
    ggsave(file_out, height = 4, width = 5.2)

    # Variable features
    print('- Variable features')
    keep_var_scrna <- rownames(CreateSeuratObject(LayerData(sc_filtered_both, assay = 'RNA', layer = 'counts'),
                                    min.cells = scrna_min_cells))
    keep_var_scrna <- intersect(keep_var_scrna, keep_gene)
    keep_var_scatac <- rownames(CreateSeuratObject(LayerData(sc_filtered_both, assay = 'ATAC', layer = 'counts'),
                                    min.cells = scatac_min_cells))
    keep_var <- c(keep_var_scrna, keep_var_scatac)
    feat_var_scrna <- intersect(VariableFeatures(FindVariableFeatures(sc_filtered_both, nfeatures = 1E4, verbose = FALSE)),
                                    keep_var)[1:2000]
    VariableFeatures(sc_filtered_both) <- feat_var_scrna
    DefaultAssay(sc_filtered_both) <- 'ATAC'
    feat_var_scatac <- intersect(VariableFeatures(FindTopFeatures(sc_filtered_both, min.cutoff = 'q0', verbose = FALSE)),
                                    keep_var)
    VariableFeatures(sc_filtered_both) <- feat_var_scatac
    DefaultAssay(sc_filtered_both) <- 'RNA'

    # Seurat workflow - RNA
    print('- Seurat workflow: RNA')
    sc_filtered_both <- NormalizeData(sc_filtered_both, assay = 'RNA', verbose = FALSE)
    sc_filtered_both <- ScaleData(sc_filtered_both, vars.to.regress = vars_regress, verbose = FALSE)
    sc_filtered_both <- RunPCA(sc_filtered_both, verbose = FALSE)
    sc_filtered_both <- RunUMAP(sc_filtered_both, dims = 1:k_pca, verbose = FALSE, reduction.name = 'umap_RNA')

    # UMAP - RNA
    DimPlot(sc_filtered_both, pt.size = 1.5, reduction = 'umap_RNA') + labs(x = 'UMAP=1', y = 'UMAP-2') +
                theme_void() + NoLegend()
    file_out <- file.path(dir_out_sample, paste0('umap_', file_save, '_filterMulti_scrna_', name_sample[ix], '.png'))
    ggsave(file_out)
    for (jx in plot_meta) {
        FeaturePlot(sc_filtered_both, features = jx, pt.size = 1.5, reduction = 'umap_RNA', order = TRUE,
                        min.cutoff = 'q5', max.cutoff = 'q95', cols = c('lightgrey', 'red')) + labs(x = 'UMAP=1', y = 'UMAP-2')
        file_out <- file.path(dir_out_sample, paste0('umap_', file_save, '_filterMulti_scrna_',
                                    name_sample[ix], '_meta_', jx, '.png'))
        ggsave(file_out)
    }

    # Seurat workflow - ATAC
    print('- Seurat workflow: ATAC')
    DefaultAssay(sc_filtered_both) <- 'ATAC'
    sc_filtered_both <- RunTFIDF(sc_filtered_both, verbose = FALSE)
    sc_filtered_both <- RunSVD(sc_filtered_both, assay = 'ATAC', reduction.key = 'LSI_',
                                    reduction.name = 'lsi', verbose = FALSE)
    cor_depth <- DepthCor(sc_filtered_both, assay = 'ATAC', reduction = 'lsi', n = k_lsi)
    file_out <- file.path(dir_out_sample, paste0('lsi_', file_save, '_filterMulti_scatac_',
                            name_sample[ix], '_cor_depth.pdf'))
    ggsave(file_out, cor_depth)
    keep_dim_lsi <- which(abs(cor_depth$data$nCount_ATAC) < 0.75)
    sc_filtered_both <- RunUMAP(sc_filtered_both, reduction = 'lsi', dims = keep_dim_lsi,
                                    reduction.name = 'umap_ATAC', verbose = FALSE)
    # UMAP - scATAC
    DimPlot(sc_filtered_both, pt.size = 1.5, reduction = 'umap_ATAC', group.by = 'orig.ident') +
                labs(x = 'UMAP=1', y = 'UMAP-2') + theme_void() + NoLegend() + ggtitle('')
    file_out <- file.path(dir_out_sample, paste0('umap_', file_save, '_filterMulti_scatac_', name_sample[ix], '.png'))
    ggsave(file_out)
    for (jx in plot_meta) {
        FeaturePlot(sc_filtered_both, features = jx, pt.size = 1.5, reduction = 'umap_ATAC', order = TRUE,
                        min.cutoff = 'q5', max.cutoff = 'q95', cols = c('lightgrey', 'red')) + labs(x = 'UMAP=1', y = 'UMAP-2')
        file_out <- file.path(dir_out_sample, paste0('umap_', file_save, '_filterMulti_scatac_',
                                    name_sample[ix], '_meta_', jx, '.png'))
        ggsave(file_out)
    }


    # Filter cells - ATAC only
    # Include all ATAC-passing cells, including cells that also pass RNA QC.
    print('- Filter cells: ATAC only')
    sc_atac <- subset(sc, cells = curr_cells_scatac$barcode)
    DefaultAssay(sc_atac) <- 'ATAC'
    sc_atac[['RNA']] <- NULL
    keep_feat_scatac <- rownames(CreateSeuratObject(LayerData(sc_atac, assay = 'ATAC', layer = 'counts'), min.cells = 1))
    sc_atac <- subset(sc_atac, features = keep_feat_scatac)

    # Scatter plot QC metrics
    df_meta$Quality_ATAC <- 'Low'
    df_meta[colnames(sc_atac), 'Quality_ATAC'] <- 'High'
    ggplot(df_meta) + geom_point(aes(x = atac_fragments, y = nCount_RNA, color = Quality_ATAC), size = 0.5) +
        theme_minimal() +
        labs(x = 'ATAC Cut Sites in Peaks', y = 'RNA UMIs', color = 'Quality') +
        theme(axis.title = element_text(size = 14, face = 'plain')) +
        theme(axis.text = element_text(size = 11, face = 'plain')) +
        theme(legend.text = element_text(size = 12)) +
        scale_y_continuous(trans = 'log10') + annotation_logticks(sides = 'l') +
        scale_x_continuous(trans = 'log10') + annotation_logticks(sides = 'b') +
        scale_color_manual(values = c('black', 'grey70')) +
        geom_vline(xintercept = c(min_atac_frag, max_atac_frag), linetype = 'dashed', color = 'red')
    file_out <- file.path(dir_out_sample, paste0('qc_', file_save, '_raw_scatter_counts_quality_',
                                                        name_sample[ix], '_filterATAC.pdf'))
    ggsave(file_out, height = 4, width = 5.2)

    # Variable features
    print('- Variable features')
    keep_var_scatac <- rownames(CreateSeuratObject(LayerData(sc_atac, assay = 'ATAC', layer = 'counts'),
                                    min.cells = scatac_min_cells))
    keep_var <- keep_var_scatac
    feat_var_scatac <- intersect(VariableFeatures(FindTopFeatures(sc_atac, min.cutoff = 'q0', verbose = FALSE)), keep_var)
    VariableFeatures(sc_atac) <- feat_var_scatac

    # Seurat workflow - ATAC
    print('- Seurat workflow: ATAC')
    sc_atac <- RunTFIDF(sc_atac, verbose = FALSE)
    sc_atac <- FindTopFeatures(sc_atac, min.cutoff = 'q0', verbose = FALSE)
    sc_atac <- RunSVD(sc_atac, assay = 'ATAC', n = k_lsi, reduction.key = 'LSI_', reduction.name = 'lsi', verbose = FALSE)
    cor_depth <- DepthCor(sc_atac, assay = 'ATAC', reduction = 'lsi', n = k_lsi)
    file_out <- file.path(dir_out_sample, paste0('lsi_', file_save, '_filterATAC_scatac_',
                                name_sample[ix], '_cor_depth.pdf'))
    ggsave(file_out, cor_depth)
    keep_dim_lsi <- which(abs(cor_depth$data$nCount_ATAC) < 0.75)
    sc_atac <- RunUMAP(sc_atac, reduction = 'lsi', dims = keep_dim_lsi, reduction.name = 'umap_ATAC', verbose = FALSE)
    DimPlot(sc_atac, pt.size = 1.5, reduction = 'umap_ATAC', group.by = 'orig.ident') + labs(x = 'UMAP=1', y = 'UMAP-2') +
                theme_void() + NoLegend() + ggtitle('')
    file_out <- file.path(dir_out_sample, paste0('umap_', file_save, '_filterATAC_scatac_', name_sample[ix], '.png'))
    ggsave(file_out)
    for (jx in plot_meta) {
        FeaturePlot(sc_atac, features = jx, pt.size = 1.5, reduction = 'umap_ATAC', order = TRUE) + labs(x = 'UMAP=1', y = 'UMAP-2')
        file_out <- file.path(dir_out_sample, paste0('umap_', file_save, '_filterATAC_scatac_',
                                    name_sample[ix], '_meta_', jx, '.png'))
        ggsave(file_out)
    }
    DimPlot(sc_atac, pt.size = 1.5, reduction = 'umap_ATAC', cells.highlight = colnames(sc_filtered_both)) +
                labs(x = 'UMAP=1', y = 'UMAP-2') + theme_void() + NoLegend() + ggtitle('')
    file_out <- file.path(dir_out_sample, paste0('umap_', file_save, '_filterATAC_scatac_',
                                name_sample[ix], '_RNAhighQuality.png'))
    ggsave(file_out)
    DimPlot(sc_atac, pt.size = 1.5, reduction = 'umap_ATAC',
                cells.highlight = setdiff(colnames(sc_atac), colnames(sc_filtered_both))) +
                labs(x = 'UMAP=1', y = 'UMAP-2') + theme_void() + NoLegend() + ggtitle('')
    file_out <- file.path(dir_out_sample, paste0('umap_', file_save, '_filterATAC_scatac_',
                                name_sample[ix], '_RNAlowQuality.png'))
    ggsave(file_out)

    # Estimate gene activities
    print('- Estimate gene activities')
    gene_activities <- GeneActivity(sc_atac, extend.upstream = gene_upstream,
                                            extend.downstream = gene_downstream, verbose = FALSE)
    gene_activities <- gene_activities[which(!duplicated(rownames(gene_activities))), ]
    gene_activities <- gene_activities[intersect(rownames(gene_activities), keep_gene), ]
    sc_atac[['ATAC_Activity']] <- CreateAssayObject(gene_activities, min.cells = 1)
    sc_atac <- NormalizeData(sc_atac, assay = 'ATAC_Activity', scale.factor = median(sc_atac$nCount_ATAC_Activity),
                                        normalization.method = 'LogNormalize', verbose = FALSE)

    # Add quality cell information to objects
    df_quality <- data.frame(row.names = colnames(sc_atac), good_rna = rep(0, ncol(sc_atac)),
                                                            good_atac = rep(1, ncol(sc_atac)))
    df_quality[colnames(sc_filtered_both), 'good_rna'] <- 1
    sc_filtered_both <- AddMetaData(sc_filtered_both, df_quality)
    sc_atac <- AddMetaData(sc_atac, df_quality)

    # Export the joint RNA/ATAC object.
    print('- Save filtered object: Multi')
    file_out <- file.path(dir_out_sample, paste0('obj_', file_save, '_', name_sample[ix], '_Multi.rds'))
    saveRDS(sc_filtered_both, file_out)

    # Save filtered ATAC object
    print('- Save filtered ATAC object')
    file_out <- file.path(dir_out_sample, paste0('obj_', file_save, '_', name_sample[ix], '_ATAC.rds'))
    saveRDS(sc_atac, file_out)

}

# Export combined cell counts ------------------------------------------------
file_out <- file.path(dir_out, paste0('ncell_', file_save, '.txt'))
write.table(df_ncell, file_out, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

# Plot cell counts and recovery ----------------------------------------------
df_ncell$Sample <- factor(df_ncell$Sample, levels = order_sample)
df_ncell$Type <- factor(df_ncell$Type, levels = c('ATAC', 'RNA', 'ATAC & RNA'))
ggplot(df_ncell, aes(x = Sample, y = Ncell, fill = Type)) +
            geom_bar(stat = 'identity', position = position_dodge(), width = 0.8) +
            theme_minimal() +
            theme(axis.text.x = element_text(size = 10, face = 'bold', angle = 90, hjust = 1, vjust = 0.5)) +
            theme(axis.text.y = element_text(size = 13, face = 'plain', color = 'black')) +
            theme(axis.title = element_text(size = 14, face = 'plain')) +
            labs(x = '', y = 'Nuclei', fill = 'Filter') +
            theme(legend.text = element_text(size = 12, face = 'bold')) +
            theme(legend.title = element_text(size = 12, face = 'bold')) +
            scale_fill_manual(values = c('steelblue2', 'gray20', 'salmon3'))
file_out <- file.path(dir_out, paste0('ncell_', file_save, '_bar.pdf'))
ggsave(file_out, height = 5, width = 10)

add_ncell <- data.frame(Sample = meta_sample$Sample, Ncell = meta_sample$nNuclei, Type = 'Loaded')
df_ncell2 <- rbind(df_ncell, add_ncell)
df_ncell2$Sample <- factor(df_ncell2$Sample, levels = order_sample)
df_ncell2$Type <- factor(df_ncell2$Type, levels = c('Loaded', 'ATAC', 'RNA', 'ATAC & RNA'))
ggplot(df_ncell2, aes(x = Sample, y = Ncell, fill = Type)) +
            geom_bar(stat = 'identity', position = position_dodge(), width = 0.8) +
            theme_minimal() +
            theme(axis.text.x = element_text(size = 10, face = 'bold', angle = 90, hjust = 1, vjust = 0.5)) +
            theme(axis.text.y = element_text(size = 13, face = 'plain', color = 'black')) +
            theme(axis.title = element_text(size = 14, face = 'plain')) +
            labs(x = '', y = 'Nuclei', fill = 'Filter') +
            theme(legend.text = element_text(size = 12, face = 'bold')) +
            theme(legend.title = element_text(size = 12, face = 'bold')) +
            scale_fill_manual(values = c('grey20', 'tomato2', 'steelblue2', 'plum'))
file_out <- file.path(dir_out, paste0('ncell_', file_save, '_bar2.pdf'))
ggsave(file_out, height = 5, width = 10)

# Quality cells as a percent of cells loaded
df_frac <- NULL
for (ix in levels(df_ncell2$Sample)) {
    curr_frac <- subset(df_ncell2, Sample == ix & Type !='Loaded')
    curr_load <- subset(df_ncell2, Sample == ix & Type == 'Loaded')$Ncell
    curr_frac$Frac <- curr_frac$Ncell/curr_load*100
    df_frac <- rbind(df_frac, curr_frac)
}
ggplot(df_frac, aes(x = Sample, y = Frac, fill = Type)) +
            geom_bar(stat = 'identity', position = position_dodge(), width = 0.8) +
            theme_minimal() +
            theme(axis.text.x = element_text(size = 10, face = 'bold', angle = 90, hjust = 1, vjust = 0.5)) +
            theme(axis.text.y = element_text(size = 13, face = 'plain', color = 'black')) +
            theme(axis.title = element_text(size = 14, face = 'plain')) +
            labs(x = '', y = '% High Quality', fill = 'Filter') +
            theme(legend.text = element_text(size = 12, face = 'bold')) +
            theme(legend.title = element_text(size = 12, face = 'bold')) +
            scale_fill_manual(values = c('tomato2', 'steelblue2', 'plum'))
file_out <- file.path(dir_out, paste0('ncell_', file_save, '_bar3.pdf'))
ggsave(file_out, height = 5, width = 10)

print('DONE')
