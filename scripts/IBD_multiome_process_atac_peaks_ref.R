# IBD_multiome_process_atac_peaks_ref.R
#
# Purpose: Combine sample-specific peaks into a common reference, count processed
# ATAC cut sites in that reference, and update each sample's ATAC Seurat object.
# Export reference peaks, updated objects, and diagnostic plots.
#
# Setup: Replace all /path/to/... placeholders with local paths.
# Required sample metadata columns: Sample and DiseaseTissueStatus.
# Sample peak files: <dir_peaks>/<Sample>/peaks_<file_save>_<Sample>.rds.
# Filtered objects: <dir_filtered>/<Sample>/<file_filter_prefix>_<Sample>_ATAC.rds.
# Processed cut-site files: <dir_fragments>/fragments_<Sample>.tsv.gz, with indexes.
# Input objects must contain the QC metadata listed in plot_meta.
# The filtered gene list contains one gene name per line.
#
# Existing dependency: the plotting section calls gg_vln_bar() from plot_utils.R,
# but this script does not source that helper. This dependency is documented
# without changing the original execution workflow.
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
    library(GenomicRanges)
    library(ggplot2)
    library(magrittr)
    library(sva)
    library(Matrix)
    library(rtracklayer)
    library(scales)
    library(DESeq2)
})
source('/path/to/helper_scripts/plot_utils.R')

# Inputs and analysis parameters ---------------------------------------------

# Output directory
dir_out <- '/path/to/output/atac_peaks_ref'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Sample metadata
file_meta <- '/path/to/sample_metadata/meta_multiomeIBD_sample.txt'

# Directory of sample-specific MACS2 peaks
dir_peaks <- '/path/to/atac_peaks_sample/sizeFrag200_FDR0p1_lambda_llocal10000'

# Directory of sample-specific ATAC objects
dir_filtered <- '/path/to/qc_filter'

# Fragment files directory
dir_fragments <- '/path/to/processed_atac_fragments'

# Gene annotation
file_annot <- '/path/to/cellranger_arc_reference/genes/genes.gtf'

# Filtered genes
file_gene <- '/path/to/gene_filter/gene_set_refdata-cellranger-arc-GRCh38-2020-A-2.0.0_filtered.txt'

# Analysis parameters
file_save <- 'multiomeIBD' # file save base
file_filter_prefix <- 'obj_multiomeIBD'
name_peaks <- 'sizeFrag200_FDR0p1_lambda_llocal10000' # reference peak name for output files
keep_chr <- paste0('chr', c(1:22, 'X', 'M')) # chromosomes to keep
order_disease <- c('CD', 'UC')
order_tissue <- c('R', 'TI')
order_status <- c('NI', 'I')
gene_upstream <- 2000 # bases upstream gene for activity estimates
gene_downstream <- 500 # bases downstream gene for activity estimates
min_cells <- 100 # minimum number of cells for feature to be considered in dimension reduction
k_lsi <- 30 # dimension reduction individual samples
cutoff_cor_lsi <- 0.75 # cutoff correlation with count depth to exclude LSI components

# Metadata to plot
plot_meta <- c('nCount_RNA', 'nFeature_RNA', 'gex_pct_mito', 'gex_pct_ribo',
                    'atac_fragments', 'nCount_ATAC', 'nFeature_ATAC', 'atac_pct_reads_in_peaks',
                    'atac_pct_mito', 'atac_pct_blacklist', 'TSS.enrichment', 'nucleosome_signal')

# Configure outputs and load reference annotation ----------------------------

# Output directory
dir_out <- file.path(dir_out, name_peaks)
dir.create(dir_out, showWarnings = FALSE)

# Load sample metadata
meta_sample <- read.delim(file_meta, header = TRUE, sep = '\t')
name_sample <- meta_sample$Sample
rownames(meta_sample) <- meta_sample$Sample

# Load filtered genes
keep_gene <- readLines(file_gene)

# Load gene annotation
annotations <- rtracklayer::import(file_annot)
annotations$gene_biotype <- annotations$gene_type
annotations <- subset(annotations, gene_name %in% keep_gene)
annotations <- keepSeqlevels(annotations, keep_chr, pruning.mode = 'coarse')
seqlevelsStyle(annotations) <- 'UCSC'

# Order samples by metadata
order_sample <- NULL
order_group <- NULL
for (ix in order_disease) {
    for (jx in order_tissue) {
        for (kx in order_status) {
            curr_group <- paste(ix, jx, kx, sep = '_')
            curr_sample <- sort(subset(meta_sample, DiseaseTissueStatus %in% curr_group)$Sample)
            order_sample <- c(order_sample, curr_sample)
            order_group <- c(order_group, curr_group)
        }
    }
}
order_group <- intersect(order_group, meta_sample$DiseaseTissueStatus)

# Summarize sample-specific peak counts and sizes -----------------------------
print('Plot peak distributions')
meta_peaks_size <- NULL
meta_peaks_n <- NULL
for (ix in 1:length(name_sample)) {
    file_peaks <- file.path(dir_peaks, name_sample[ix], paste0('peaks_', file_save, '_', name_sample[ix], '.rds'))
    curr_peaks <- as.data.frame(readRDS(file_peaks))[, 1:3]
    curr_size <- data.frame(Size = abs(curr_peaks[, 3]-curr_peaks[, 2]), Sample = name_sample[ix])
    meta_peaks_size <- rbind(meta_peaks_size, curr_size)
    curr_n <- data.frame(N = nrow(curr_peaks), Sample = name_sample[ix])
    meta_peaks_n <- rbind(meta_peaks_n, curr_n)
}
meta_peaks_size$Sample <- factor(meta_peaks_size$Sample, levels = order_sample)
meta_peaks_n$Sample <- factor(meta_peaks_n$Sample, levels = order_sample)
meta_peaks_n$DiseaseTissueStatus <- meta_sample[as.character(meta_peaks_n$Sample), 'DiseaseTissueStatus']
meta_peaks_n$DiseaseTissueStatus <- factor(meta_peaks_n$DiseaseTissueStatus, levels = order_group)
# Peak size
p <- gg_vln_bar(meta_peaks_size, x = 'Sample', y = 'Size', fill = 'Sample', type = 'vln',
                    flag_box = TRUE, box_width = 0.2,
                    flag_leg = FALSE, lab_y = 'Peak Length [bp]',
                    fill_col = rep('white', length(name_sample)), angle_x = 90, size_text = 10) +
                    theme(axis.text.x = element_text(vjust = 0.5))
file_out <- file.path(dir_out, paste0('peaks_vln_', name_peaks, '_size.pdf'))
ggsave(file_out, plot = p, width = 7, height = 5)
p <- gg_vln_bar(meta_peaks_size, x = 'Sample', y = 'Size', fill = 'Sample', type = 'vln',
                    flag_box = TRUE, box_width = 0.2,
                    flag_leg = FALSE, lab_y = 'Peak Length [bp]', flag_log10 = TRUE,
                    fill_col = rep('white', length(name_sample)), angle_x = 90, size_text = 10) +
                    theme(axis.text.x = element_text(vjust = 0.5))
file_out <- file.path(dir_out, paste0('peaks_vln_', name_peaks, '_size_log10.pdf'))
ggsave(file_out, plot = p, width = 7, height = 5)
# Number of peaks per sample
ggplot(meta_peaks_n, aes(x = Sample, y = N, fill = DiseaseTissueStatus)) +
            geom_bar(stat = 'identity', position = position_dodge(), width = 0.8) +
            theme_minimal() +
            theme(axis.text.x = element_text(size = 10, face = 'bold', angle = 90, hjust = 1, vjust = 0.5)) +
            theme(axis.text.y = element_text(size = 13, face = 'plain', color = 'black')) +
            theme(axis.title = element_text(size = 14, face = 'plain')) +
            labs(x = '', y = 'Peaks', fill = '') +
            theme(legend.text = element_text(size = 12, face = 'bold')) +
            theme(legend.title = element_text(size = 12, face = 'bold'))
file_out <- file.path(dir_out, paste0('peaks_bar_', name_peaks, '_Npeaks.pdf'))
ggsave(file_out, height = 5, width = 8)

# Combine and reduce sample-specific peaks -----------------------------------
print('Reference peaks')
peaks_ref <- c()
for (ix in 1:length(name_sample)) {
    file_peaks <- file.path(dir_peaks, name_sample[ix], paste0('peaks_', file_save, '_', name_sample[ix], '.rds'))
    curr_peaks <- readRDS(file_peaks)
    peaks_ref <- append(peaks_ref, curr_peaks)
}
peaks_ref <- GenomicRanges::reduce(peaks_ref)
peaks_ref_bed <- as.data.frame(peaks_ref)[, 1:3]

# Export reference peaks using the original coordinate handling.
file_out <- file.path(dir_out, paste0('peaks_ref_', file_save, '_', name_peaks, '.rds'))
saveRDS(peaks_ref, file_out)
file_out <- file.path(dir_out, paste0('peaks_ref_', file_save, '_', name_peaks, '.bed'))
write.table(peaks_ref_bed, file_out, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

# Map each sample to the reference and update its ATAC object -----------------
print('Map fragments to reference peaks')
for (ix in 1:length(name_sample)) {

    # Load object
    print(paste0('Sample ', ix, '/', length(name_sample), ': ', name_sample[ix]))
    dir_out_sample <- file.path(dir_out, name_sample[ix])
    dir.create(dir_out_sample, showWarnings = FALSE)
    file_data <- file.path(dir_filtered, name_sample[ix], paste0(file_filter_prefix, '_', name_sample[ix], '_ATAC.rds'))
    curr_scatac <- readRDS(file_data)
    DefaultAssay(curr_scatac) <- 'ATAC'

    # Replace the fragment reference with the processed cut-site file.
    Fragments(curr_scatac@assays$ATAC) <- NULL
    file_scatac_frag <- file.path(dir_fragments, paste0('fragments_', name_sample[ix], '.tsv.gz'))
    new_frag <- CreateFragmentObject(path = file_scatac_frag, cells = colnames(curr_scatac), validate.fragments = TRUE)
    Fragments(curr_scatac@assays$ATAC) <- new_frag

    # Peak counts matrix with reference peaks
    print('- Peak count matrix')
    mat <- FeatureMatrix(fragments = Fragments(curr_scatac), features = peaks_ref,
                            cells = colnames(curr_scatac), process_n = 10000, verbose = FALSE)

    # Create object with reference peaks
    new_assay <- CreateChromatinAssay(mat, fragments = Fragments(curr_scatac), annotation = annotations)
    curr_scatac[['ATAC']] <- new_assay

    # Recompute TF-IDF, LSI, and UMAP using the reference-peak counts.
    print('- Dimension reduction')
    curr_scatac <- RunTFIDF(curr_scatac, verbose = FALSE)
    curr_scatac <- FindTopFeatures(curr_scatac, min.cutoff = min_cells, verbose = FALSE)
    curr_scatac <- RunSVD(curr_scatac, assay = 'ATAC', n = k_lsi, reduction.key = 'LSI_', reduction.name = 'lsi', verbose = FALSE)
    cor_depth <- DepthCor(curr_scatac, assay = 'ATAC', reduction = 'lsi', n = k_lsi)
    file_out <- file.path(dir_out_sample, paste0('lsi_ATAC_', name_sample[ix], '_cor_depth.pdf'))
    ggsave(file_out, cor_depth)
    keep_dim_lsi <- which(abs(cor_depth$data$nCount_ATAC) < cutoff_cor_lsi)
    curr_scatac <- RunUMAP(curr_scatac, reduction = 'lsi', dims = keep_dim_lsi, reduction.name = 'umap_ATAC', verbose = FALSE)
    DimPlot(curr_scatac, pt.size = 1.5, reduction = 'umap_ATAC', group.by = 'orig.ident') + labs(x = 'UMAP=1', y = 'UMAP-2') +
                theme_void() + NoLegend() + ggtitle('')
    file_out <- file.path(dir_out_sample, paste0('umap_ATAC_', name_sample[ix], '.png'))
    ggsave(file_out)
    for (jx in plot_meta) {
        FeaturePlot(curr_scatac, features = jx, pt.size = 1.5, reduction = 'umap_ATAC', order = TRUE,
                        cols = c('lightgrey', 'red')) + labs(x = 'UMAP=1', y = 'UMAP-2')
        file_out <- file.path(dir_out_sample, paste0('umap_ATAC_', name_sample[ix], '_meta_', jx, '.png'))
        ggsave(file_out)
    }

    # Recompute and normalize gene activity scores.
    print('- Gene activities')
    gene_activities <- GeneActivity(curr_scatac, extend.upstream = gene_upstream,
                                            extend.downstream = gene_downstream, verbose = FALSE)
    gene_activities <- gene_activities[which(!duplicated(rownames(gene_activities))), ]
    gene_activities <- gene_activities[intersect(rownames(gene_activities), keep_gene), ]
    curr_scatac[['ATAC_Activity']] <- CreateAssayObject(gene_activities, min.cells = 1)
    curr_scatac <- NormalizeData(curr_scatac, assay = 'ATAC_Activity',
                                        scale.factor = median(curr_scatac$nCount_ATAC_Activity),
                                        normalization.method = 'LogNormalize', verbose = FALSE)

    # Save object
    print('- Save object')
    file_out <- file.path(dir_out_sample, paste0('obj_', file_save, '_', name_sample[ix], '_ATAC.rds'))
    saveRDS(curr_scatac, file_out)

}

print('DONE')
