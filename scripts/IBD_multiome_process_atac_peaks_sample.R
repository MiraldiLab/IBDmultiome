# IBD_multiome_process_atac_peaks_sample.R
#
# Purpose: Cluster QC-filtered ATAC cells within each sample, merge clusters
# with insufficient fragment depth, and call MACS2 peaks per resulting cluster.
# Export filtered combined peaks, marker coverage tracks, and diagnostic plots.
#
# Setup: Replace all /path/to/... placeholders, including macs2.path below.
# Required sample metadata column: Sample.
# Filtered objects: <dir_filtered>/<Sample>/<file_filter_prefix>_<Sample>_ATAC.rds.
# Objects must contain an ATAC assay, lsi and umap_ATAC reductions,
# atac_fragments metadata, and gene annotation for marker coverage plots.
# Processed cut-site files: <dir_fragments>/fragments_<Sample>.tsv.gz, with indexes.
# The blacklist is a three-column BED file.
# Helpers must provide bedtools_sort() and single_cell_seurat_pseudobulk_v5().
# BEDTools and igvtools must be available on PATH; MACS2 uses the configured path.
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
    library(ggpubr)
})
source('/path/to/helper_scripts/bedtools_utils.R')
source('/path/to/helper_scripts/scRNA_utils.R')

# Inputs and peak-calling parameters ------------------------------------------

# Output directory
dir_out <- '/path/to/output/atac_peaks_sample'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Sample metadata
file_meta <- '/path/to/sample_metadata/meta_multiomeIBD_sample.txt'

# Filtered objects directory
dir_filtered <- '/path/to/qc_filter'

# Fragment files directory
dir_fragments <- '/path/to/processed_atac_fragments'

# Blacklist bed file
file_blacklist <- '/path/to/blacklist/blacklist_hg38_ENCODE_maxATAC_extended.bed'

# Analysis parameters
file_save <- 'multiomeIBD' # file save base
file_filter_prefix <- 'obj_multiomeIBD'
clust_res <- 0.3 # sample cluster resolution
chr_keep <- c(paste0('chr', 1:22), 'chrX') # chromosomes to keep
macs_min_frags <- 3E6 # min number of fragments per cluster for MACS peak calls
macs_sig_type <- 'FDR' # MACS2 type of significance cutoff (options: 'FDR', 'P')
macs_sig_cutoff <- 0.001 # Significance cutoff for MACS2 peak calling
macs_ext <- 200 # MACS2 extension
macs_shift <- -macs_ext/2 # MACS2 shift
macs_genome <- 2.9e9 # MACS2 effective genome size (hg38=2.9e9)
macs_lambda <- TRUE # use local lambda for background
macs_win_llocal <- 10000 # large local window size (MACS2 default: 10000)
genome_name <- 'hg38'
k_lsi <- 30 # dimension reduction individual samples
cutoff_cor_lsi <- 0.75 # cutoff correlation with count depth to exclude LSI components
markers <- c('CD3E', 'CD19', 'CD14', 'MUC2', 'ACTA2', 'PDGFRA', 'APOA1', 'ACE2', 'SELL') # Plot marker tracks

# Configure output directory and load inputs ---------------------------------

# Output directory
dir_out <- file.path(dir_out, paste0('sizeFrag', macs_ext, '_', macs_sig_type,
                                    ifelse(macs_sig_type == 'FDR', gsub('[:.:]', 'p', 100*macs_sig_cutoff),
                                    -log10(macs_sig_cutoff)), '_',
                                    ifelse(macs_lambda, 'lambda', 'nolambda'),
                                    ifelse(macs_lambda, paste0('_llocal', macs_win_llocal), '')))
dir.create(dir_out, showWarnings = FALSE)

# Load sample metadata
meta_sample <- read.delim(file_meta, header = TRUE, sep = '\t')
name_sample <- meta_sample$Sample

# Load blacklist
genome_blacklist <- read.delim(file_blacklist, header = FALSE, sep = '\t')
colnames(genome_blacklist) <- c('chr', 'start', 'end')
genome_blacklist <- subset(genome_blacklist, chr %in% chr_keep)
genome_blacklist <- makeGRangesFromDataFrame(genome_blacklist)

# Process samples and call peaks per cluster ---------------------------------
print('Call MACS2 peaks per cluster')
for (ix in 1:length(name_sample)) {

    # Load scATAC
    print(paste0('Processing: ', name_sample[ix]))
    dir_out_sample <- file.path(dir_out, name_sample[ix])
    dir.create(dir_out_sample, showWarnings = FALSE)
    file_scatac_obj <- file.path(dir_filtered, name_sample[ix], paste0(file_filter_prefix, '_',
                                                                        name_sample[ix], '_ATAC.rds'))
    curr_scatac <- readRDS(file_scatac_obj)
    DefaultAssay(curr_scatac) <- 'ATAC'

    # Cluster cells
    print('- Cluster')
    cor_depth <- DepthCor(curr_scatac, assay = 'ATAC', reduction = 'lsi', n = k_lsi)
    keep_dim_lsi <- which(abs(cor_depth$data$nCount_ATAC) < cutoff_cor_lsi)
    curr_scatac <- FindNeighbors(curr_scatac, reduction = 'lsi', dims = keep_dim_lsi, verbose = FALSE)
    curr_scatac <- FindClusters(curr_scatac, resolution = clust_res, verbose = FALSE)
    DimPlot(curr_scatac, label = TRUE, label.size = 6, reduction = 'umap_ATAC') + labs(x = 'UMAP-1', y = 'UMAP-2') + NoLegend()
    file_out <- file.path(dir_out_sample, paste0('umap_', name_sample[ix], '_cluster_Res', clust_res, '.png'))
    ggsave(file_out)

    # Merge clusters below the fragment-depth threshold.
    curr_scatac[['currClust']] <- Idents(curr_scatac)
    size_clust_frag <- aggregate(atac_fragments ~ currClust,
                                    curr_scatac@meta.data[, c('atac_fragments', 'currClust')], sum)
    size_clust_frag <- setNames(size_clust_frag$atac_fragments, as.character(size_clust_frag$currClust))
    while (min(size_clust_frag) < macs_min_frags) {
        if (length(levels(curr_scatac)) == 2) {
            Idents(curr_scatac) <- 0
            size_clust_frag <- macs_min_frags + 1
        } else {
            # Merge the smallest cluster using pseudobulk Spearman correlation.
            curr_counts <- single_cell_seurat_pseudobulk_v5(curr_scatac)
            curr_meta <- data.frame(row.names = as.character(levels(curr_scatac)), Cluster = levels(curr_scatac))
            curr_meta <- curr_meta[as.character(colnames(curr_counts)), , drop = FALSE]
            dds <- DESeqDataSetFromMatrix(countData = curr_counts, colData = curr_meta, design = ~ 1)
            counts_norm <- varianceStabilizingTransformation(dds)
            counts_norm <- t(scale(t(assay(counts_norm))))
            mat_cor <- cor(counts_norm, method = 'spearman')
            min_ident <- names(which(size_clust_frag == min(size_clust_frag)))[1]
            sort_cor <- sort(mat_cor[as.character(min_ident), ], decreasing = TRUE)
            Idents(curr_scatac, cells = WhichCells(curr_scatac, idents = as.numeric(min_ident))) <-
                                        as.numeric(names(sort_cor)[2])
            curr_scatac[['currClust']] <- Idents(curr_scatac)
            size_clust_frag <- aggregate(atac_fragments ~ currClust,
                                            curr_scatac@meta.data[, c('atac_fragments', 'currClust')], sum)
            size_clust_frag <- setNames(size_clust_frag$atac_fragments, as.character(size_clust_frag$currClust))
        }
    }
    DimPlot(curr_scatac, label = TRUE, label.size = 6, reduction = 'umap_ATAC') + labs(x = 'UMAP-1', y = 'UMAP-2') + NoLegend()
    file_out <- file.path(dir_out_sample, paste0('umap_', name_sample[ix], '_cluster_Res', clust_res, '_merged.png'))
    ggsave(file_out)

    # Save cluster idents
    curr_scatac[['cluster_MACS']] <- Idents(curr_scatac)

    # Replace the fragment reference with the processed cut-site file.
    Fragments(curr_scatac@assays$ATAC) <- NULL
    file_scatac_frag <- file.path(dir_fragments, paste0('fragments_', name_sample[ix], '.tsv.gz'))
    new_frag <- CreateFragmentObject(path = file_scatac_frag, cells = colnames(curr_scatac), validate.fragments = TRUE)
    Fragments(curr_scatac@assays$ATAC) <- new_frag

    # Call MACS2 peaks
    print('- Call peaks')
    dir_out_macs <- file.path(dir_out_sample, 'MACS2')
    dir.create(dir_out_macs, showWarnings = FALSE)
    peaks <- CallPeaks(curr_scatac, group.by = 'cluster_MACS',
                        macs2.path = '/path/to/macs2/bin/macs2',
                        extsize = macs_ext, shift = macs_shift, effective.genome.size = macs_genome,
                        additional.args = paste0('-', ifelse(macs_sig_type == 'FDR', 'q', 'p'), ' ',
                            macs_sig_cutoff, ' --nomodel --keep-dup all -B --SPMR --verbose 0',
                            ifelse(macs_lambda, '', ' --nolambda'),
                            ifelse(macs_lambda, paste0(' --llocal ', macs_win_llocal), '')),
                        verbose = FALSE, outdir = dir_out_macs, cleanup = FALSE, combine.peaks = TRUE)

    # Convert bdg to tdf and clean up outputs
    print('- Process MACS2 outs')
    for (jx in levels(curr_scatac)) {
        system(paste0('igvtools toTDF -z 3 ', dir_out_macs, '/', jx, '_treat_pileup.bdg '
                                            , dir_out_macs, '/', jx, '_treat_pileup.tdf ', genome_name))
        system(paste0('rm ', dir_out_macs, '/', jx, '_treat_pileup.bdg'))
        system(paste0('rm ', dir_out_macs, '/', jx, '_control_lambda.bdg'))
        system(paste0('rm ', dir_out_macs, '/', jx, '_peaks.xls'))
        system(paste0('rm ', dir_out_macs, '/', jx, '_summits.bed'))
    }
    system('rm igv.log')

    # Filter chromosomes
    print('- Filter chromosomes')
    peaks <- peaks[as.character(peaks@seqnames) %in% chr_keep]

    # Remove blacklist regions
    print('- Remove blacklist regions')
    peaks <- subsetByOverlaps(x = peaks, ranges = genome_blacklist, invert = TRUE)

    # Marker tracks
    print('- Marker tracks')
    Idents(curr_scatac) <- factor(Idents(curr_scatac), levels = sort(levels(curr_scatac)))
    curr_scatac@assays$ATAC@annotation$tx_id <- curr_scatac@assays$ATAC@annotation$transcript_id
    dir_out_markers <- file.path(dir_out_sample, 'markers')
    dir.create(dir_out_markers, showWarnings = FALSE)
    list_p_markers <- list()
    for (jx in markers) {
        curr_region <- LookupGeneCoords(curr_scatac, jx)
        list_p_markers[[jx]] <- CoveragePlot(curr_scatac, group.by = 'cluster_MACS', region = curr_region, ranges = peaks,
                                                extend.upstream = 5000, extend.downstream = 5000, ranges.title = 'MACS2')
        file_out <- file.path(dir_out_markers, paste0('track_marker_', file_save, '_', name_sample[ix], '_', jx, '.pdf'))
        ggsave(file_out, plot = list_p_markers[[jx]], height = 5, width = 6)
    }
    ggarrange(plotlist = list_p_markers, nrow = 3, ncol = 3)
    file_out <- file.path(dir_out_sample, paste0('markers_', file_save, '_', name_sample[ix], '_peaks.pdf'))
    ggsave(file_out, width = 15, height = 15)

    # Summarize chromosome-filtered narrowPeak calls before blacklist removal.
    df_peaks_n <- NULL
    df_peaks_size <- NULL
    for (jx in levels(curr_scatac)) {
        file_in <- file.path(dir_out_macs, paste0(jx, '_peaks.narrowPeak'))
        curr_peaks <- read.delim(file_in, header = FALSE)
        curr_peaks <- subset(curr_peaks, V1 %in% chr_keep)
        df_peaks_n <- rbind(df_peaks_n, data.frame(Cluster = jx, Npeaks = nrow(curr_peaks)))
        df_peaks_size <- rbind(df_peaks_size, data.frame(Cluster = jx, Lpeaks = abs(curr_peaks$V3-curr_peaks$V2)))
    }
    df_ncell <- as.data.frame(table(Idents(curr_scatac)))
    colnames(df_ncell) <- c('Cluster', 'Ncell')
    df_nfrag <- aggregate(atac_fragments ~ cluster_MACS,
                                    curr_scatac@meta.data[, c('cluster_MACS', 'atac_fragments')], sum)
    colnames(df_nfrag) <- c('Cluster', 'Nfrag')
    df_peaks_n <- merge(merge(df_peaks_n, df_ncell, by = 'Cluster'), df_nfrag, by = 'Cluster')

    # Plot summary of peak calling
    print('- Plot summary of peak calling')
    p_umap <- DimPlot(curr_scatac, label = TRUE, label.size = 6, reduction = 'umap_ATAC') + theme_void() + NoLegend()
    p_ncell <- ggplot(data = df_peaks_n, aes(x = Cluster, y = Ncell)) + geom_bar(stat = 'identity', width = 0.7) +
                        theme_minimal() +
                        labs(x = 'Cluster', y = 'Cells') +
                        theme(axis.title = element_text(size = 12, face = 'plain')) +
                        theme(axis.text = element_text(size = 11, face = 'plain'))
    p_nfrag <- ggplot(data = df_peaks_n, aes(x = Cluster, y = Nfrag)) +
                        geom_hline(yintercept = macs_min_frags, lwd = 0.5, linetype = 'dashed', colour = 'red') +
                        geom_bar(stat = 'identity', width = 0.7) +
                        theme_minimal() +
                        labs(x = 'Cluster', y = 'Fragments') +
                        theme(axis.title = element_text(size = 12, face = 'plain')) +
                        theme(axis.text = element_text(size = 11, face = 'plain'))
    p_npeaks <- ggplot(data = df_peaks_n, aes(x = Cluster, y = Npeaks)) + geom_bar(stat = 'identity', width = 0.7) +
                        theme_minimal() +
                        labs(x = 'Cluster', y = 'Peaks') +
                        theme(axis.title = element_text(size = 12, face = 'plain')) +
                        theme(axis.text = element_text(size = 11, face = 'plain'))
    p_size <- ggplot(df_peaks_size, aes(x = Cluster, y = Lpeaks)) +
                        geom_violin(width = 0.8) + geom_boxplot(width = 0.3) +
                        theme_minimal() +
                        labs(x = '', y = 'Peak size [bp]') +
                        theme(axis.title = element_text(size = 12, face = 'plain')) +
                        theme(axis.text = element_text(size = 11, face = 'plain'))
                        # scale_y_continuous(trans='log10') + annotation_logticks(sides='l')
    p_summary <- ggarrange(p_umap, ggarrange(plotlist = list(p_ncell, p_nfrag, p_npeaks, p_size), nrow = 4, ncol = 1), nrow = 1)
    annotate_figure(p_summary, top = text_grob(paste0('Sample: ', name_sample[ix]), color = 'black', face = 'plain', size = 16))
    file_out <- file.path(dir_out_sample, paste0('summary_', file_save, '_', name_sample[ix], '_peaks.pdf'))
    ggsave(file_out, width = 11, height = 7)

    # Export the combined peaks after chromosome and blacklist filtering.
    # Preserve the original coordinate handling in the BED-style table.
    print('- Save peaks')
    peaks_bed <- as.data.frame(peaks)[, 1:3]
    peaks_bed <- bedtools_sort(peaks_bed)
    file_out <- file.path(dir_out_sample, paste0('peaks_', file_save, '_', name_sample[ix], '.rds'))
    saveRDS(peaks, file_out)
    file_out <- file.path(dir_out_sample, paste0('peaks_', file_save, '_', name_sample[ix], '.bed'))
    write.table(peaks_bed, file_out, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

}

print('DONE')
