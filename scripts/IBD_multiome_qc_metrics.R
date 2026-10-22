# IBD_multiome_qc_metrics.R
#
# Purpose: Calculate and visualize RNA and ATAC quality-control metrics for
# Cell Ranger ARC-called cells, then export metadata for cells passing RNA,
# ATAC, or joint QC criteria.
#
# Setup: Replace all /path/to/... placeholders with local paths.
# The plot_utils.R helper must provide gg_vln_bar().
# Required sample metadata columns:
#   Sample, dir_cellranger_arc, dir_cellbender, DiseaseTissueStatus.
# DiseaseTissueStatus is used to order samples (e.g., CD_R_NI).
# Cell Ranger ARC directories must contain outs/filtered_feature_bc_matrix.h5,
# outs/atac_fragments.tsv.gz (with its index), and outs/per_barcode_metrics.csv.
# Each CellBender directory must contain <directory_basename>.h5.
# Doublet tables must contain Sample and Cell columns. The filtered gene list
# contains one gene name per line; the blacklist is a three-column BED file.
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
    library(ggpubr)
})
source('/path/to/helper_scripts/plot_utils.R')

# Inputs and QC parameters ---------------------------------------------------

# Output directory
dir_out <- '/path/to/output/qc_metrics'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Sample metadata
file_meta <- '/path/to/sample_metadata/meta_multiomeIBD_sample.txt'

# Cell doublets (RNA)
file_doublets_scrna <- '/path/to/rna_doublets/cells_multiomeIBD_qc_scrna_DoubletFinder.txt'

# Cell doublets (ATAC)
file_doublets_scatac <- '/path/to/atac_doublets/cells_multiomeIBD_qc_scatac_ArchR_doublets_ratio0p5.txt'

# Gene annotation
file_annot <- '/path/to/cellranger_arc_reference/genes/genes.gtf'

# Filtered genes
file_gene <- '/path/to/gene_filter/gene_set_refdata-cellranger-arc-GRCh38-2020-A-2.0.0_filtered.txt'

# Blacklist bed file
file_blacklist <- '/path/to/blacklist/blacklist_hg38_ENCODE_maxATAC_extended.bed'

# Analysis and plotting parameters
file_save <- 'multiomeIBD'
meta_group <- c('Disease', 'Tissue', 'Status')
order_disease <- c('CD', 'UC')
order_tissue <- c('R', 'TI')
order_status <- c('NI', 'I')
genome <- 'hg38'
keep_chr <- paste0('chr', c(1:22, 'X', 'M')) # chromosomes to keep
pattern_pct_mito <- '^MT-' # gene pattern for percent mitochondrial expression
pattern_pct_ribo <- '^RP[SL]' # gene pattern for percent ribosome expression
min_atac_frag <- 2000 # minimum ATAC fragments
max_atac_frag <- 1E5 # maximum ATAC fragments
min_atac_tss_enrich <- 3 # minimum ATAC TSS enrichment
max_atac_pct_mito <- 10 # maximum ATAC percent mitochondrial fragments
max_atac_pct_blacklist <- 3 # maximum ATAC percent blacklist fragments
min_atac_nucsig <- 0.2 # minimum nucleosome signal (ratio: mono-nucleosome/nucleosome-free)
max_atac_nucsig <- 2 # maximum nucleosome signal (ratio: mono-nucleosome/nucleosome-free)
min_atac_pct_reads_in_peaks <- 20 # minimum ATAC percent reads in peaks
min_nCount_RNA <- 1000 # minimum RNA UMI counts
max_nCount_RNA <- 20000 # maximum RNA UMI counts
min_nFeature_RNA <- 200 # minimum RNA genes
max_nFeature_RNA <- 8000 # maximum RNA genes
max_gex_pct_mito <- 25 # maximum RNA percent mitochondrial gene
meta_cell_remove <- c('dir_cellranger_arc', 'dir_cellranger_count', 'dir_cellbender') # remove from cell-level metadata

# QC metrics, display labels, and threshold lines -----------------------------
plot_meta <- c('nCount_RNA', 'nFeature_RNA', 'gex_pct_mito', 'gex_pct_ribo',
                    'atac_fragments', 'nCount_ATAC',
                    'nFeature_ATAC', 'atac_pct_reads_in_peaks',
                    'atac_pct_mito', 'atac_pct_blacklist', 'TSS.enrichment', 'nucleosome_signal')
label_meta <- c('RNA: Total UMIs', 'RNA: Genes', 'RNA: % Mitochondrial', 'RNA: % Ribosome',
                    'ATAC: Total Fragments', 'ATAC: Cut Sites in Peaks',
                    'ATAC: Peaks', 'ATAC: % Reads in Peaks',
                    'ATAC: % Mitochondrial', 'ATAC: % Blacklist', 'TSS Enrichment Score',
                    'Mononucleosome / Nucleosome-free')
plot_log <- c(TRUE, TRUE, FALSE, FALSE,
                    TRUE, TRUE,
                    TRUE, FALSE,
                    FALSE, FALSE, FALSE, FALSE) # plot log10 values
plot_meta_min <- c(min_nCount_RNA, min_nFeature_RNA, NA, NA,
                    min_atac_frag, NA,
                    NA, min_atac_pct_reads_in_peaks,
                    NA, NA, min_atac_tss_enrich, min_atac_nucsig)
plot_meta_max <- c(max_nCount_RNA, max_nFeature_RNA, max_gex_pct_mito, NA,
                    max_atac_frag, NA,
                    NA, NA,
                    max_atac_pct_mito, max_atac_pct_blacklist, NA, max_atac_nucsig)

# Load reference data and sample metadata ------------------------------------

# Load sample metadata
meta_sample <- read.delim(file_meta, header = TRUE, sep = '\t')
name_sample <- meta_sample$Sample
meta_cell <- meta_sample[, setdiff(colnames(meta_sample), meta_cell_remove)]

# Load filtered genes
keep_gene <- readLines(file_gene)

# Load gene annotation
annotations <- rtracklayer::import(file_annot)
annotations$gene_biotype <- annotations$gene_type
annotations <- subset(annotations, gene_name %in% keep_gene)
annotations <- keepSeqlevels(annotations, keep_chr, pruning.mode = 'coarse')
seqlevelsStyle(annotations) <- 'UCSC'

# Load blacklist
genome_blacklist <- read.delim(file_blacklist, header = FALSE, sep = '\t')
colnames(genome_blacklist) <- c('chr', 'start', 'end')
genome_blacklist <- subset(genome_blacklist, chr %in% keep_chr)
genome_blacklist <- makeGRangesFromDataFrame(genome_blacklist)

# Load doublets RNA
cells_doublets_scrna <- read.delim(file_doublets_scrna, header = TRUE, sep = '\t')

# Load doublets ATAC
cells_doublets_scatac <- read.delim(file_doublets_scatac, header = TRUE, sep = '\t')

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

# Calculate and plot metrics for each sample ---------------------------------
print('*** Process metadata ***')
meta_all <- NULL
for (ix in 1:length(name_sample)) {

    # Output directory
    print(paste0('Sample ', ix, '/', length(name_sample), ': ', name_sample[ix]))
    dir_out_sample <- file.path(dir_out, name_sample[ix])
    dir.create(dir_out_sample, showWarnings = FALSE)

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

    # Seurat object
    sc <- CreateSeuratObject(sc_data_scrna, min.cells = 1)
    sc[['gex_pct_mito']] <- PercentageFeatureSet(sc, pattern = pattern_pct_mito)
    sc[['gex_pct_ribo']] <- PercentageFeatureSet(sc, pattern = pattern_pct_ribo)
    sc[['ATAC']] <- chrom_assay

    # TSS Enrichment
    print('- TSS Enrichment')
    sc <- TSSEnrichment(sc, fast = FALSE, assay = 'ATAC', verbose = FALSE, process_n = 10000)
    Idents(sc) <- ' '
    TSSPlot(sc, assay = 'ATAC') + ggtitle('TSS Enrichment') + theme(plot.title = element_text(hjust = 0.5))
    file_out <- file.path(dir_out_sample, paste0('qc_', file_save, '_raw_', name_sample[ix], '_TSSenrich.pdf'))
    ggsave(file_out, width = 4, height = 4)

    # Nucleosome Signal
    print('- Nucleosome Signal')
    sc <- NucleosomeSignal(sc, assay = 'ATAC', verbose = FALSE)
    FragmentHistogram(sc, assay = 'ATAC', region = 'chr1-1-10000000') + ggtitle('Nucleosome Signal') +
                            theme(plot.title = element_text(hjust = 0.5))
    file_out <- file.path(dir_out_sample, paste0('qc_', file_save, '_raw_', name_sample[ix], '_FragHist.pdf'))
    ggsave(file_out, width = 4, height = 4)

    # Barcode metadata
    file_meta_sc <- paste0(meta_sample$dir_cellranger_arc[ix], '/outs/per_barcode_metrics.csv')
    meta_sc <- read.delim(file_meta_sc, header = TRUE, sep = ',')
    rownames(meta_sc) <- meta_sc$barcode
    meta_sc <- meta_sc[colnames(sc), ]
    meta_sc <- meta_sc[, c('barcode', 'gex_mapped_reads', 'atac_fragments',
                                    'atac_mitochondrial_reads', 'atac_peak_region_fragments')]
    meta_sc$atac_pct_mito <- signif((meta_sc$atac_mitochondrial_reads/
                                    (meta_sc$atac_mitochondrial_reads+meta_sc$atac_fragments))*100, 4)

    # Fraction of ATAC counts in blacklist regions (reported as a percentage).
    meta_sc$atac_pct_blacklist <- signif(FractionCountsInRegion(sc, assay = 'ATAC', regions = genome_blacklist)*100, 4)

    # Percentage of ATAC fragments in peak regions from Cell Ranger ARC metrics.
    meta_sc$atac_pct_reads_in_peaks <- signif(meta_sc$atac_peak_region_fragments/meta_sc$atac_fragments*100, 4)

    # Annotate RNA and ATAC doublet calls independently.
    curr_doublets_scatac <- subset(cells_doublets_scatac, Sample == name_sample[ix])
    meta_sc$isDoubletATAC <- 0
    if (length(intersect(curr_doublets_scatac$Cell, rownames(meta_sc))) > 0) {
        meta_sc[intersect(curr_doublets_scatac$Cell, rownames(meta_sc)), 'isDoubletATAC'] <- 1
    }
    curr_doublets_scrna <- subset(cells_doublets_scrna, Sample == name_sample[ix])
    meta_sc$isDoubletRNA <- 0
    if (length(intersect(curr_doublets_scrna$Cell, rownames(meta_sc))) > 0) {
        meta_sc[intersect(curr_doublets_scrna$Cell, rownames(meta_sc)), 'isDoubletRNA'] <- 1
    }

    # Metadata
    curr_meta <- cbind(sc@meta.data[, c('nCount_RNA', 'nFeature_RNA', 'gex_pct_mito', 'gex_pct_ribo',
                                        'nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'nucleosome_signal')], meta_sc)
    curr_meta$Sample <- name_sample[ix]
    curr_meta <- merge(curr_meta, meta_cell, by = 'Sample')
    rownames(curr_meta) <- curr_meta$barcode
    meta_all <- rbind(meta_all, curr_meta)

    # Save metadata
    file_out <- file.path(dir_out_sample, paste0('meta_', file_save, '_', name_sample[ix], '.txt'))
    write.table(curr_meta, file_out, quote = FALSE, sep = '\t', col.names = NA, row.names = TRUE)

    # Scatter plot metadata
    p_meta <- list()
    for (jx in 1:length(plot_meta)) {
        ggplot(curr_meta) +
            geom_point(aes_string(x = 'atac_fragments', y = 'nCount_RNA', color = plot_meta[jx]), size = 0.5) +
            theme_minimal() +
            labs(x = 'ATAC: Total Fragments', y = 'RNA: Total UMIs', color = label_meta[jx]) +
            theme(axis.title = element_text(size = 14, face = 'plain')) +
            theme(axis.text = element_text(size = 11, face = 'plain')) +
            theme(legend.text = element_text(size = 12)) +
            scale_y_continuous(trans = 'log10') + annotation_logticks(sides = 'l') +
            scale_x_continuous(trans = 'log10') + annotation_logticks(sides = 'b') +
            scale_color_gradient(low = 'grey70', high = 'red')
        file_out <- file.path(dir_out_sample,
                                paste0('qc_', file_save, '_raw_scatter_', name_sample[ix], '_', plot_meta[jx], '.pdf'))
        ggsave(file_out, height = 4, width = 5.5)

        p_meta[[jx]] <- ggplot(curr_meta, aes_string(x = 'Sample', y = plot_meta[jx])) +
                geom_violin(width = 0.8) + geom_boxplot(width = 0.3, fill = 'white') +
                labs(x = '', y = label_meta[jx]) +
                theme(axis.text = element_text(size = 12, face = 'plain')) +
                theme(axis.title = element_text(size = 14, face = 'plain')) +
                theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
                theme(legend.position = 'none')
        if (plot_log[jx]) {
            p_meta[[jx]] <- p_meta[[jx]] + scale_y_continuous(trans = 'log10') + annotation_logticks(sides = 'l')
        }
        file_out <- file.path(dir_out_sample,
                                paste0('qc_', file_save, '_raw_bar_', name_sample[ix], '_', plot_meta[jx], '.pdf'))
        ggsave(file_out, plot = p_meta[[jx]], height = 5, width = 2)
    }

    # Summary of all QC metrics
    dir_summary <- file.path(dir_out, 'summary')
    dir.create(dir_summary, showWarnings = FALSE)
    p_tss <- TSSPlot(sc, assay = 'ATAC') + ggtitle('TSS Enrichment') + theme(plot.title = element_text(hjust = 0.5)) +
                theme(plot.title = element_text(hjust = 0.5, size = 16)) +
                theme(axis.text = element_text(size = 12, face = 'plain')) +
                theme(axis.title = element_text(size = 14, face = 'plain'))
    p_nucsig <- FragmentHistogram(sc, assay = 'ATAC', region = 'chr1-1-10000000') + ggtitle('Nucleosome Signal') +
                                    theme(plot.title = element_text(hjust = 0.5, size = 16)) +
                                    theme(axis.text = element_text(size = 12, face = 'plain')) +
                                    theme(axis.title = element_text(size = 14, face = 'plain'))
    p_summary <- ggarrange(ggarrange(p_tss, p_nucsig, nrow = 1), ggarrange(plotlist = p_meta, nrow = 2, ncol = 6),
                                nrow = 2, heights = c(1, 2))
    annotate_figure(p_summary, top = text_grob(paste0('Sample: ', name_sample[ix]), color = 'black', face = 'plain', size = 20))
    file_out <- file.path(dir_summary, paste0('qc_', file_save, '_raw_summary_', name_sample[ix], '.pdf'))
    ggsave(file_out, width = 11, height = 15)

}

# Round selected metrics before export and filtering -------------------------
meta_all$gex_pct_mito <- signif(meta_all$gex_pct_mito, digits = 4)
meta_all$gex_pct_ribo <- signif(meta_all$gex_pct_ribo, digits = 4)
meta_all$TSS.enrichment <- signif(meta_all$TSS.enrichment, digits = 4)
meta_all$nucleosome_signal <- signif(meta_all$nucleosome_signal, digits = 4)
meta_all$atac_pct_mito <- signif(meta_all$atac_pct_mito, digits = 4)
meta_all$atac_pct_blacklist <- signif(meta_all$atac_pct_blacklist, digits = 4)

# Export combined metadata ---------------------------------------------------
file_out <- file.path(dir_out, paste0('meta_', file_save, '_All.txt'))
write.table(meta_all, file_out, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

# Plot QC distributions across samples ---------------------------------------
print('Plot metadata')
meta_all$Sample <- factor(meta_all$Sample, levels = order_sample)
for (ix in 1:length(plot_meta)) {
    p <- gg_vln_bar(meta_all, x = 'Sample', y = plot_meta[ix], fill = 'Sample', type = 'vln',
                        flag_box = TRUE, box_width = 0.2,
                        flag_leg = FALSE, lab_y = label_meta[ix],
                        fill_col = rep('white', length(name_sample)), angle_x = 90, size_text = 10) +
                        theme(axis.text.x = element_text(vjust = 0.5)) +
                        geom_hline(yintercept = plot_meta_min[ix], linetype = 'dashed', color = 'grey50') +
                        geom_hline(yintercept = plot_meta_max[ix], linetype = 'dashed', color = 'grey50')
    file_out <- file.path(dir_out, paste0('qc_', file_save, '_raw_vln_', plot_meta[ix], '.pdf'))
    ggsave(file_out, plot = p, width = 9, height = 5)
    p <- gg_vln_bar(meta_all, x = 'Sample', y = plot_meta[ix], fill = 'Sample', type = 'vln',
                        flag_box = TRUE, box_width = 0.2,
                        flag_leg = FALSE, lab_y = label_meta[ix], flag_log10 = TRUE,
                        fill_col = rep('white', length(name_sample)), angle_x = 90, size_text = 10) +
                        theme(axis.text.x = element_text(vjust = 0.5)) +
                        geom_hline(yintercept = plot_meta_min[ix], linetype = 'dashed', color = 'grey50') +
                        geom_hline(yintercept = plot_meta_max[ix], linetype = 'dashed', color = 'grey50')
    file_out <- file.path(dir_out, paste0('qc_', file_save, '_raw_vln_', plot_meta[ix], '_log10.pdf'))
    ggsave(file_out, plot = p, width = 9, height = 5)

    p <- gg_vln_bar(meta_all, x = 'Sample', y = plot_meta[ix], fill = 'Sample', type = 'box',
                        flag_box = TRUE, box_width = 0.4,
                        flag_leg = FALSE, lab_y = label_meta[ix],
                        fill_col = rep('white', length(name_sample)), angle_x = 90, size_text = 10) +
                        theme(axis.text.x = element_text(vjust = 0.5)) +
                        geom_hline(yintercept = plot_meta_min[ix], linetype = 'dashed', color = 'grey50') +
                        geom_hline(yintercept = plot_meta_max[ix], linetype = 'dashed', color = 'grey50')
    file_out <- file.path(dir_out, paste0('qc_', file_save, '_raw_box_', plot_meta[ix], '.pdf'))
    ggsave(file_out, plot = p, width = 9, height = 5)
    p <- gg_vln_bar(meta_all, x = 'Sample', y = plot_meta[ix], fill = 'Sample', type = 'box',
                        flag_box = TRUE, box_width = 0.4,
                        flag_leg = FALSE, lab_y = label_meta[ix], flag_log10 = TRUE,
                        fill_col = rep('white', length(name_sample)), angle_x = 90, size_text = 10) +
                        theme(axis.text.x = element_text(vjust = 0.5)) +
                        geom_hline(yintercept = plot_meta_min[ix], linetype = 'dashed', color = 'grey50') +
                        geom_hline(yintercept = plot_meta_max[ix], linetype = 'dashed', color = 'grey50')
    file_out <- file.path(dir_out, paste0('qc_', file_save, '_raw_box_', plot_meta[ix], '_log10.pdf'))
    ggsave(file_out, plot = p, width = 9, height = 5)
}

# Select and export cells passing QC ------------------------------------------
# Joint QC excludes cells called doublets by both modalities.
# Each modality-specific set combines its own QC-passing, non-doublet cells
# with all joint-QC cells, then removes duplicate rows.
cells_both <- subset(meta_all, nCount_RNA >= min_nCount_RNA & nCount_RNA <= max_nCount_RNA &
                                nFeature_RNA >= min_nFeature_RNA & nFeature_RNA <= max_nFeature_RNA &
                                gex_pct_mito <= max_gex_pct_mito &
                                atac_fragments >= min_atac_frag & atac_fragments <= max_atac_frag &
                                TSS.enrichment >= min_atac_tss_enrich &
                                atac_pct_mito <= max_atac_pct_mito & atac_pct_blacklist <= max_atac_pct_blacklist &
                                nucleosome_signal >= min_atac_nucsig & nucleosome_signal <= max_atac_nucsig &
                                atac_pct_reads_in_peaks >= min_atac_pct_reads_in_peaks &
                                (isDoubletRNA == 0 | isDoubletATAC == 0))
cells_scatac <- rbind(subset(meta_all, atac_fragments >= min_atac_frag & atac_fragments <= max_atac_frag &
                                TSS.enrichment >= min_atac_tss_enrich &
                                atac_pct_mito <= max_atac_pct_mito & atac_pct_blacklist <= max_atac_pct_blacklist &
                                nucleosome_signal >= min_atac_nucsig & nucleosome_signal <= max_atac_nucsig &
                                atac_pct_reads_in_peaks >= min_atac_pct_reads_in_peaks &
                                isDoubletATAC == 0), cells_both) %>% distinct()
cells_scrna <- rbind(subset(meta_all, nCount_RNA >= min_nCount_RNA & nCount_RNA <= max_nCount_RNA &
                                nFeature_RNA >= min_nFeature_RNA & nFeature_RNA <= max_nFeature_RNA &
                                gex_pct_mito <= max_gex_pct_mito &
                                isDoubletRNA == 0), cells_both) %>% distinct()
file_out <- file.path(dir_out, paste0('meta_', file_save, '_filtered_scrna.txt'))
write.table(cells_scrna, file_out, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
file_out <- file.path(dir_out, paste0('meta_', file_save, '_filtered_scatac.txt'))
write.table(cells_scatac, file_out, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
file_out <- file.path(dir_out, paste0('meta_', file_save, '_filtered_both.txt'))
write.table(cells_both, file_out, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

print('DONE')
