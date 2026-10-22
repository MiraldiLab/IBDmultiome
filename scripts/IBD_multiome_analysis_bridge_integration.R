# IBD_multiome_analysis_bridge_integration.R
# load R/4.4.0-R9 bedtools/2.30.0 hdf5/1.8.13 igvtools/2.4.0 openssl/1.0.2
# Annotate ATAC query cells by bridge integration with an annotated multiome.
# Use the multiome RNA assay as the reference and paired RNA/ATAC as the bridge.
#
# Replace the /path/to placeholders before running. The multiome object must
# contain RNA and ATAC assays, pca and harmony_RNA reductions, a modeled
# umap_WNN reduction, and CellTypeL1 annotations. The ATAC query must contain
# the ATAC assay and lsi, harmony_ATAC, and umap_ATAC reductions.
# tissue_subset labels outputs; it does not filter either input object.
# The active mapping transfers CellTypeL1 only.
# An existing stop() ends execution after the annotated ATAC object is saved;
# the subsequent comparison block is retained but does not run by default.
# https://satijalab.org/seurat/articles/seurat5_integration_bridge.html
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
    library(DESeq2)
})

# ---- Input paths and integration parameters ------------------------------------

# Output directory
dir_out <- '/path/to/output/bridge_integration'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Multiome object (reference & bridge)
file_multi <- '/path/to/annotated_multiome/TI/obj_IBD_multiome_Multi_harmony_TI_annotated.rds'

# Query object
file_atac <- '/path/to/integrate_atac/TI/obj_IBD_multiome_ATAC_harmony_TI.rds'

# Analysis parameters
file_save <- 'IBD_multiome' # file save base
tissue_subset <- 'TI' # name subset of data
k_ref <- 30 # reference dimension reduction
k_query <- 30 # query dimension reduction

# ---- Load and prepare reference, bridge, and query objects ----------------------

# Output directory
dir_out <- file.path(dir_out, tissue_subset)
dir.create(dir_out, showWarnings = FALSE)

# Load multiome object
print('Load multiome object')
obj_multi <- readRDS(file_multi)
DefaultAssay(obj_multi) <- 'ATAC'
keep_peak <- rownames(obj_multi)
DefaultAssay(obj_multi) <- 'RNA'

# Create the RNA reference, pairing Harmony embeddings with PCA loadings.
obj_ref <- DietSeurat(obj_multi, assays = 'RNA', dimreducs = c('pca', 'harmony_RNA', 'umap_WNN'))
obj_ref[['harmony_RNA2']] <- CreateDimReducObject(embeddings = obj_ref[['harmony_RNA']]@cell.embeddings,
                                    key = 'harmonyRNA2_', loadings = obj_ref[['pca']]@feature.loadings, assay = 'RNA')
obj_ref <- RenameCells(obj_ref, add.cell.id = 'ref')
obj_multi <- RenameCells(obj_multi, add.cell.id = 'bridge')

# Load ATAC object
print('Load ATAC object')
obj_atac <- readRDS(file_atac)
DefaultAssay(obj_atac) <- 'ATAC'

# Retain query peaks shared with the multiome and prepare query reductions.
keep_peak <- intersect(rownames(obj_atac), keep_peak)
obj_atac <- DietSeurat(obj_atac, assays = 'ATAC', features = keep_peak,
                            dimreducs = c('lsi', 'harmony_ATAC', 'umap_ATAC'))
obj_atac[['harmony_RNA2']] <- CreateDimReducObject(embeddings = obj_atac[['harmony_ATAC']]@cell.embeddings,
                                    key = 'harmonyATAC2_', loadings = obj_atac[['lsi']]@feature.loadings, assay = 'RNA')
obj_atac[['lsi']] <- obj_atac[['harmony_RNA2']]

# Harmony compatibility with Seurat
# Add PCA feature loadings to harmony reduction
obj_multi[['harmony_RNA2']] <- CreateDimReducObject(embeddings = obj_multi[['harmony_RNA']]@cell.embeddings,
                                    key = 'harmonyRNA2_', loadings = obj_multi[['pca']]@feature.loadings, assay = 'RNA')

# Randomly subsample the bridge to 40,000 cells when needed.
if (ncol(obj_multi) > 4E4) {
  cell_subsample <- intersect(colnames(obj_multi), colnames(obj_multi)[sample(ncol(obj_multi))[1:4E4]])
  obj_multi_sub <- subset(obj_multi, cells = cell_subsample)
} else {
  obj_multi_sub <- obj_multi
}

# ---- Prepare the bridge reference and map the query -----------------------------

# Prepare the extended reference.
print('Prepare bridge reference')
obj_ref_ext <- PrepareBridgeReference(reference = obj_ref, bridge = obj_multi_sub,
                                        reference.reduction = 'harmony_RNA2', reference.dims = 1:k_ref,
                                        normalization.method = 'LogNormalize')

# Find anchors between the extended reference and query object
print('Find anchors')
DefaultAssay(obj_ref_ext) <- 'ATAC'
anchor_bridge <- FindBridgeTransferAnchors(extended.reference = obj_ref_ext, query = obj_atac,
                                            reduction = 'lsiproject', dims = 2:k_query)

# Map query dataset onto reference
print('Map query')
# obj_atac <- MapQuery(anchorset=anchor_bridge, reference=obj_ref_ext,
#                         query=obj_atac, reduction.model='umap_WNN.1',
#                         refdata=list(CellTypeL0='CellTypeL0', CellTypeL1='CellTypeL1',
#                                         CellTypeL2='CellTypeL2'))
obj_atac <- MapQuery(anchorset = anchor_bridge, reference = obj_ref_ext,
                        query = obj_atac, reduction.model = 'umap_WNN.1',
                        refdata = list(CellTypeL1 = 'CellTypeL1'))

# ---- Visualize reference and predicted query labels ----------------------------

# Multiome UMAP
DefaultAssay(obj_multi) <- 'RNA'
Idents(obj_multi) <- 'CellTypeL1'
color_celltype <- setNames(hue_pal()(length(levels(obj_multi))), levels(obj_multi))
color_celltype['ref'] <- 'grey90'
DimPlot(obj_multi, reduction = 'umap_WNN', label = TRUE, repel = TRUE, cols = color_celltype[levels(obj_multi)]) +
        theme_void() + NoLegend()
file_out <- file.path(dir_out, 'umap_multiome_CellTypeL1.png')
ggsave(file_out)
DimPlot(obj_multi, reduction = 'umap_WNN', label = FALSE, cols = rep('grey70', length(levels(obj_multi))), pt.size = 1.1) +
            theme_void() + NoLegend()
file_out <- file.path(dir_out, 'umap_multiome_blank.png')
ggsave(file_out)

# ATAC UMAP
Idents(obj_atac) <- 'predicted.CellTypeL1'
DimPlot(obj_atac, reduction = 'ref.umap', label = TRUE, repel = TRUE, cols = color_celltype[levels(obj_atac)]) +
        theme_void() + NoLegend()
file_out <- file.path(dir_out, 'umap_atac_CellTypeL1_predicted.png')
ggsave(file_out)

# Overlay mapped ATAC cells on the multiome reference coordinates.
umap_ref <- as.data.frame(obj_multi@reductions$umap_WNN@cell.embeddings)
colnames(umap_ref) <- c('UMAP1', 'UMAP2')
umap_ref$CellType <- 'ref'
umap_atac <- as.data.frame(obj_atac@reductions$ref.umap@cell.embeddings)
colnames(umap_atac) <- c('UMAP1', 'UMAP2')
umap_atac$CellType <- obj_atac@meta.data[rownames(umap_atac), 'predicted.CellTypeL1']
umap_all <- rbind(umap_ref, umap_atac)
umap_all$CellType <- factor(umap_all$CellType, levels = c(intersect(levels(obj_multi), unique(umap_all$CellType)), 'ref'))
ggplot(umap_all, aes(x = UMAP1, y = UMAP2, color = CellType)) + geom_point(size = 0.5) + theme_classic() +
            scale_color_manual(values = color_celltype[levels(umap_all$CellType)]) +
            theme(axis.line = element_blank(), axis.text.x = element_blank(),
                    axis.text.y = element_blank(), axis.ticks = element_blank(),
                    axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = 'none')
file_out <- file.path(dir_out, 'umap_atac_CellTypeL1_predicted_overlay.png')
ggsave(file_out)

# ---- Save the annotated ATAC query ---------------------------------------------
print('Save annotated ATAC object')
file_out <- file.path(dir_out, gsub('.rds', '_annotated.rds', basename(file_atac)))
saveRDS(obj_atac, file_out)

# Preserve the original execution boundary after saving the mapped query.
stop()

# ---- Inactive metadata comparison block ----------------------------------------

# This block requires CellTypeL0/L1/L2 and their predicted counterparts.
# The active MapQuery call above transfers only CellTypeL1.
meta_multi <- obj_multi@meta.data[, c('Sample', 'Patient', 'Tissue', 'Disease', 'Status',
                                        'nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal',
                                        'atac_fragments', 'atac_mitochondrial_reads',
                                        'atac_peak_region_fragments', 'atac_pct_mito',
                                        'atac_pct_blacklist', 'atac_pct_reads_in_peaks',
                                        'CellTypeL0', 'CellTypeL1', 'CellTypeL2')]
meta_atac <- obj_atac@meta.data[, c('Sample', 'Patient', 'Tissue', 'Disease', 'Status',
                                        'nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal',
                                        'atac_fragments', 'atac_mitochondrial_reads',
                                        'atac_peak_region_fragments', 'atac_pct_mito',
                                        'atac_pct_blacklist', 'atac_pct_reads_in_peaks',
                                        'predicted.CellTypeL0', 'predicted.CellTypeL1',
                                        'predicted.CellTypeL2')]
colnames(meta_atac) <- colnames(meta_multi)
meta_multi <- subset(meta_multi, CellTypeL1 %in% meta_atac$CellTypeL1)
meta_atac <- subset(meta_atac, CellTypeL1 %in% meta_multi$CellTypeL1)
meta_multi$Type <- 'ATAC+RNA'
meta_atac$Type <- 'ATAC'
meta_all <- rbind(meta_multi, meta_atac)
order_celltype <- setdiff(names(sort(table(meta_all$CellTypeL1), decreasing = TRUE)), 'Low_Quality')
meta_all$CellTypeL1 <- factor(meta_all$CellTypeL1, levels = order_celltype)
meta_all$Type <- factor(meta_all$Type, levels = c('ATAC+RNA', 'ATAC'))

for (ix in c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal',
            'atac_fragments', 'atac_mitochondrial_reads',
            'atac_peak_region_fragments', 'atac_pct_mito',
            'atac_pct_blacklist', 'atac_pct_reads_in_peaks')) {

  curr_meta <- meta_all
  curr_meta$Meta <- curr_meta[, ix]
  ggplot(curr_meta, aes(CellTypeL1, Meta, fill = Type)) + geom_violin(position = position_dodge(), width = 0.8) +
          labs(x = '', y = ix) +
          theme(axis.text.x = element_text(size = 10, face = 'plain', angle = 90, hjust = 1, vjust = 0.5))
  file_out <- file.path(dir_out, paste0('vln_compare_Quality_', ix, '.pdf'))
  ggsave(file_out, height = 5, width = 12)

  ggplot(curr_meta, aes(CellTypeL1, Meta, fill = Type)) + geom_violin(position = position_dodge(), width = 0.8) +
          labs(x = '', y = ix) +
          theme(axis.text.x = element_text(size = 10, face = 'plain', angle = 90, hjust = 1, vjust = 0.5)) +
          scale_y_continuous(trans = 'log10') + annotation_logticks(sides = 'l')
  file_out <- file.path(dir_out, paste0('vln_compare_Quality_', ix, '_log10.pdf'))
  ggsave(file_out, height = 5, width = 12)

}

print('DONE')
