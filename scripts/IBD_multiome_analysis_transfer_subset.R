# IBD_multiome_analysis_transfer_subset.R
# module load R/4.2.3 bedtools/2.30.0 hdf5/1.8.13 igvtools/2.4.0
# Transfer reference cell-type labels to an integrated multiome subset.
# Export prediction scores and filtered labels, with UMAP visualizations.
#
# Replace the /path/to placeholders before running, including commented
# reference paths if those references are enabled.
# The query object must contain the configured assay and reduction, plus
# orig.ident metadata. References require the selected assay and annotation
# column; PCA feature loadings are used when variable features are absent.
# Keep file_transfer, name_transfer, transfer_meta, and ref_assay aligned.
# tissue_subset and query_subset label outputs; they do not filter cells.
#
rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(42)
suppressPackageStartupMessages({
    library(Seurat)
    library(Signac)
    library(ggplot2)
    library(magrittr)
    library(sva)
    library(RColorBrewer)
    library(scales)
    library(dplyr)
})

# ---- Input paths and label-transfer parameters ---------------------------------

# Output directory
dir_out <- '/path/to/output/transfer_subset'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Integrated query object
file_query <- '/path/to/integrate_subset/TI/Stromal/integrated_stromal_multiome.rds'
query_subset <- 'Stromal'
query_assay <- 'RNA' # query assay
query_reduction <- 'umap_WNN' # query dimension reduction

# Reference atlas objects
file_transfer <- c(
    # '/path/to/reference_atlases/scrna_Elmentaite2021_Hs_atlas_gut.rds'
    # '/path/to/reference_atlases/scrna_Smillie2019_Hs_atlas_UC_subMax2000.rds'
    # '/path/to/reference_atlases/clustered_colon_object.rds'
    # '/path/to/reference_atlases/clustered_ileum_object.rds
    # '/path/to/reference_atlases/clustered_immune_object.rds'
    '/path/to/reference_atlases/clustered_stromal_object.rds'
    # '/path/to/reference_atlases/scrna_Elmentaite2021_Hs_atlas_gut_Tcell.rds'
)
# Reference names used in metadata columns and output filenames
name_transfer <- c(
    # 'Elmentaite2021'
    # 'Smillie2019'
    # 'Hickey2023_Epithelial_colon'
    # 'Hickey2023_Epithelial_ileum'
    # 'Hickey2023_Immune'
    'Hickey2023_Stromal'
    # 'Elmentaite2021_Tcell'
)
# Reference annotation columns
transfer_meta <- c(
    # 'annotation_V2'
    # 'Cluster'
    # 'CellType'
    # 'CellType'
    # 'CellType'
    'CellType'
    # 'annotation'
)
# Reference assays
ref_assay <- c(
    # 'RNA'
    # 'RNA'
    # 'decontXcounts'
    # 'decontXcounts'
    # 'decontXcounts'
    'decontXcounts'
    # 'RNA'
)

# Analysis parameters
file_save <- 'multiomeIBD' # file save base
tissue_subset <- 'TI'
thresh_score <- 0.2 # Predictions below this score are labeled Unclassified.

# ---- Load the query subset -----------------------------------------------------

# Output directory
dir_out <- file.path(dir_out, tissue_subset)
dir.create(dir_out, showWarnings = FALSE)

# Load query data
print('Load query data')
dir_out_transfer <- file.path(dir_out, query_subset)
dir.create(dir_out_transfer, showWarnings = FALSE)
scrna <- readRDS(file_query)
DefaultAssay(scrna) <- query_assay
Idents(scrna) <- 'orig.ident'

# ---- Transfer labels from each reference atlas ---------------------------------
print('- Label transfer')
all_transfer_labels <- NULL
for (ix in 1:length(file_transfer)) {

    # Reference output directory
    curr_name_transfer <- name_transfer[ix]
    print(paste0('  - Dataset: ', curr_name_transfer))
    dir_out_transfer_ref <- file.path(dir_out_transfer, curr_name_transfer)
    dir.create(dir_out_transfer_ref, showWarnings = FALSE)

    # Load reference
    data_transfer <- readRDS(file_transfer[ix])
    Idents(data_transfer) <- transfer_meta[ix]
    curr_assay <- ref_assay[ix]

    # Find transfer anchors
    print('  - find transfer anchors')
    if (length(VariableFeatures(data_transfer)) == 0) {
        VariableFeatures(data_transfer) <- rownames(data_transfer@reductions$pca@feature.loadings)
    }
    transfer_anchors <- FindTransferAnchors(reference = data_transfer, query = scrna,
                            reference.assay = curr_assay, query.assay = query_assay, reduction = 'pcaproject',
                            features = intersect(VariableFeatures(data_transfer), rownames(scrna)))

    # Predict labels and append reference-specific metadata columns.
    print('  - predict labels')
    celltype_labels <- TransferData(anchorset = transfer_anchors, refdata = Idents(data_transfer),
                                        weight.reduction = 'pcaproject')
    colnames(celltype_labels) <- paste0('Transfer_', curr_name_transfer, '_', colnames(celltype_labels))
    scrna <- AddMetaData(scrna, metadata = celltype_labels)

    # Identify cells whose maximum prediction score is below the threshold.
    name_pred_max <- paste0('Transfer_', curr_name_transfer, '_prediction.score.max')
    name_pred_id <- paste0('Transfer_', curr_name_transfer, '_predicted.id')
    idx_unclassified <- which(scrna@meta.data[, name_pred_max] < thresh_score)
    unclassified_cells <- rownames(scrna@meta.data)[idx_unclassified]

    # Store filtered labels, with Unclassified as the final factor level.
    Idents(scrna) <- name_pred_id
    scrna <- SetIdent(scrna, cells = unclassified_cells, value = 'Unclassified')
    Idents(scrna) <- factor(Idents(scrna), levels = c(setdiff(levels(Idents(scrna)), 'Unclassified'), 'Unclassified'))
    scrna[[paste0(name_pred_id, '_filtered')]] <- Idents(scrna)
    all_transfer_labels <- c(all_transfer_labels, colnames(celltype_labels), paste0(name_pred_id, '_filtered'))

    # UMAP plots
    DimPlot(scrna, reduction = query_reduction) + labs(x = 'UMAP-1', y = 'UMAP-2')
    file_out <- file.path(dir_out_transfer_ref, paste0('umap_CellType_Transfer_', curr_name_transfer, '.png'))
    ggsave(file_out, height = 6, width = 7.5)
    DimPlot(scrna, label = TRUE, label.size = 5, repel = TRUE, shuffle = TRUE, reduction = query_reduction) +
                NoLegend() + labs(x = 'UMAP-1', y = 'UMAP-2')
    file_out <- file.path(dir_out_transfer_ref, paste0('umap_CellType_Transfer_', curr_name_transfer, '_label.png'))
    ggsave(file_out)

    # Plot prediction scores for each label and the maximum score.
    meta_pred_score <- colnames(scrna@meta.data)[which(grepl(paste0(curr_name_transfer, '_prediction.score'),
                                                                            colnames(scrna@meta.data)))]
    for (kx in meta_pred_score) {
        label_pred_score <- gsub('[:.:]', '_', kx)
        FeaturePlot(scrna, features = kx, reduction = query_reduction, min.cutoff = 'q5', max.cutoff = 'q95', order = TRUE) +
                        labs(x = 'UMAP-1', y = 'UMAP-2')
        file_out <- file.path(dir_out_transfer_ref, paste0('umap_CellType_', label_pred_score, '.png'))
        ggsave(file_out)
    }
    for (kx in unique(Idents(scrna))) {
        label_ident <- gsub('/', '_', kx)
        DimPlot(scrna, cells.highlight = WhichCells(scrna, idents = kx), cols.highlight = 'black',
                    reduction = query_reduction) + NoLegend() + labs(x = 'UMAP-1', y = 'UMAP-2')
        file_out <- file.path(dir_out_transfer_ref,
                                paste0('umap_CellType_Transfer_', curr_name_transfer, '_', label_ident, '.png'))
        ggsave(file_out)
    }

    Idents(scrna) <- 'orig.ident'
}

# ---- Export prediction scores and filtered labels ------------------------------
meta_transfer <- scrna@meta.data[, all_transfer_labels]
file_out <- file.path(dir_out_transfer, paste0('meta_transfer_', file_save, '_', query_subset, '_', query_assay, '.txt'))
write.table(meta_transfer, file_out, quote = FALSE, col.names = NA, row.names = TRUE, sep = '\t')

print('DONE')
