# IBD_multiome_analysis_transfer.R
# module load R/4.2.3 bedtools/2.30.0 hdf5/1.8.13 igvtools/2.4.0
# Transfer reference cell-type labels to individual multiome samples.
# Export prediction scores and filtered labels, with UMAP visualizations.
#
# Replace the /path/to placeholders before running.
# Sample metadata must include Sample. Query objects must contain the
# configured query assay and reduction, plus orig.ident metadata.
# Reference objects must contain the selected assay, variable features,
# and the annotation column specified in transfer_meta.
# Keep file_transfer, name_transfer, transfer_meta, and ref_assay aligned:
# entries at the same position describe the same reference atlas.
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
dir_out <- '/path/to/output/label_transfer'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Sample metadata
file_meta <- '/path/to/sample_metadata/meta_multiomeIBD_sample.txt'

# Query data directory
dir_query <- '/path/to/qc_filter'
query_type <- 'Multi' # query data (Multi or ATAC)
query_assay <- 'RNA' # query assay
query_reduction <- 'umap_RNA' # query dimension reduction

# Reference atlas objects
file_transfer <- c(
'/path/to/reference_atlases/scrna_Elmentaite2021_Hs_atlas_gut.rds',
    '/path/to/reference_atlases/scrna_Smillie2019_Hs_atlas_UC_subMax2000.rds'
)
# Reference names used in metadata columns and output filenames
name_transfer <- c(
'Elmentaite2021',
    'Smillie2019'
)
# Reference annotation columns
transfer_meta <- c(
'annotation_V2',
    'Cluster'
)
# Reference assays
ref_assay <- c(
'RNA',
    'RNA'
)

# Analysis parameters
file_save <- 'multiomeIBD' # file save base
file_filter_prefix <- 'obj_multiomeIBD'
thresh_score <- 0.2 # Predictions below this score are labeled Unclassified.

# ---- Load sample metadata -----------------------------------------------------

# Load sample metadata
meta_sample <- read.delim(file_meta, header = TRUE, sep = '\t')
name_sample <- meta_sample$Sample

# ---- Transfer labels for each sample -------------------------------------------
for (ix in 1:length(name_sample)) {

    # Label transfer output directory
    print(paste0('Sample ', ix, '/', length(name_sample), ': ', name_sample[ix]))
    dir_out_transfer <- file.path(dir_out, name_sample[ix], query_type)
    dir.create(dir_out_transfer, showWarnings = FALSE, recursive = TRUE)

    # Load data
    print('- Load data')
    file_in_data <- file.path(dir_query, name_sample[ix], paste0(file_filter_prefix, '_', name_sample[ix],
                                '_', query_type, '.rds'))
    scrna <- readRDS(file_in_data)
    DefaultAssay(scrna) <- query_assay
    Idents(scrna) <- 'orig.ident'

    # Process each reference atlas in the configured order.
    print('- Label transfer')
    n_transfer <- length(file_transfer)
    all_transfer_labels <- NULL
    for (jx in 1:n_transfer) {

        # Reference output directory
        curr_name_transfer <- name_transfer[jx]
        print(paste0('  - Dataset: ', curr_name_transfer))
        dir_out_transfer_ref <- file.path(dir_out_transfer, curr_name_transfer)
        dir.create(dir_out_transfer_ref, showWarnings = FALSE)

        # Load reference
        data_transfer <- readRDS(file_transfer[jx])
        Idents(data_transfer) <- transfer_meta[jx]
        curr_assay <- ref_assay[jx]

        # Find transfer anchors
        print('  - find transfer anchors')
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

    # Export all reference predictions and filtered labels for this sample.
    meta_transfer <- scrna@meta.data[, all_transfer_labels]
    file_out <- file.path(dir_out, paste0('meta_transfer_', file_save, '_', name_sample[ix], '_',
                                                query_type, '_', query_assay, '.txt'))
    write.table(meta_transfer, file_out, quote = FALSE, col.names = NA, row.names = TRUE, sep = '\t')

}

print('DONE')
