# scRNA_utils.R
# Helpers for Seurat processing, pseudobulk aggregation, batch correction,
# and cellular detection-rate summaries.
#
# Dependencies: Seurat, Matrix, sva (ComBat), and ggplot2 (CDR plots).
# Load Seurat before calling pseudobulk helpers directly; package-loading
# behavior within each function is preserved from the original script.
# Both slot-based and Seurat v5 layer-based pseudobulk helpers are retained.
# No local file paths require configuration.


# scRNA_seurat_workflow ----------------------------------------------------
# Run standard Seurat workflow starting with raw counts matrix
# or 10X directory. Returns processed Seurat object.
# Inputs:
# counts = raw counts matrix or 10X directory
# min_cells = minimum number of cells to keep gene (default: 1)
# min_features = minimum number of features to keep cell (default: 1)
# name_project = object project name (default: 'scRNA')
# meta = cell metadata, rownames=cells, colnames=metadata (default: NULL)
# n_var = number of variable features (default: 2000)
# verbosity = print progress/updates (default: FALSE)
# vars_regress = variables to regress off during scaling (default: NULL)
# k_pca = number of PC dimensions (default: 50)
# k_tsne = number of PCs as input to tSNE (default: 30)
# k_umap = number of PCs as input to UMAP (default: 30)
# res_cluster = clustering resolution parameter, can be vector (default: 0.5)
# flag_tsne = perform tSNE? (default: TRUE)
# flag_cluster = perform clustering? (default: TRUE)
# batch = metadata column identifying batches for ComBat (default: NULL)
scRNA_seurat_workflow <- function(
    counts,
    min_cells = 1,
    min_features = 1,
    name_project = 'scRNA',
    meta = NULL,
    n_var = 2000,
    verbosity = FALSE,
    vars_regress = NULL,
    k_pca = 50,
    k_tsne = 30,
    k_umap = 30,
    res_cluster = 0.5,
    flag_tsne = TRUE,
    flag_cluster = TRUE,
    batch = NULL
) {

    library(Seurat)

    # Load 10X data when a directory path is supplied.
    if (is.character(counts)) {
        if (file.exists(counts)) {
            counts <- Read10X(counts)
        } else {
            stop('10X directory does not exist!')
        }
    }

    # Create the Seurat object.
    obj <- CreateSeuratObject(counts, min.cells = min_cells, min.features = min_features,
                                            project = name_project, assay = 'RNA', meta.data = meta, names.delim = '____')

    # Select variable features.
    obj <- FindVariableFeatures(obj, selection.method = 'vst', nfeatures = n_var, verbose = verbosity)

    # Normalize counts.
    obj <- NormalizeData(obj, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = verbosity)

    # Apply optional batch correction.
    if (!is.null(batch)) {
        mat_data <- as.matrix(GetAssayData(obj, slot = 'data'))
        batch_id <- obj@meta.data[, batch]
        mat_correct <- batch_correct(mat_data, sample_batch = batch_id, type = 'ComBat')
        obj <- SetAssayData(obj, slot = 'data', mat_correct)
    }

    # Scale expression values.
    obj <- ScaleData(obj, features = VariableFeatures(obj), vars.to.regress = vars_regress, verbose = verbosity)

    # Run PCA.
    obj <- RunPCA(obj, npcs = k_pca, features = VariableFeatures(obj), reduction.name = 'pca', verbose = verbosity)

    # Run optional t-SNE.
    if (flag_tsne) {
        obj <- RunTSNE(obj, dims = 1:k_tsne, reduction = 'pca', verbose = verbosity)
    }

    # Run UMAP.
    obj <- RunUMAP(obj, dims = 1:k_umap, reduction = 'pca', verbose = verbosity)

    # Run optional clustering at each requested resolution.
    if (flag_cluster) {
        obj <- FindNeighbors(obj, dims = 1:k_umap, verbose = verbosity)
        for (ix in res_cluster) {
            obj <- FindClusters(obj, resolution = ix, verbose = verbosity)
        }
    }

    return(obj)
}

# scRNA_seurat_pseudobulk --------------------------------------------------
# Create a pseudobulk matrix from a Seurat object
# Inputs:
# obj = Seurat object
# name_ident = metadata column to define identity (default: NULL for current idents)
# type_norm = character, type of normalization (default: NULL returns raw counts)
#   options: 'TPM', 'Log2TPM', 'zscore', 'Log2FC'
# The historical 'TPM' label denotes counts scaled to 1 million per column;
# no gene-length correction is applied. Other modes transform these values.
# Output: a feature-by-identity matrix with all-zero features removed.
scRNA_seurat_pseudobulk <- function(
    obj,
    name_ident = NULL,
    type_norm = NULL
) {
    library(Matrix)

    # Set cell identities.
    if (!is.null(name_ident)) {
        Idents(obj) <- name_ident
    }

    # Sum counts by identity.
    all_idents <- as.character(unique(Idents(obj)))
    counts <- matrix(0, nrow(obj), length(all_idents))
    rownames(counts) <- rownames(obj)
    colnames(counts) <- all_idents
    for (ix in all_idents) {
        curr_counts <- Matrix::rowSums(GetAssayData(subset(obj, idents = ix), slot = 'counts'))
        counts[names(curr_counts), ix] <- as.numeric(curr_counts)
    }
    counts <- counts[which(rowSums(counts) > 0), , drop = FALSE]

    # Apply the requested normalization.
    if (is.null(type_norm)) {
        mat_bulk <- counts
    } else if (type_norm == 'TPM') {
        mat_bulk <- 1E6*sweep(counts, MARGIN = 2, colSums(counts), '/')
    } else if (type_norm == 'Log2TPM') {
        mat_bulk <- log2(1E6*sweep(counts, MARGIN = 2, colSums(counts), '/') + 1)
    } else if (type_norm == 'zscore') {
        mat_bulk <- t(scale(t(log2(1E6*sweep(counts, MARGIN = 2, colSums(counts), '/') + 1))))
    } else if (type_norm == 'Log2FC') {
        mat_tpm <- 1E6*sweep(counts, MARGIN = 2, colSums(counts), '/') + 1
        mat_bulk <- log2(sweep(mat_tpm, MARGIN = 1, rowMeans(mat_tpm), '/'))
    }
    return(mat_bulk)
}

# single_cell_seurat_pseudobulk --------------------------------------------
# Create a pseudobulk matrix from a Seurat object.
# Inputs:
# obj = Seurat object
# name_ident = metadata column to define identity (default: NULL for current idents)
# assay = assay name for pseudobulk (default: NULL for current assay)
# Output: summed raw counts by identity, with all-zero features removed.
single_cell_seurat_pseudobulk <- function(
    obj,
    name_ident = NULL,
    assay = NULL
) {
    library(Matrix)

    # Set assay
    if (!is.null(assay)) {
        DefaultAssay(obj) <- assay
    }

    # Set idents
    if (!is.null(name_ident)) {
        Idents(obj) <- name_ident
    }

    # Slim down object
    obj <- DietSeurat(obj, counts = TRUE, data = TRUE, scale.data = FALSE, assays = DefaultAssay(obj))

    # Sum raw counts
    all_idents <- as.character(unique(Idents(obj)))
    counts <- matrix(0, nrow(obj), length(all_idents))
    rownames(counts) <- rownames(obj)
    colnames(counts) <- all_idents
    for (ix in all_idents) {
        curr_counts <- Matrix::rowSums(GetAssayData(subset(obj, idents = ix), slot = 'counts'))
        counts[names(curr_counts), ix] <- as.numeric(curr_counts)
    }
    counts <- counts[which(rowSums(counts) > 0), , drop = FALSE]
    return(counts)
}

# single_cell_seurat_pseudobulk_v5 -----------------------------------------
# Create pseudobulk matrix from Seurat (v5) object.
# Inputs:
# obj = Seurat object
# name_ident = metadata column to define identity (default: NULL for current idents)
# assay = assay name for pseudobulk (default: NULL for current assay)
# Output: summed raw counts by identity, with all-zero features removed.
single_cell_seurat_pseudobulk_v5 <- function(
    obj,
    name_ident = NULL,
    assay = NULL
) {
    library(Matrix)

    # Set assay
    if (!is.null(assay)) {
        DefaultAssay(obj) <- assay
    }

    # Set idents
    if (!is.null(name_ident)) {
        Idents(obj) <- name_ident
    }

    # Slim down object
    obj <- DietSeurat(obj, assays = DefaultAssay(obj), layers = 'counts')

    # Sum raw counts
    all_idents <- as.character(unique(Idents(obj)))
    counts <- matrix(0, nrow(obj), length(all_idents))
    rownames(counts) <- rownames(obj)
    colnames(counts) <- all_idents
    for (ix in all_idents) {
        curr_counts <- Matrix::rowSums(LayerData(subset(obj, idents = ix), layer = 'counts'))
        counts[names(curr_counts), ix] <- as.numeric(curr_counts)
    }
    counts <- counts[which(rowSums(counts) > 0), , drop = FALSE]
    return(counts)
}

# batch_correct ------------------------------------------------------------
# Batch correct gene expression data.
# Inputs:
# data_matrix = gene expression matrix of the form [genes x samples]
# sample_batch = batch covariate character vector (default: NULL)
# type = string specifying type of batch correction (default: 'ComBat')
# 	options:
# 	'ComBat': ComBat batch correction
# Outputs:
# data_corrected = batch corrected data
batch_correct <- function(
    data_matrix,
    sample_batch = NULL,
    type = 'ComBat') {

    if (type == 'ComBat') {

        # Run ComBat.
        library(sva)
        data_corrected <- ComBat(data_matrix, batch = sample_batch)

    } else {
        stop('Not valid batch correction type')
    }

    return(data_corrected)
}


# scRNA_cdr ----------------------------------------------------------------
# scRNA-seq cellular detection rate (CDR).
# CDR is the fraction of genes with nonzero counts in each cell.
# Inputs:
# counts = counts matrix, rows=genes, cols=cells
# digits = number of significant digits to keep (default: 5)
# Output:
# cdr = cellular detection rate
scRNA_cdr <- function(
    counts,
    digits = 5
) {

    # Calculate cellular detection rate.
    n_gene <- dim(counts)[1]
    cdr <- signif(as.numeric(colSums(sign(abs(counts))))/n_gene, digits = digits)

    return(cdr)
}


# plot_scRNA_cdr -----------------------------------------------------------
# Plot cellular detection rate of scRNA-seq dataset.
# Inputs:
# counts = counts matrix, rows=genes, cols=cells
# flag_hist = plot histogram (default: TRUE)
# flag_dens = overlay density on the histogram (default: FALSE)
# When flag_hist is FALSE, a density plot is returned regardless of flag_dens.
# Outputs:
# plot_cdr = cdr plot
plot_scRNA_cdr <- function(
    counts,
    flag_hist = TRUE,
    flag_dens = FALSE
) {

    library(ggplot2)

    # Calculate cellular detection rate.
    cdr <- scRNA_cdr(counts = counts)

    # Prepare per-cell detection rates for plotting.
    df_cdr <- data.frame(CDR = cdr, row.names = colnames(counts))

    # Plot the detection-rate distribution.
    if (flag_hist) {
        if (flag_dens) {
            plot_cdr <- ggplot(df_cdr, aes(x = CDR)) +
                                geom_histogram(aes(y = ..density..), colour = 'black', fill = 'grey') +
                                geom_density(alpha = .2, fill = "#FF6666") +
                                labs(x = 'Cellular Detection Rate', y = 'Density')
        } else {
            plot_cdr <- ggplot(df_cdr, aes(x = CDR)) + geom_histogram() +
                                labs(x = 'Cellular Detection Rate', y = 'Cells')
        }
    } else {
        plot_cdr <- ggplot(df_cdr, aes(x = CDR)) +
                                geom_density(alpha = .2, fill = "#FF6666") +
                                labs(x = 'Cellular Detection Rate', y = 'Density')
    }

    return(plot_cdr)
}
