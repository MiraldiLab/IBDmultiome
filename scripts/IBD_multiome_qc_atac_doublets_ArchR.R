# IBD_multiome_qc_atac_doublets_ArchR.R
#
# Purpose: Identify putative doublets from ATAC data using ArchR and export
# cell lists before and after doublet filtering, plus the excluded doublets.
#
# Setup: Replace the /path/to/... placeholders below with local paths.
# The tab-delimited metadata must contain Sample and dir_cellranger_arc columns.
# Each dir_cellranger_arc entry must point to a Cell Ranger ARC sample directory
# containing outs/filtered_feature_bc_matrix.h5 and outs/atac_fragments.tsv.gz.
#
# Original environment:
# module load R/4.2.3 bedtools/2.30.0 hdf5/1.8.13 igvtools/2.4.0

# Session setup and packages --------------------------------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(42)
suppressPackageStartupMessages({
    library(ArchR)
    library(Seurat)
    library(Signac)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(pheatmap)
    library(parallel)
})

# Inputs and analysis parameters ---------------------------------------------

# Output directory
dir_out <- '/path/to/output/atac_doublets_ArchR'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Sample metadata
file_meta <- '/path/to/sample_metadata/meta_multiomeIBD_sample.txt'

# Analysis parameters
file_save <- 'multiomeIBD' # output file name
n_threads <- 16 # number of threads
genome <- 'hg38' # genome
bin_size <- 500 # genome bin size for count matrix
min_tss_score <- 3 # min TSS qc score
min_frags <- 2000 # min number of fragments
max_frags <- 1E5 # max number of fragments
doublet_filter_ratio <- 0.5 # filter ratio for ArchR doublet removal
param_lsi <- 2 # LSI method type: 1="tf-logidf", 2="log(tf-idf)" (Seurat default), 3="logtf-logidf"

# Sample metadata and input files --------------------------------------------

# Load sample metadata
meta_sample <- read.delim(file_meta, header = TRUE, sep = '\t')
name_sample <- meta_sample$Sample

# Restrict eligible barcodes to the filtered Cell Ranger ARC RNA matrix.
valid_barcodes <- list()
for (ix in 1:length(name_sample)) {
    file_data <- paste0(meta_sample$dir_cellranger_arc[ix], '/outs/filtered_feature_bc_matrix.h5')
    curr_sc <- Read10X_h5(file_data)
    valid_barcodes[[name_sample[ix]]] <- colnames(curr_sc[['Gene Expression']])
}

# Fragment file list
file_scatac_frag <- list()
for (ix in 1:length(name_sample)) {
    file_scatac_frag[[name_sample[ix]]] <- paste0(meta_sample$dir_cellranger_arc[ix], '/outs/atac_fragments.tsv.gz')
}

# Configure ArchR ------------------------------------------------------------
print('Set up ArchR')
addArchRThreads(threads = n_threads)
addArchRGenome(genome)

# Create Arrow files ---------------------------------------------------------
# For each sample, this step will:
# * Read accessible fragments from the provided input files
# * Calculate quality control information for each cell (ie, TSS enrichment scores and nucleosome info)
# * Filter cells based on quality control parameters
# * Create a genome-wide TileMatrix using specified bin size (default: 500 bp).
# * Create a GeneScoreMatrix using the annotation configured by addArchRGenome().
print('Create Arrow files')
file_arrow <- createArrowFiles(
    inputFiles = as.character(file_scatac_frag),
    sampleNames = names(file_scatac_frag),
    validBarcodes = valid_barcodes,
    minTSS = min_tss_score, minFrags = min_frags, maxFrags = max_frags,
    addTileMat = TRUE, addGeneScoreMat = TRUE,
    TileMatParams = list(tileSize = bin_size, excludeChr = c('chrM', 'chrY')),
    GeneScoreMatParams = list(tileSize = bin_size, excludeChr = c('chrM', 'chrY')),
    promoterRegion = c(2000, 500), QCDir = file.path(dir_out, 'QualityControl')
)

# Create ArchRProject
print('Create ArchRProject')
proj_ArchR <- ArchRProject(ArrowFiles = file_arrow, outputDirectory = dir_out, copyArrows = TRUE)

# Score putative doublets ----------------------------------------------------
print('Doublet detection')
dir_dd <- file.path(dir_out, 'DoubletDetection')
dir.create(dir_dd, showWarnings = FALSE)
proj_ArchR <- addDoubletScores(
    input = proj_ArchR,
    k = 10, # Number of neighbors used for doublet scoring
    knnMethod = 'UMAP', # embedding to use for nearest neighbor search
    LSIMethod = param_lsi,
    outDir = dir_dd
)

# Record all cells retained in the ArchR project before doublet filtering.
df_cells_raw <- data.frame(Name = proj_ArchR$cellNames)
df_cells_raw$Sample <- proj_ArchR$Sample
df_cells_raw$Cell <- gsub('^.*#', '', df_cells_raw$Name)
rownames(df_cells_raw) <- paste(df_cells_raw$Sample, df_cells_raw$Cell, sep = '_')
df_cells_raw <- df_cells_raw[, c('Sample', 'Cell')]

# Filter doublets
print('Filter doublets')
proj_ArchR <- filterDoublets(ArchRProj = proj_ArchR, filterRatio = doublet_filter_ratio)

# Record cells retained after doublet filtering.
df_cells_filtered <- data.frame(Name = proj_ArchR$cellNames)
df_cells_filtered$Sample <- proj_ArchR$Sample
df_cells_filtered$Cell <- gsub('^.*#', '', df_cells_filtered$Name)
rownames(df_cells_filtered) <- paste(df_cells_filtered$Sample, df_cells_filtered$Cell, sep = '_')
df_cells_filtered <- df_cells_filtered[, c('Sample', 'Cell')]

# Identify excluded doublets by comparing cell identifiers.
df_cells_doublet <- df_cells_raw[setdiff(rownames(df_cells_raw), rownames(df_cells_filtered)), ]

# Export cell lists ----------------------------------------------------------
print('Save cell info')
label_ratio <- gsub('[:.:]', 'p', doublet_filter_ratio)
file_out <- file.path(dir_out, paste0('cells_', file_save, '_qc_scatac_ArchR_raw.txt'))
write.table(df_cells_raw, file_out, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
file_out <- file.path(dir_out, paste0('cells_', file_save, '_qc_scatac_ArchR_filtered.txt'))
write.table(df_cells_filtered, file_out, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
file_out <- file.path(dir_out, paste0('cells_', file_save, '_qc_scatac_ArchR_doublets_ratio', label_ratio, '.txt'))
write.table(df_cells_doublet, file_out, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

# Clean up ArchR files -------------------------------------------------------
# Preserve the original cleanup of working-directory files and copied Arrows.
print('Clean up ArchR files')
for (ix in name_sample) {
    curr_file <- paste0(ix, '.arrow')
    system(paste0('rm ', curr_file))
}
system('rm -r ArchRLogs')
system('rm -r tmp')
system(paste0('rm -r ', file.path(dir_out, 'ArrowFiles')))

print('DONE')
