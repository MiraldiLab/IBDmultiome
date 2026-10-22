# IBD_multiome_process_atac_fragments.R
#
# Purpose: Retain ATAC fragments from QC-passing cells, remove fragments
# overlapping blacklist regions, and split each retained fragment into two
# one-base cut-site intervals. Sort, compress, and index each sample's output.
#
# Setup: Replace all /path/to/... placeholders with local paths.
# Required sample metadata columns: Sample and dir_cellranger_arc.
# Each Cell Ranger ARC directory must contain outs/atac_fragments.tsv.gz.
# The QC cell table must contain Sample and barcode columns.
# The blacklist is a BED file. The helper script must provide bedtools_subtract()
# and bedtools_sort(); BEDTools, bgzip, and tabix must be available on PATH.
# n_skip_frag specifies the number of header lines skipped in each input file.
#
# Original environment:
# module load R/4.2.3 bedtools/2.30.0 hdf5/1.8.13 igvtools/2.4.0

# Session setup and packages --------------------------------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(42)
suppressPackageStartupMessages({
    library(magrittr)
    library(foreach)
    library(doParallel)
})
source('/path/to/helper_scripts/bedtools_utils.R')

# Inputs and processing parameters -------------------------------------------

# Output directory
dir_out <- '/path/to/output/atac_cut_sites'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Sample metadata
file_meta <- '/path/to/sample_metadata/meta_multiomeIBD_sample.txt'

# Cells passing QC
file_cells_qc <- '/path/to/qc_metrics/meta_multiomeIBD_filtered_scatac.txt'

# Blacklist bed file
file_blacklist <- '/path/to/blacklist/blacklist_hg38_ENCODE_maxATAC_extended.bed'

# Processing parameters
n_skip_frag <- 52 # lines to skip in the fragment file (header)
keep_chr <- paste0('chr', c(1:22, 'X', 'Y', 'M')) # chromosomes to keep

# Load sample metadata, QC cell lists, and blacklist -------------------------

# Load sample metadata
meta_sample <- read.delim(file_meta, header = TRUE, sep = '\t')
name_sample <- meta_sample$Sample

# Load cells passing QC
cells_qc <- read.delim(file_cells_qc, header = TRUE, sep = '\t')

# Load blacklist
bed_blacklist <- read.delim(file_blacklist, header = FALSE, sep = '\t')

# Process each sample --------------------------------------------------------
for (ix in 1:length(name_sample)) {

    # Load original fragment bed
    print(paste0('Process fragments: ', name_sample[ix]))
    curr_file <- paste0(meta_sample$dir_cellranger_arc[ix], '/outs/atac_fragments.tsv.gz')
    frag_old <- read.delim(gzfile(curr_file), header = FALSE, sep = '\t', skip = n_skip_frag)
    colnames(frag_old) <- c('Chr', 'Start', 'End', 'Barcode', 'Count')

    # Filter quality cells
    curr_cells <- subset(cells_qc, Sample == name_sample[ix])$barcode
    frag_old <- subset(frag_old, Barcode %in% curr_cells & Chr %in% keep_chr)

    # Remove entire fragments that overlap blacklist regions (subtract -A).
    print(' - Remove blacklist')
    frag_old <- bedtools_subtract(frag_old, bed_blacklist)
    colnames(frag_old) <- c('Chr', 'Start', 'End', 'Barcode', 'Count')

    # Split each fragment into [Start, Start + 1) and [End - 1, End) intervals.
    # Preserve the original barcode and count on both entries; apply no shift.
    frag_left <- data.frame(Chr = frag_old$Chr, Start = frag_old$Start, End = frag_old$Start + 1,
                            Barcode = frag_old$Barcode, Count = frag_old$Count)
    frag_right <- data.frame(Chr = frag_old$Chr, Start = frag_old$End - 1, End = frag_old$End,
                            Barcode = frag_old$Barcode, Count = frag_old$Count)
    frag_new <- rbind(frag_left, frag_right)

    # Sort the combined cut-site intervals.
    print(' - Sort')
    frag_new <- bedtools_sort(frag_new)

    # Export cut-site intervals using the original fragments_<Sample> basename.
    print(' - Save new fragment bed')
    file_out <- file.path(dir_out, paste0('fragments_', name_sample[ix], '.tsv'))
    write.table(frag_new, file_out, quote = FALSE, sep = '\t', col.names = FALSE, row.names = FALSE)

    # Compress and index the cut-site file.
    # The data is block-gzipped to allow indexing and to save disk space
    # The tabix index is created with --preset=bed
    print(' - Zip')
    system(paste0('bgzip -f ', file_out, ' > ', gsub('.tsv', '.tsv.gz', file_out)))
    print(' - Index')
    system(paste0('tabix -p bed ', gsub('.tsv', '.tsv.gz', file_out)))

}

print('DONE')
