# IBD_multiome_data_gene_filter.R
#
# Purpose: Extract and filter Cell Ranger GTF gene annotations for downstream
# analyses. Export filtered gene bodies, gene lists, and an Ensembl-to-gene key.
#
# Setup: Replace each /path/to/... placeholder with the appropriate local path.
# The helper script must provide bedtools_sort() and bedtools_merge().
# Original environment:
# module load R/4.2.3 bedtools/2.30.0 hdf5/1.8.13 igvtools/2.4.0

# Session setup ---------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)
source('/path/to/helper_scripts/bedtools_utils.R')

# Inputs and filtering parameters --------------------------------------------

# Output directory
dir_out <- '/path/to/output/gene_filter'
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# Output file name
file_save <- 'refdata-cellranger-arc-GRCh38-2020-A-2.0.0'

# Gene annotation GTF file
file_gtf <- '/path/to/cellranger_arc_reference/genes/genes.gtf'

# Pseudogene list
file_pseudogene <- '/path/to/pseudogene_list/gene_pseudogene_hg38_Ensembl_v109_GRCh38_p13_processed.txt'

# Filtering parameters (preserved from the original analysis)
keep_chr <- paste0('chr', c(1:22, 'X', 'M')) # chromosomes to keep
keep_genetype <- c('protein_coding', 'lncRNA', 'TR_C_gene', 'IG_C_gene') # gene types to keep
del_gene_pattern <- '^HBA-|^HBB|^HBD|^HBG|^HBM|^HBQ1|^HBZ' # genes to delete from analysis

# Read and annotate genes ----------------------------------------------------

# Read GTF file
print('Read GTF file')
gtf <- read.delim(file_gtf, header = FALSE, skip = 5, sep = '\t')
gtf <- gtf[which(gtf[, 3] == 'gene'), ]
colnames(gtf) <- c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes')

# Load pseudogenes
gene_pseudo <- readLines(file_pseudogene)

# Function to extract attributes from GTF file
# Source: https://www.biostars.org/p/272889/
extract_attributes <- function(gtf_attributes, att_of_interest) {
    att <- strsplit(gtf_attributes, "; ")
    att <- gsub("\"", "", unlist(att))
    if (!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))) {
        return(unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
    } else {
        return(NA)
    }
}
gtf$gene_name <- unlist(lapply(gtf$attributes, extract_attributes, 'gene_name'))
gtf$gene_type <- unlist(lapply(gtf$attributes, extract_attributes, 'gene_type'))
gtf$id <- unlist(lapply(gtf$attributes, extract_attributes, 'gene_id'))

# Filter gene bodies ---------------------------------------------------------
# Retain the original GTF coordinates without conversion.
print('Create gene body bed file')
bed_all <- gtf[, c('chr', 'start', 'end', 'gene_name', 'gene_type')]
gene_y <- sort(unique(subset(bed_all, chr == 'chrY')[, 'gene_name'])) # chrY genes
gene_mt <- sort(unique(subset(bed_all, chr == 'chrM')[, 'gene_name'])) # mito genes
bed_all <- subset(bed_all, chr %in% keep_chr)
bed_all <- subset(bed_all, gene_type %in% keep_genetype)

# Remove pseudogenes and custom genes
bed_all <- subset(bed_all, gene_name %in% setdiff(bed_all$gene_name, gene_pseudo))
bed_all <- subset(bed_all, gene_name %in% setdiff(bed_all$gene_name,
                                        bed_all$gene_name[which(grepl(del_gene_pattern, bed_all$gene_name))]))
bed_all <- bedtools_sort(bed_all)

# Merge duplicate gene regions -----------------------------------------------
print('Merge duplicate gene regions')
dup_gene <- bed_all$gene_name[which(duplicated(bed_all$gene_name))]
dup_bed <- subset(bed_all, gene_name %in% dup_gene)
for (ix in dup_gene) {
    curr_bed <- subset(bed_all, gene_name == ix)
    add_bed <- bedtools_merge(curr_bed[, 1:3])
    add_bed$gene_name <- curr_bed$gene_name[1]
    add_bed$gene_type <- curr_bed$gene_type[1]
    bed_all <- subset(bed_all, gene_name != ix)
    bed_all <- rbind(bed_all, add_bed)
}
bed_all <- bedtools_sort(bed_all)

# Export results -------------------------------------------------------------
print('Save gene sets')
# Save gene body bed file
file_out <- file.path(dir_out, paste0('gene_body_', file_save, '_filtered.bed'))
write.table(bed_all, file_out, quote = FALSE, sep = '\t', col.names = FALSE, row.names = FALSE)

# Save all filtered genes
file_out <- file.path(dir_out, paste0('gene_set_', file_save, '_filtered.txt'))
writeLines(sort(unique(bed_all$gene_name)), file_out)

# Save chromosome Y and mitochondrial genes from the unfiltered annotation.
file_out <- file.path(dir_out, paste0('gene_set_', file_save, '_chrY.txt'))
writeLines(gene_y, file_out)
file_out <- file.path(dir_out, paste0('gene_set_', file_save, '_chrM.txt'))
writeLines(gene_mt, file_out)

# Save the gene-to-Ensembl key from the unfiltered annotation.
df_id <- gtf[, c('id', 'gene_name')]
colnames(df_id) <- c('ID', 'Gene')
file_out <- file.path(dir_out, paste0('gene_key_', file_save, '_Ensembl.txt'))
write.table(df_id, file_out, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')

print('DONE')
