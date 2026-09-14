rm(list=ls())
options(stringsAsFactors=FALSE)
set.seed(42)
suppressPackageStartupMessages({
library(Seurat)
library(ggplot2)
library(magrittr)
library(sva)
library(fgsea)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(matrixStats)
library(dplyr)
library(patchwork)
})
project_dir <- "."
script_dir <- file.path(project_dir, "scripts")
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results","GSEA_RNA_Status")
file_df_gsea_CytoSig<-file.path(results_dir,paste('Downsample_statWald/CytoSig/UC_R/CytoSig_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_df_CytoSig_RNAEnrich<-file.path(results_dir,paste('Downsample_RNAEnrich/CytoSig/UC_R/CytoSig_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_df_cytokine_family<-file.path(data_dir,"gene_sets","CytoSig","CytoSig_cytokine_family_class.txt")
file_df_cytokine_order<-file.path(results_dir,'Downsampled_fgsea_RNAEnrich_co_viz/CD_R/sig_sets_CytoSig_cytokine_family_order_by_max_log10Padj.txt')
file_keep_context<-file.path(data_dir,"keep_context","keep_context_UC_R_Nsample_3_min_filter.txt")
data_assay <- 'RNA' # data assay
order_disease_tissue<-c('UC_R')
#####
# Params
file_save <- 'IBD_multiome' # file save base
# output directories
dir_out <- file.path(results_dir,'Downsampled_fgsea_RNAEnrich_co_viz/UC_R')
dir.create(dir_out, showWarnings=FALSE)
keep_context<-readLines(file_keep_context)
print('load enrichment analysis')
df<-read.table(file_df_gsea_CytoSig,sep='\t',header=T,check.names=F)
df<-df[df$Context%in%keep_context,]
sig_sets_CytoSig<-unique(df[abs(df$`-log10Padj`)>1,'Pathway'])
df_CytoSig_RNAEnrich<-read.table(file_df_CytoSig_RNAEnrich,sep='\t',header=T,check.names=F)
df_CytoSig_RNAEnrich<-df_CytoSig_RNAEnrich[df_CytoSig_RNAEnrich$Context%in%keep_context,]
sig_sets_CytoSig_RNAEnrich<-unique(df_CytoSig_RNAEnrich[abs(df_CytoSig_RNAEnrich$`-log10Padj`)>1,'Pathway'])
sig_sets_CytoSig<-intersect(sig_sets_CytoSig,sig_sets_CytoSig_RNAEnrich)
df<-df[df$Pathway%in%sig_sets_CytoSig,]
#cytokine order from CD_R
df_cytokine_order<-read.table(file_df_cytokine_order,sep='\t',header=T)
df_cytokine_order<-df_cytokine_order[,c('Pathway','Family')]
df_cytokine_order<-df_cytokine_order[df_cytokine_order$Pathway%in%sig_sets_CytoSig,]
extra_cytokine<-setdiff(sig_sets_CytoSig,df_cytokine_order$Pathway)
#missing cytokine added to the related family
df_cytokine_family<-read.table(file_df_cytokine_family,sep='\t',header=T)
df_cytokine_family<-df_cytokine_family[df_cytokine_family$Pathway%in%extra_cytokine,]
for (fam in unique(df_cytokine_family$Family)) {

  add_rows <- df_cytokine_family %>% filter(Family == fam)

  last_idx <- max(which(df_cytokine_order$Family == fam))

  df_cytokine_order <- bind_rows(
    df_cytokine_order[1:last_idx, ],
    add_rows,
    df_cytokine_order[(last_idx + 1):nrow(df_cytokine_order), ]
  )
} 
file_out<-file.path(dir_out,paste0('sig_sets_CytoSig_cytokine_family_order.txt'))
write.table(df_cytokine_order,file_out,sep='\t',col.names=T,row.names=F,quote=F)