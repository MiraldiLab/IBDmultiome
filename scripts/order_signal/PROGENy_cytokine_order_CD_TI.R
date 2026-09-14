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
file_df_gsea_PROGENy<-file.path(results_dir,paste('Downsample_statWald/PROGENy_CD64_ImmDict/CD_TI/PROGENy_CD64_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_df_PROGENy_RNAEnrich<-file.path(results_dir,paste('Downsample_RNAEnrich/PROGENy_CD64_ImmDict/CD_TI/PROGENy_CD64_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_keep_context<-file.path(data_dir,'keep_context/keep_context_CD_TI_Nsample_3_min_filter.txt')
file_PROGENy_order<-file.path(results_dir,'Downsampled_fgsea_RNAEnrich_co_viz/CD_R/sig_sets_PROGENy_cytokine_order_by_max_log10Padj.txt')
data_assay <- 'RNA' # data assay
#####
# Params
file_save <- 'IBD_multiome' # file save base
# output directories
dir_out <- file.path(results_dir,'Downsampled_fgsea_RNAEnrich_co_viz/CD_TI')
dir.create(dir_out, showWarnings=FALSE)
keep_context<-readLines(file_keep_context)
print('load enrichment analysis')
df_gsea_PROGENy<-read.table(file_df_gsea_PROGENy,sep='\t',header=T,check.names=F)
df_gsea_PROGENy<-df_gsea_PROGENy[df_gsea_PROGENy$Context%in%keep_context,]
sig_sets_PROGENy<-unique(df_gsea_PROGENy[abs(df_gsea_PROGENy$`-log10Padj`)>1,'Pathway'])

df_PROGENy_RNAEnrich<-read.table(file_df_PROGENy_RNAEnrich,sep='\t',header=T,check.names=F)
df_PROGENy_RNAEnrich<-df_PROGENy_RNAEnrich[df_PROGENy_RNAEnrich$Context%in%keep_context,]
sig_sets_PROGENy_RNAEnrich<-unique(df_PROGENy_RNAEnrich[abs(df_PROGENy_RNAEnrich$`-log10Padj`)>1,'Pathway'])
sig_sets_PROGENy<-intersect(sig_sets_PROGENy,sig_sets_PROGENy_RNAEnrich)
#signal order in CD, rectum
PROGENy_order<-read.table(file_PROGENy_order,sep='\t',header=T)
PROGENy_order<-PROGENy_order$Pathway
PROGENy_order<-intersect(PROGENy_order,sig_sets_PROGENy)
#add signals if missed from CD, rectum
PROGENy_order<-c(PROGENy_order,setdiff(sig_sets_PROGENy,PROGENy_order))
file_out<-file.path(dir_out,paste0('sig_sets_PROGENy_cytokine_order_by_max_log10Padj.txt'))
write.table(df_max_padj_pathway,file_out,sep='\t',col.names=T,row.names=F,quote=F)