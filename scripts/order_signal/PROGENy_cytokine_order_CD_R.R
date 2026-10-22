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
file_df_gsea_PROGENy<-file.path(results_dir,paste('Downsample_statWald/PROGENy_CD64_ImmDict/CD_R/PROGENy_CD64_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_df_PROGENy_RNAEnrich<-file.path(results_dir,paste('Downsample_RNAEnrich/PROGENy_CD64_ImmDict/CD_R/PROGENy_CD64_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_keep_context<-file.path(data_dir,'keep_context/keep_context_CD_R_Nsample_3_min_filter.txt')
data_assay <- 'RNA' # data assay
#####
# Params
file_save <- 'IBD_multiome' # file save base
# output directories
dir_out <- file.path(results_dir,'Downsampled_fgsea_RNAEnrich_co_viz/CD_R')
dir.create(dir_out, showWarnings=FALSE)
keep_context<-readLines(file_keep_context)
print('load enrichment analysis')
df<-read.table(file_df_gsea_PROGENy,sep='\t',header=T,check.names=F)
df<-df[df$Context%in%keep_context,]
sig_sets_PROGENy<-unique(df[abs(df$`-log10Padj`)>1,'Pathway'])
df_PROGENy_RNAEnrich<-read.table(file_df_PROGENy_RNAEnrich,sep='\t',header=T,check.names=F)
df_PROGENy_RNAEnrich<-df_PROGENy_RNAEnrich[df_PROGENy_RNAEnrich$Context%in%keep_context,]
sig_sets_PROGENy_RNAEnrich<-unique(df_PROGENy_RNAEnrich[abs(df_PROGENy_RNAEnrich$`-log10Padj`)>1,'Pathway'])
sig_sets_PROGENy<-intersect(sig_sets_PROGENy,sig_sets_PROGENy_RNAEnrich)
df<-df[df$Pathway%in%sig_sets_PROGENy,]
df_max_padj_pathway<-NULL
for (ix in sig_sets_PROGENy){
  curr_df<-df[df$Pathway==ix,]
  curr_abs_log10Padj_max<-max(abs(curr_df$`-log10Padj`))
  curr_Nsig_context<-length(unique(curr_df[abs(curr_df$`-log10Padj`)>=1,'Context']))
  curr_df_max_padj_pathway<-data.frame(Pathway=ix,`max_log10Padj`=curr_abs_log10Padj_max,
                                               Nsig_Context=curr_Nsig_context)
  df_max_padj_pathway<-rbind(df_max_padj_pathway,curr_df_max_padj_pathway)                                        
}
df_max_padj_pathway<-df_max_padj_pathway[rev(order(df_max_padj_pathway$max_log10Padj)),]
file_out<-file.path(dir_out,paste0('sig_sets_PROGENy_cytokine_order_by_max_log10Padj.txt'))
write.table(df_max_padj_pathway,file_out,sep='\t',col.names=T,row.names=F,quote=F)