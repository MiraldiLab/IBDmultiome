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
# input date object
file_df_gsea_CytoSig<-file.path(results_dir,paste('Downsample_statWald/CytoSig/CD_R/CytoSig_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_df_CytoSig_RNAEnrich<-file.path(results_dir,paste('Downsample_RNAEnrich/CytoSig/CD_R/CytoSig_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_df_cytokine_family<-file.path(data_dir,"gene_sets","CytoSig","CytoSig_cytokine_family_class.txt")
file_keep_context<-file.path(data_dir,"keep_context","keep_context_CD_R_Nsample_3_min_filter.txt")
data_assay <- 'RNA' # data assay
order_disease_tissue<-c('CD_R')
order_celltype<-c('CD4_T_cell','CD4_T_Eff','CD8_T_cell','Naive_T_cell','gdT_cell','NK_cell',
                  'B_cell','Plasma_cell','Macrophage','DC','Fibroblast','BEST4_Colonocyte','Colonocyte','Immature_Colonocyte',
                  'Enterocyte','Immature_Enterocyte','Enteroendocrine','Endothelial','Goblet','Immature_Goblet','TA','Cycling_TA','Stem','Tuft','Paneth','M_cell')

#####
# Params
file_save <- 'IBD_multiome' # file save base
# output directories
dir_out <- file.path(results_dir,'Downsampled_fgsea_RNAEnrich_co_viz/CD_R')
dir.create(dir_out, showWarnings=FALSE)
order_clustersubtype<-NULL
for (ix in order_celltype){
  for (jx in order_disease_tissue){
    order_clustersubtype<-c(order_clustersubtype,paste(ix,jx,sep='_'))
  }
}
keep_context<-readLines(file_keep_context)
order_clustersubtype<-intersect(order_clustersubtype,keep_context)
print('load enrichment analysis')
df<-read.table(file_df_gsea_CytoSig,sep='\t',header=T,check.names=F)
df<-df[df$Context%in%order_clustersubtype,]
sig_sets_CytoSig<-unique(df[abs(df$`-log10Padj`)>1,'Pathway'])
df_CytoSig_RNAEnrich<-read.table(file_CytoSig_RNAEnrich,sep='\t',header=T,check.names=F)
df_CytoSig_RNAEnrich<-df_CytoSig_RNAEnrich[df_CytoSig_RNAEnrich$Context%in%keep_context,]
sig_sets_CytoSig_RNAEnrich<-unique(df_CytoSig_RNAEnrich[abs(df_CytoSig_RNAEnrich$`-log10Padj`)>1,'Pathway'])
sig_sets_CytoSig<-intersect(sig_sets_CytoSig,sig_sets_CytoSig_RNAEnrich)
df<-df[df$Pathway%in%sig_sets_CytoSig,]
df_cytokine_family<-read.table(file_df_cytokine_family,sep='\t',header=T)
df<-merge(df,df_cytokine_family,by='Pathway',all.x=T)
df_max_padj_cytokine_family<-NULL
for (ix in sig_sets_CytoSig){
  curr_df<-df[df$Pathway==ix,]
  curr_family<-unique(curr_df$Family)
  curr_abs_log10Padj_max<-max(abs(curr_df$`-log10Padj`))
  curr_Nsig_context<-length(unique(curr_df[abs(curr_df$`-log10Padj`)>=1,'Context']))
  curr_df_max_padj_cytokine_family<-data.frame(Pathway=ix,
                                               Family=curr_family,`max_log10Padj`=curr_abs_log10Padj_max,
                                               Nsig_Context=curr_Nsig_context)
  df_max_padj_cytokine_family<-rbind(df_max_padj_cytokine_family,curr_df_max_padj_cytokine_family)                                        
}
df_max_padj_cytokine_family<-df_max_padj_cytokine_family[rev(order(df_max_padj_cytokine_family$max_log10Padj)),]
order_cytokine_family<-readLines(file.path(dir_out,paste0('sig_sets_CytoSig_cytokine_family_order.txt')))
df_max_padj_cytokine_family <- df_max_padj_cytokine_family %>%
  mutate(Family = factor(Family, levels = order_cytokine_family)) %>%
  arrange(Family, desc(max_log10Padj))
file_out<-file.path(dir_out,paste0('sig_sets_CytoSig_cytokine_family_order_by_max_log10Padj.txt'))
write.table(df_max_padj_cytokine_family,file_out,sep='\t',col.names=T,row.names=F,quote=F)