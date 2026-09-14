rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(42)
suppressPackageStartupMessages({
  library(ggplot2)
  library(magrittr)
  library(sva)
  library(RColorBrewer)
  library(scales)
  library(dplyr)
  library(reshape2)
  library(DESeq2)
  library(ggrepel)
  library(scales)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(tidyverse)
  library(glue)
  library(fs)
  library(colorspace)
})
project_dir <- "."
script_dir <- file.path(project_dir, "scripts")
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results","GSEA_RNA_Status")
dir_out<-file.path(results_dir,'Downsampled_fgsea_RNAEnrich_co_viz/PROGENy_CD64_ImmDict/ScatterPlot_comparison')
file_gsea_res_ImmDict_CD_R<-file.path(results_dir,'Downsample_statWald/PROGENy_CD64_ImmDict/CD_R/ImmDict_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt')
file_gsea_res_ImmDict_CD_TI<-file.path(results_dir,'Downsample_statWald/PROGENy_CD64_ImmDict/CD_TI/ImmDict_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt')
file_gsea_res_ImmDict_UC_R<-file.path(results_dir,'Downsample_statWald/PROGENy_CD64_ImmDict/UC_R/ImmDict_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt')

file_gsea_res_ImmDict_CD_R_RNAEnrich<-file.path(results_dir,'Downsample_RNAEnrich/PROGENy_CD64_ImmDict/CD_R/ImmDict_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt')
file_gsea_res_ImmDict_CD_TI_RNAEnrich<-file.path(results_dir,'Downsample_RNAEnrich/PROGENy_CD64_ImmDict/CD_TI/ImmDict_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt')
file_gsea_res_ImmDict_UC_R_RNAEnrich<-file.path(results_dir,'Downsample_RNAEnrich/PROGENy_CD64_ImmDict/UC_R/ImmDict_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt')

file_gsea_res_PROGENy_CD_R<-file.path(results_dir,'Downsample_statWald/PROGENy_CD64_ImmDict/CD_R/PROGENy_CD64_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt')
file_gsea_res_PROGENy_CD_TI<-file.path(results_dir,'Downsample_statWald/PROGENy_CD64_ImmDict/CD_TI/PROGENy_CD64_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt')
file_gsea_res_PROGENy_UC_R<-file.path(results_dir,'Downsample_statWald/PROGENy_CD64_ImmDict/UC_R/PROGENy_CD64_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt')

file_gsea_res_PROGENy_CD_R_RNAEnrich<-file.path(results_dir,'Downsample_RNAEnrich/PROGENy_CD64_ImmDict/CD_R/PROGENy_CD64_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt')
file_gsea_res_PROGENy_CD_TI_RNAEnrich<-file.path(results_dir,'Downsample_RNAEnrich/PROGENy_CD64_ImmDict/CD_TI/PROGENy_CD64_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt')
file_gsea_res_PROGENy_UC_R_RNAEnrich<-file.path(results_dir,'Downsample_RNAEnrich/PROGENy_CD64_ImmDict/UC_R/PROGENy_CD64_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt')
#load results ImmDict
gsea_res_ImmDict_CD_R<-read.table(file_gsea_res_ImmDict_CD_R,sep='\t',header=T,check.names=F)
gsea_res_ImmDict_CD_TI<-read.table(file_gsea_res_ImmDict_CD_TI,sep='\t',header=T,check.names=F)
gsea_res_ImmDict_UC_R<-read.table(file_gsea_res_ImmDict_UC_R,sep='\t',header=T,check.names=F)
gsea_res_ImmDict_CD_R_RNAEnrich<-read.table(file_gsea_res_ImmDict_CD_R_RNAEnrich,sep='\t',header=T,check.names=F)
gsea_res_ImmDict_CD_TI_RNAEnrich<-read.table(file_gsea_res_ImmDict_CD_TI_RNAEnrich,sep='\t',header=T,check.names=F)
gsea_res_ImmDict_UC_R_RNAEnrich<-read.table(file_gsea_res_ImmDict_UC_R_RNAEnrich,sep='\t',header=T,check.names=F)
#modify the fgsea results based on RNAEnrich CD R
gsea_res_ImmDict_CD_R$ContextPathway<-paste(gsea_res_ImmDict_CD_R$Context,gsea_res_ImmDict_CD_R$Pathway,sep='_')
gsea_res_ImmDict_CD_R_RNAEnrich$ContextPathway<-paste(gsea_res_ImmDict_CD_R_RNAEnrich$Context,gsea_res_ImmDict_CD_R_RNAEnrich$Pathway,sep='_')
sig_up_context_pathway_fgsea_CD_R<-gsea_res_ImmDict_CD_R[gsea_res_ImmDict_CD_R$`-log10Padj`>1,'ContextPathway']
sig_dn_context_pathway_fgsea_CD_R<-gsea_res_ImmDict_CD_R[gsea_res_ImmDict_CD_R$`-log10Padj`<(-1),'ContextPathway']
sig_up_context_pathway_RNAEnrich_CD_R<-gsea_res_ImmDict_CD_R_RNAEnrich[gsea_res_ImmDict_CD_R_RNAEnrich$`-log10Padj`>1,'ContextPathway']
sig_dn_context_pathway_RNAEnrich_CD_R<-gsea_res_ImmDict_CD_R_RNAEnrich[gsea_res_ImmDict_CD_R_RNAEnrich$`-log10Padj`<(-1),'ContextPathway']
sig_up_disagree_CD_R<-setdiff(sig_up_context_pathway_fgsea_CD_R,sig_up_context_pathway_RNAEnrich_CD_R)
sig_dn_disagree_CD_R<-setdiff(sig_dn_context_pathway_fgsea_CD_R,sig_dn_context_pathway_RNAEnrich_CD_R)
sig_disagree_CD_R<-c(sig_up_disagree_CD_R,sig_dn_disagree_CD_R)
gsea_res_ImmDict_CD_R[gsea_res_ImmDict_CD_R$ContextPathway%in%sig_disagree_CD_R,'-log10Padj']<-0
#modify the fgsea results based on RNAEnrich CD TI
gsea_res_ImmDict_CD_TI$ContextPathway<-paste(gsea_res_ImmDict_CD_TI$Context,gsea_res_ImmDict_CD_TI$Pathway,sep='_')
gsea_res_ImmDict_CD_TI_RNAEnrich$ContextPathway<-paste(gsea_res_ImmDict_CD_TI_RNAEnrich$Context,gsea_res_ImmDict_CD_TI_RNAEnrich$Pathway,sep='_')
sig_up_context_pathway_fgsea_CD_TI<-gsea_res_ImmDict_CD_TI[gsea_res_ImmDict_CD_TI$`-log10Padj`>1,'ContextPathway']
sig_dn_context_pathway_fgsea_CD_TI<-gsea_res_ImmDict_CD_TI[gsea_res_ImmDict_CD_TI$`-log10Padj`<(-1),'ContextPathway']
sig_up_context_pathway_RNAEnrich_CD_TI<-gsea_res_ImmDict_CD_TI_RNAEnrich[gsea_res_ImmDict_CD_TI_RNAEnrich$`-log10Padj`>1,'ContextPathway']
sig_dn_context_pathway_RNAEnrich_CD_TI<-gsea_res_ImmDict_CD_TI_RNAEnrich[gsea_res_ImmDict_CD_TI_RNAEnrich$`-log10Padj`<(-1),'ContextPathway']
sig_up_disagree_CD_TI<-setdiff(sig_up_context_pathway_fgsea_CD_TI,sig_up_context_pathway_RNAEnrich_CD_TI)
sig_dn_disagree_CD_TI<-setdiff(sig_dn_context_pathway_fgsea_CD_TI,sig_dn_context_pathway_RNAEnrich_CD_TI)
sig_disagree_CD_TI<-c(sig_up_disagree_CD_TI,sig_dn_disagree_CD_TI)
gsea_res_ImmDict_CD_TI[gsea_res_ImmDict_CD_TI$ContextPathway%in%sig_disagree_CD_TI,'-log10Padj']<-0
#modify the fgsea results based on RNAEnrich CD R
gsea_res_ImmDict_UC_R$ContextPathway<-paste(gsea_res_ImmDict_UC_R$Context,gsea_res_ImmDict_UC_R$Pathway,sep='_')
gsea_res_ImmDict_UC_R_RNAEnrich$ContextPathway<-paste(gsea_res_ImmDict_UC_R_RNAEnrich$Context,gsea_res_ImmDict_UC_R_RNAEnrich$Pathway,sep='_')
sig_up_context_pathway_fgsea_UC_R<-gsea_res_ImmDict_UC_R[gsea_res_ImmDict_UC_R$`-log10Padj`>1,'ContextPathway']
sig_dn_context_pathway_fgsea_UC_R<-gsea_res_ImmDict_UC_R[gsea_res_ImmDict_UC_R$`-log10Padj`<(-1),'ContextPathway']
sig_up_context_pathway_RNAEnrich_UC_R<-gsea_res_ImmDict_UC_R_RNAEnrich[gsea_res_ImmDict_UC_R_RNAEnrich$`-log10Padj`>1,'ContextPathway']
sig_dn_context_pathway_RNAEnrich_UC_R<-gsea_res_ImmDict_UC_R_RNAEnrich[gsea_res_ImmDict_UC_R_RNAEnrich$`-log10Padj`<(-1),'ContextPathway']
sig_up_disagree_UC_R<-setdiff(sig_up_context_pathway_fgsea_UC_R,sig_up_context_pathway_RNAEnrich_UC_R)
sig_dn_disagree_UC_R<-setdiff(sig_dn_context_pathway_fgsea_UC_R,sig_dn_context_pathway_RNAEnrich_UC_R)
sig_disagree_UC_R<-c(sig_up_disagree_UC_R,sig_dn_disagree_UC_R)
gsea_res_ImmDict_UC_R[gsea_res_ImmDict_UC_R$ContextPathway%in%sig_disagree_UC_R,'-log10Padj']<-0
#load ImmDict results
gsea_res_PROGENy_CD_R<-read.table(file_gsea_res_PROGENy_CD_R,sep='\t',header=T,check.names=F)
gsea_res_PROGENy_CD_TI<-read.table(file_gsea_res_PROGENy_CD_TI,sep='\t',header=T,check.names=F)
gsea_res_PROGENy_UC_R<-read.table(file_gsea_res_PROGENy_UC_R,sep='\t',header=T,check.names=F)
gsea_res_PROGENy_CD_R_RNAEnrich<-read.table(file_gsea_res_PROGENy_CD_R_RNAEnrich,sep='\t',header=T,check.names=F)
gsea_res_PROGENy_CD_TI_RNAEnrich<-read.table(file_gsea_res_PROGENy_CD_TI_RNAEnrich,sep='\t',header=T,check.names=F)
gsea_res_PROGENy_UC_R_RNAEnrich<-read.table(file_gsea_res_PROGENy_UC_R_RNAEnrich,sep='\t',header=T,check.names=F)
#modify the fgsea results based on RNAEnrich CD R
gsea_res_PROGENy_CD_R$ContextPathway<-paste(gsea_res_PROGENy_CD_R$Context,gsea_res_PROGENy_CD_R$Pathway,sep='_')
gsea_res_PROGENy_CD_R_RNAEnrich$ContextPathway<-paste(gsea_res_PROGENy_CD_R_RNAEnrich$Context,gsea_res_PROGENy_CD_R_RNAEnrich$Pathway,sep='_')
sig_up_context_pathway_fgsea_CD_R<-gsea_res_PROGENy_CD_R[gsea_res_PROGENy_CD_R$`-log10Padj`>1,'ContextPathway']
sig_dn_context_pathway_fgsea_CD_R<-gsea_res_PROGENy_CD_R[gsea_res_PROGENy_CD_R$`-log10Padj`<(-1),'ContextPathway']
sig_up_context_pathway_RNAEnrich_CD_R<-gsea_res_PROGENy_CD_R_RNAEnrich[gsea_res_PROGENy_CD_R_RNAEnrich$`-log10Padj`>1,'ContextPathway']
sig_dn_context_pathway_RNAEnrich_CD_R<-gsea_res_PROGENy_CD_R_RNAEnrich[gsea_res_PROGENy_CD_R_RNAEnrich$`-log10Padj`<(-1),'ContextPathway']
sig_up_disagree_CD_R<-setdiff(sig_up_context_pathway_fgsea_CD_R,sig_up_context_pathway_RNAEnrich_CD_R)
sig_dn_disagree_CD_R<-setdiff(sig_dn_context_pathway_fgsea_CD_R,sig_dn_context_pathway_RNAEnrich_CD_R)
sig_disagree_CD_R<-c(sig_up_disagree_CD_R,sig_dn_disagree_CD_R)
gsea_res_PROGENy_CD_R[gsea_res_PROGENy_CD_R$ContextPathway%in%sig_disagree_CD_R,'-log10Padj']<-0
#modify the fgsea results based on RNAEnrich CD TI
gsea_res_PROGENy_CD_TI$ContextPathway<-paste(gsea_res_PROGENy_CD_TI$Context,gsea_res_PROGENy_CD_TI$Pathway,sep='_')
gsea_res_PROGENy_CD_TI_RNAEnrich$ContextPathway<-paste(gsea_res_PROGENy_CD_TI_RNAEnrich$Context,gsea_res_PROGENy_CD_TI_RNAEnrich$Pathway,sep='_')
sig_up_context_pathway_fgsea_CD_TI<-gsea_res_PROGENy_CD_TI[gsea_res_PROGENy_CD_TI$`-log10Padj`>1,'ContextPathway']
sig_dn_context_pathway_fgsea_CD_TI<-gsea_res_PROGENy_CD_TI[gsea_res_PROGENy_CD_TI$`-log10Padj`<(-1),'ContextPathway']
sig_up_context_pathway_RNAEnrich_CD_TI<-gsea_res_PROGENy_CD_TI_RNAEnrich[gsea_res_PROGENy_CD_TI_RNAEnrich$`-log10Padj`>1,'ContextPathway']
sig_dn_context_pathway_RNAEnrich_CD_TI<-gsea_res_PROGENy_CD_TI_RNAEnrich[gsea_res_PROGENy_CD_TI_RNAEnrich$`-log10Padj`<(-1),'ContextPathway']
sig_up_disagree_CD_TI<-setdiff(sig_up_context_pathway_fgsea_CD_TI,sig_up_context_pathway_RNAEnrich_CD_TI)
sig_dn_disagree_CD_TI<-setdiff(sig_dn_context_pathway_fgsea_CD_TI,sig_dn_context_pathway_RNAEnrich_CD_TI)
sig_disagree_CD_TI<-c(sig_up_disagree_CD_TI,sig_dn_disagree_CD_TI)
gsea_res_PROGENy_CD_TI[gsea_res_PROGENy_CD_TI$ContextPathway%in%sig_disagree_CD_TI,'-log10Padj']<-0
#modify the fgsea results based on RNAEnrich UC R
gsea_res_PROGENy_UC_R$ContextPathway<-paste(gsea_res_PROGENy_UC_R$Context,gsea_res_PROGENy_UC_R$Pathway,sep='_')
gsea_res_PROGENy_UC_R_RNAEnrich$ContextPathway<-paste(gsea_res_PROGENy_UC_R_RNAEnrich$Context,gsea_res_PROGENy_UC_R_RNAEnrich$Pathway,sep='_')
sig_up_context_pathway_fgsea_UC_R<-gsea_res_PROGENy_UC_R[gsea_res_PROGENy_UC_R$`-log10Padj`>1,'ContextPathway']
sig_dn_context_pathway_fgsea_UC_R<-gsea_res_PROGENy_UC_R[gsea_res_PROGENy_UC_R$`-log10Padj`<(-1),'ContextPathway']
sig_up_context_pathway_RNAEnrich_UC_R<-gsea_res_PROGENy_UC_R_RNAEnrich[gsea_res_PROGENy_UC_R_RNAEnrich$`-log10Padj`>1,'ContextPathway']
sig_dn_context_pathway_RNAEnrich_UC_R<-gsea_res_PROGENy_UC_R_RNAEnrich[gsea_res_PROGENy_UC_R_RNAEnrich$`-log10Padj`<(-1),'ContextPathway']
sig_up_disagree_UC_R<-setdiff(sig_up_context_pathway_fgsea_UC_R,sig_up_context_pathway_RNAEnrich_UC_R)
sig_dn_disagree_UC_R<-setdiff(sig_dn_context_pathway_fgsea_UC_R,sig_dn_context_pathway_RNAEnrich_UC_R)
sig_disagree_UC_R<-c(sig_up_disagree_UC_R,sig_dn_disagree_UC_R)
gsea_res_PROGENy_UC_R[gsea_res_PROGENy_UC_R$ContextPathway%in%sig_disagree_UC_R,'-log10Padj']<-0
#RPOGENy intraceullar pathway without clear ligand receptor information set as 1.
file_df_PROGENy_LR<-file.path(data_dir, "gene_sets", "PROGENy_CD64", "PROGENy_CD64_LR_pairs_by_cellchatDB.txt")
df_PROGENy_LR<-read.table(file_df_PROGENy_LR,header=T,sep='\t')
PROGENy_pathway<-unique(c(unique(gsea_res_PROGENy_CD_R$Pathway),
                          unique(gsea_res_PROGENy_CD_TI$Pathway),
                          unique(gsea_res_PROGENy_UC_R$Pathway)))
pathway_R<-unique(df_PROGENy_LR[df_PROGENy_LR$Receptor!='','Pathway'])
pathway_miss_R<-setdiff(PROGENy_pathway,pathway_R)
pathway_L<-unique(df_PROGENy_LR[df_PROGENy_LR$Ligand!='','Pathway'])
pathway_miss_L<-setdiff(PROGENy_pathway,pathway_L)
pathway_miss_L_R<-intersect(pathway_miss_R,pathway_miss_L)
gsea_res_PROGENy_CD_R[gsea_res_PROGENy_CD_R$Pathway%in%pathway_miss_R,'Receptor_pct']<-1
gsea_res_PROGENy_CD_R[gsea_res_PROGENy_CD_R$Pathway%in%pathway_miss_L,'Ligand_pct']<-1
gsea_res_PROGENy_CD_R[gsea_res_PROGENy_CD_R$Pathway%in%pathway_miss_L_R,'Receptor_pct']<-0.5
gsea_res_PROGENy_CD_R[gsea_res_PROGENy_CD_R$Pathway%in%pathway_miss_L_R,'Ligand_pct']<-0
gsea_res_PROGENy_CD_TI[gsea_res_PROGENy_CD_TI$Pathway%in%pathway_miss_R,'Receptor_pct']<-1
gsea_res_PROGENy_CD_TI[gsea_res_PROGENy_CD_TI$Pathway%in%pathway_miss_L,'Ligand_pct']<-1
gsea_res_PROGENy_CD_TI[gsea_res_PROGENy_CD_TI$Pathway%in%pathway_miss_L_R,'Receptor_pct']<-0.5
gsea_res_PROGENy_CD_TI[gsea_res_PROGENy_CD_TI$Pathway%in%pathway_miss_L_R,'Ligand_pct']<-0
gsea_res_PROGENy_UC_R[gsea_res_PROGENy_UC_R$Pathway%in%pathway_miss_R,'Receptor_pct']<-1
gsea_res_PROGENy_UC_R[gsea_res_PROGENy_UC_R$Pathway%in%pathway_miss_L,'Ligand_pct']<-1
gsea_res_PROGENy_UC_R[gsea_res_PROGENy_UC_R$Pathway%in%pathway_miss_L_R,'Receptor_pct']<-0.5
gsea_res_PROGENy_UC_R[gsea_res_PROGENy_UC_R$Pathway%in%pathway_miss_L_R,'Ligand_pct']<-0
#rename PROGENy
gsea_res_PROGENy_CD_R$Pathway<-paste0(gsea_res_PROGENy_CD_R$Pathway,'(P)')
gsea_res_PROGENy_CD_R$Pathway<-gsub('CD64(P)','CD64(C)',gsea_res_PROGENy_CD_R$Pathway)
gsea_res_PROGENy_CD_TI$Pathway<-paste0(gsea_res_PROGENy_CD_TI$Pathway,'(P)')
gsea_res_PROGENy_CD_TI$Pathway<-gsub('CD64(P)','CD64(C)',gsea_res_PROGENy_CD_TI$Pathway)
gsea_res_PROGENy_UC_R$Pathway<-paste0(gsea_res_PROGENy_UC_R$Pathway,'(P)')
gsea_res_PROGENy_UC_R$Pathway<-gsub('CD64(P)','CD64(C)',gsea_res_PROGENy_UC_R$Pathway)
#combine PROGENy and ImmDict co-viz
gsea_res_CD_R<-rbind(gsea_res_PROGENy_CD_R[,-1],gsea_res_ImmDict_CD_R[,-1])
gsea_res_CD_TI<-rbind(gsea_res_PROGENy_CD_TI[,-1],gsea_res_ImmDict_CD_TI[,-1])
gsea_res_UC_R<-rbind(gsea_res_PROGENy_UC_R[,-1],gsea_res_ImmDict_UC_R[,-1])
df_celltype_colors<-read.delim(file.path(dir_data,"color_palette","IBD_main_color_palette.txt"),header=T,row.names=1,sep='\t')
#contrast group CD R vs CD TI
df_contrast<-data.frame(Ident1=c('Macrophage','CD4_T_cell','Colonocyte'),Ident2=c('DC','CD4_T_cell','Enterocyte'),celltype_colors=c('Myeloid','CD4_T_cell','Epithelial'))
for (ix in 1:nrow(df_contrast)){
  curr_ident1<-df_contrast[ix,'Ident1']
  curr_ident2<-df_contrast[ix,'Ident2']
  df_viz_1<-gsea_res_CD_R[gsea_res_CD_R$Context==paste0(curr_ident1,'_CD_R'),c('Pathway','-log10Padj','Receptor_pct','Ligand_pct')]
  colnames(df_viz_1)<-c('Pathway','-log10Padj(CD_R)','Receptor_pct(CD_R)','Ligand_pct(CD_R)')
  df_viz_2<-gsea_res_CD_TI[gsea_res_CD_TI$Context==paste0(curr_ident2,'_CD_TI'),c('Pathway','-log10Padj','Receptor_pct','Ligand_pct')]
  colnames(df_viz_2)<-c('Pathway','-log10Padj(CD_TI)','Receptor_pct(CD_TI)','Ligand_pct(CD_TI)')
  df_viz<-merge(df_viz_1,df_viz_2,by='Pathway')
  df_viz$Ligand_pct<-0
  df_viz[df_viz$`Ligand_pct(CD_R)`==1&df_viz$`Ligand_pct(CD_TI)`==1,'Ligand_pct']<-1
  df_viz$Receptor_pct<-0
  df_viz[df_viz$`Receptor_pct(CD_R)`==1&df_viz$`Receptor_pct(CD_TI)`==1,'Receptor_pct']<-1
  df_viz[df_viz$`Receptor_pct(CD_R)`==0.5&df_viz$`Receptor_pct(CD_TI)`==0.5,'Receptor_pct']<-0.5
  df_viz$Receptor_pct <- factor(df_viz$Receptor_pct)
  fit <- lm(`-log10Padj(CD_TI)` ~ `-log10Padj(CD_R)`, data = df_viz)
  r2 <- summary(fit)$r.squared
  cort <- cor.test(
    df_viz$`-log10Padj(CD_R)`,
    df_viz$`-log10Padj(CD_TI)`,
    method = "pearson"
  )
  rho <- unname(cort$estimate)
  pval <- cort$p.value    
  label_pval_exp <- floor(log10(pval))
  mant <- pval / 10^label_pval_exp
  label_stats <- sprintf(
    "atop(rho == %.2f, p == %.2f %%*%% 10^{%d})",
    rho, mant, label_pval_exp
  )
  curr_color_celltype<-df_contrast[ix,'celltype_colors']
  curr_color<-df_celltype_colors[curr_color_celltype,'Color']
  df_label<-df_viz[abs(df_viz$`-log10Padj(CD_R)`)>=1|abs(df_viz$`-log10Padj(CD_TI)`)>=1,]
  df_label$fontface <- ifelse(
    df_label$Ligand_pct == 1,
    "bold",
    "plain")
  if (curr_ident1!=curr_ident2){
    curr_title<-paste0(curr_ident1,' CD Rectum vs ', curr_ident2,' CD TI enrichment')
  }else if (curr_ident1==curr_ident2){
    curr_title<-paste0('Treg/Tfh CD Rectum vs CD TI enrichment')
  }
  ggplot(df_viz, aes(x = `-log10Padj(CD_R)`, y = `-log10Padj(CD_TI)`)) +
    geom_point(aes(fill = Receptor_pct), 
               size = 2.3, alpha = 1, color = curr_color, shape = 21, stroke = 1.2) +
    scale_fill_manual(
      values = c("0" = "white", "1" = curr_color, "0.5" = lighten(curr_color, 0.5))
    ) + 
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.8,
      color = "grey45",
      fill = "grey80"
    ) +
    annotate("text",
      x = -Inf,y = Inf,
      parse = TRUE,
      label = label_stats,
      hjust = -1.2,
      vjust = 1.1,
      size = 4) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +              # y-axis at 0
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +              # x-axis at 0
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey20") +     # threshold x = ±1
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey20") +     # threshold y = ±1
    geom_text_repel(data = df_label,
                    aes(label = Pathway,
                        fontface = fontface),
                    size = 3.6,
                    max.overlaps = Inf,
                    box.padding = 0.25,
                    point.padding = 0.25,
                    min.segment.length = 0.5)+
    labs(x = expression(CD~rectum ~ "["-log[10]("Padj")~"*"~"sign(NES)]"),
         y = expression(CD~TI~ "["-log[10]("Padj")~"*"~"sign(NES)]"),
         title=curr_title)+
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 14,hjust=0.5),
          plot.subtitle = element_text(size = 11),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          legend.position = "none")
  file_out<-file.path(dir_out,paste0(curr_ident1,'_CD_R_vs_',curr_ident2,
                                     '_CD_TI_ImmDict_PROGENy_CD64_cytokine_response_enrichment_scatter_plot.pdf'))
  ggsave(file_out,height=5.5,width=5.8)
}

#CD R vs UC R Colonocyte
df_viz_1<-gsea_res_CD_R[gsea_res_CD_R$Context=='Colonocyte_CD_R',c('Pathway','-log10Padj','Receptor_pct','Ligand_pct')]
colnames(df_viz_1)<-c('Pathway','-log10Padj(CD_R)','Receptor_pct(CD_R)','Ligand_pct(CD_R)')
df_viz_2<-gsea_res_UC_R[gsea_res_UC_R$Context=='Colonocyte_UC_R',c('Pathway','-log10Padj','Receptor_pct','Ligand_pct')]
colnames(df_viz_2)<-c('Pathway','-log10Padj(UC_R)','Receptor_pct(UC_R)','Ligand_pct(UC_R)')
df_viz<-merge(df_viz_1,df_viz_2,by='Pathway')
df_viz$Ligand_pct<-0
df_viz[df_viz$`Ligand_pct(CD_R)`==1&df_viz$`Ligand_pct(UC_R)`==1,'Ligand_pct']<-1
df_viz$Receptor_pct<-0
df_viz[df_viz$`Receptor_pct(CD_R)`==1&df_viz$`Receptor_pct(UC_R)`==1,'Receptor_pct']<-1
df_viz[df_viz$`Receptor_pct(CD_R)`==0.5&df_viz$`Receptor_pct(UC_R)`==0.5,'Receptor_pct']<-0.5
df_viz$Receptor_pct <- factor(df_viz$Receptor_pct)
cort <- cor.test(
  df_viz$`-log10Padj(CD_R)`,
  df_viz$`-log10Padj(UC_R)`,
  method = "pearson"
)
rho <- unname(cort$estimate)
pval <- cort$p.value    
label_pval_exp <- floor(log10(pval))
mant <- pval / 10^label_pval_exp
label_stats <- sprintf(
  "atop(rho == %.2f, p == %.2f %%*%% 10^{%d})",
  rho, mant, label_pval_exp
)
curr_color<-df_celltype_colors['Colonocyte','Color']
df_label<-df_viz[abs(df_viz$`-log10Padj(CD_R)`)>=1|abs(df_viz$`-log10Padj(UC_R)`)>=1,]
df_label$fontface <- ifelse(
  df_label$Ligand_pct == 1,
  "bold",
  "plain")
ggplot(df_viz, aes(x = `-log10Padj(CD_R)`, y = `-log10Padj(UC_R)`)) +
  geom_point(aes(fill = Receptor_pct), 
             size = 2.3, alpha = 1, color = curr_color, shape = 21, stroke = 1.2) +
  scale_fill_manual(
    values = c("0" = "white", "1" = curr_color, "0.5" = lighten(curr_color, 0.5))
  ) + 
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    color = "grey45",
    fill = "grey80"
  ) +
  annotate("text",
           x = -Inf,y = Inf,
           parse = TRUE,
           label = label_stats,
           hjust = -1.2,
           vjust = 1.1,
           size = 4) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +              # y-axis at 0
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +              # x-axis at 0
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey20") +     # threshold x = ±1
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey20") +     # threshold y = ±1
  geom_text_repel(data = df_label,
                  aes(label = Pathway,
                      fontface = fontface),
                  size = 3.6,
                  max.overlaps = Inf,
                  box.padding = 0.25,
                  point.padding = 0.25,
                  min.segment.length = 0.5)+
  labs(x = expression(CD~rectum ~ "["-log[10]("Padj")~"*"~"sign(NES)]"),
       y = expression(UC~rectum~ "["-log[10]("Padj")~"*"~"sign(NES)]"),
       title='Colonocyte CD Rectum vs UC Rectum enrichment')+
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 14,hjust=0.5),
        plot.subtitle = element_text(size = 11),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        legend.position = "none")
file_out<-file.path(dir_out,paste0('Colonocyte_CD_R_vs_UC_R_ImmDict_PROGENy_CD64_cytokine_response_enrichment_scatter_plot.pdf'))
ggsave(file_out,height=5.5,width=5.8)

#CD R vs UC R Macrophage
df_viz_1<-gsea_res_CD_R[gsea_res_CD_R$Context=='Macrophage_CD_R',c('Pathway','-log10Padj','Receptor_pct','Ligand_pct')]
colnames(df_viz_1)<-c('Pathway','-log10Padj(CD_R)','Receptor_pct(CD_R)','Ligand_pct(CD_R)')
df_viz_2<-gsea_res_UC_R[gsea_res_UC_R$Context=='Macrophage_UC_R',c('Pathway','-log10Padj','Receptor_pct','Ligand_pct')]
colnames(df_viz_2)<-c('Pathway','-log10Padj(UC_R)','Receptor_pct(UC_R)','Ligand_pct(UC_R)')
df_viz<-merge(df_viz_1,df_viz_2,by='Pathway')
df_viz$Ligand_pct<-0
df_viz[df_viz$`Ligand_pct(CD_R)`==1&df_viz$`Ligand_pct(UC_R)`==1,'Ligand_pct']<-1
df_viz$Receptor_pct<-0
df_viz[df_viz$`Receptor_pct(CD_R)`==1&df_viz$`Receptor_pct(UC_R)`==1,'Receptor_pct']<-1
df_viz[df_viz$`Receptor_pct(CD_R)`==0.5&df_viz$`Receptor_pct(UC_R)`==0.5,'Receptor_pct']<-0.5
df_viz$Receptor_pct <- factor(df_viz$Receptor_pct)
cort <- cor.test(
  df_viz$`-log10Padj(CD_R)`,
  df_viz$`-log10Padj(UC_R)`,
  method = "pearson"
)
rho <- unname(cort$estimate)
pval <- cort$p.value    
label_pval_exp <- floor(log10(pval))
mant <- pval / 10^label_pval_exp
label_stats <- sprintf(
  "atop(rho == %.2f, p == %.2f %%*%% 10^{%d})",
  rho, mant, label_pval_exp
)
curr_color<-df_celltype_colors['Macrophage','Color']
df_label<-df_viz[abs(df_viz$`-log10Padj(CD_R)`)>=1|abs(df_viz$`-log10Padj(UC_R)`)>=1,]
df_label$fontface <- ifelse(
  df_label$Ligand_pct == 1,
  "bold",
  "plain")
ggplot(df_viz, aes(x = `-log10Padj(CD_R)`, y = `-log10Padj(UC_R)`)) +
  geom_point(aes(fill = Receptor_pct), 
             size = 2.3, alpha = 1, color = curr_color, shape = 21, stroke = 1.2) +
  scale_fill_manual(
    values = c("0" = "white", "1" = curr_color, "0.5" = lighten(curr_color, 0.5))
  ) + 
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    color = "grey45",
    fill = "grey80"
  ) +
  annotate("text",
           x = -Inf,y = Inf,
           parse = TRUE,
           label = label_stats,
           hjust = -1.2,
           vjust = 1.1,
           size = 4) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +              # y-axis at 0
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +              # x-axis at 0
  geom_vline(xintercept =  1, linetype = "dashed", color = "grey20") +     # threshold x = ±1
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey20") +     # threshold y = ±1
  geom_text_repel(data = df_label,
                  aes(label = Pathway,
                      fontface = fontface),
                  size = 3.6,
                  max.overlaps = Inf,
                  box.padding = 0.25,
                  point.padding = 0.25,
                  min.segment.length = 0.5)+
  labs(x = expression(CD~rectum ~ "["-log[10]("Padj")~"*"~"sign(NES)]"),
       y = expression(UC~rectum~ "["-log[10]("Padj")~"*"~"sign(NES)]"),
       title='Macrophage CD Rectum vs UC Rectum enrichment')+
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 14,hjust=0.5),
        plot.subtitle = element_text(size = 11),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        legend.position = "none")
file_out<-file.path(dir_out,paste0('Macrophage_CD_R_vs_UC_R_ImmDict_PROGENy_CD64_cytokine_response_enrichment_scatter_plot.pdf'))
ggsave(file_out,height=5.5,width=5.8)