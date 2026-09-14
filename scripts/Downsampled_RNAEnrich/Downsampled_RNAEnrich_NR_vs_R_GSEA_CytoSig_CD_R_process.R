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
file_gsea<-file.path(results_dir,paste('Downsample/1Factor_CellType_Downsample_RNAEnrich/CD_R/CytoSig/gsea_IBD_multiome_CytoSig_I_vs_NI_Disease_Tissue_All.rds'))
file_sig_sets<-file_gsea<-file.path(results_dir,paste('Downsample/1Factor_CellType_Downsample_RNAEnrich/CD_R/CytoSig/sig_setsIBD_multiome_CytoSig_All.txt'))
file_df_frac<-file.path(data_dir,paste('nominal_expression/df_5pct_min_exp_by_CellType_Disease_Tissue_Status.rds'))
data_assay <- 'RNA' # data assay
order_disease_tissue<-c('CD_R')
order_celltype<-c('CD4_T_cell','CD4_T_Eff','CD8_T_cell','Naive_T_cell','gdT_cell','NK_cell',
'B_cell','Plasma_cell','Macrophage','DC','Fibroblast','BEST4_Colonocyte','Colonocyte','Immature_Colonocyte',
'Enterocyte','Immature_Enterocyte','Enteroendocrine','Endothelial','Goblet','Immature_Goblet','TA','Cycling_TA','Stem','Tuft','Paneth','M_cell')

#####
# Params
file_save <- 'IBD_multiome' # file save base
# output directories
dir_out <- file.path(project_dir, "results","GSEA_RNA_Status", "Downsample_RNAEnrich", "CytoSig", "CD_R")
dir.create(dir_out, showWarnings=FALSE)
order_clustersubtype<-NULL
for (ix in order_celltype){
    for (jx in order_disease_tissue){
        order_clustersubtype<-c(order_clustersubtype,paste(ix,jx,sep='_'))
    }
}
df_frac<-readRDS(file_df_frac)
df_frac<-df_frac[df_frac$DiseaseTissue=='CD_R',]
df_frac<-df_frac[df_frac$min5pct==1,]
print('load enrichment analysis')
sig_sets <- readLines(file_sig_sets)
gsea<-readRDS(file_gsea)
order_clustersubtype<-intersect(order_clustersubtype,names(gsea))
# pval enrichment heatmaps - pathway vs. time point
#all pathway not just sig pathway 
df<-NULL
for (ix in order_clustersubtype){
    curr_gsea<-gsea[[ix]]
    if (nrow(curr_gsea)>0){
      curr_df <- data.frame(pathway=curr_gsea[,'pathway'],NGene=curr_gsea[,'n.genes'],
                          Padj=-log10(curr_gsea[,'padj'])*sign(curr_gsea[,'coeff']), 
                          coeff=curr_gsea[,'coeff'], odds.ratio=curr_gsea[,'odds.ratio'],Context=ix)
      df<-rbind(df,curr_df)
    }
}
colnames(df)<-c('Pathway','NGene','-log10Padj','coeff','odds.ratio','Context')
df$`-log10Padj`[which(is.na(df$`-log10Padj`))] <- 0
file_CytoSig_LR <- file.path(
  data_dir, "gene_sets", "CytoSig", "CytoSig_LR_pairs_by_cellchatDB.txt"
)
CytoSig_LR<-read.table(file_CytoSig_LR,header=T,sep='\t')
df<-merge(df,CytoSig_LR,by='Pathway')
df$Context<-factor(df$Context,levels=order_clustersubtype)
df <- df %>%
  mutate(Status = ifelse(coeff> 0, "I", ifelse(coeff < 0, "NI", NA)))
df_frac$Context<-paste(df_frac$CellType,df_frac$DiseaseTissue,sep='_')
df_frac_ligand_disease_tissue_status<-df_frac %>%
  group_by(CellType,DiseaseTissue, Status) %>%
  summarise(Genes = list(unique(Gene)),.groups = "drop")
df_frac_ligand_disease_tissue_status<-as.data.frame(df_frac_ligand_disease_tissue_status)
for (ix in 1:nrow(df)){
  curr_context<-as.character(df[ix,'Context'])
  curr_status<-df[ix,'Status']
  curr_receptor<-unlist(strsplit(df[ix,'Receptor'],'_'))
  curr_min_5pct_gene_celltype_context<-df_frac[df_frac$Context%in%curr_context&df_frac$Status%in%curr_status,'Gene']
  df[ix,'Receptor_5pct_min']<-ifelse(all(curr_receptor %in% curr_min_5pct_gene_celltype_context), 1, 0)
  curr_ligand<-unlist(strsplit(df[ix,'Ligand'],'_'))
  curr_min_5pct_gene_context<-df_frac_ligand_disease_tissue_status[df_frac_ligand_disease_tissue_status$Status%in%curr_status,]
  df[ix,'Ligand_5pct_min']<-ifelse(any(sapply(curr_min_5pct_gene_context$Genes,function(g) all(curr_ligand %in% g))), 1, 0)
}
df<-as.data.frame(df)
file_out<-file.path(dir_out,paste0('CytoSig_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R.txt'))
write.table(df,file_out,sep='\t',col.names=T,row.names=F,quote=F)
df_pathway<-df %>%
  group_by(Pathway, Context,`-log10Padj`, NGene,coeff, odds.ratio,Status) %>%
  summarise(Receptor_pct = max(Receptor_5pct_min, na.rm = TRUE), Ligand_pct=max(Ligand_5pct_min,na.rm=TRUE),.groups = 'drop')
df_pathway<-as.data.frame(df_pathway)
file_out<-file.path(dir_out,paste0('CytoSig_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
write.table(df_pathway,file_out,sep='\t',col.names=T,row.names=F,quote=F)