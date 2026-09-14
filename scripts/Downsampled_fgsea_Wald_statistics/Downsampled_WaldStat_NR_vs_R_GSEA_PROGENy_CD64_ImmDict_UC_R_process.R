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
# input data object
project_dir <- "."
script_dir <- file.path(project_dir, "scripts")
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results","GSEA_RNA_Status")
file_gsea_ImmDict<-file.path(results_dir,paste('Downsample/1Factor_CellType_Downsample_statWald/UC_R/final_ImmDict_max200/gsea_IBD_multiome_final_ImmDict_max200_I_vs_NI_Disease_Tissue_All.rds'))
file_sig_sets_ImmDict<-file.path(results_dir,paste('Downsample/1Factor_CellType_Downsample_statWald/UC_R/final_ImmDict_max200/sig_setsIBD_multiome_final_ImmDict_max200_All.txt'))
file_gsea_PROGENy_CD64<-file.path(results_dir,paste('Downsample/1Factor_CellType_Downsample_statWald/UC_R/final_PROGENy_CD64_max200/gsea_IBD_multiome_final_PROGENy_CD64_max200_I_vs_NI_Disease_Tissue_All.rds'))
file_sig_sets_PROGENy_CD64<-file.path(results_dir,paste('Downsample/1Factor_CellType_Downsample_statWald/UC_R/final_PROGENy_CD64_max200/sig_setsIBD_multiome_final_PROGENy_CD64_max200_All.txt'))
file_df_frac<-file.path(data_dir,paste('nominal_expression/df_5pct_min_exp_by_CellType_Disease_Tissue_Status.rds'))
data_assay <- 'RNA' # data assay
order_disease_tissue<-c('UC_R')
order_celltype<-c('CD4_T_cell','CD4_T_Eff','CD8_T_cell','Naive_T_cell','gdT_cell','NK_cell',
'B_cell','Plasma_cell','Macrophage','DC','Fibroblast','BEST4_Colonocyte','Colonocyte','Immature_Colonocyte',
'Enterocyte','Immature_Enterocyte','Enteroendocrine','Endothelial','Goblet','Immature_Goblet','TA','Cycling_TA','Stem','Tuft','Paneth','M_cell')

#####
# Params
file_save <- 'IBD_multiome' # file save base
# output directories
dir_out <- file.path(project_dir, "results","GSEA_RNA_Status", "Downsample_statWald", "PROGENy_CD64_ImmDict", "UC_R")
dir.create(dir_out, showWarnings=FALSE)
order_clustersubtype<-NULL
for (ix in order_celltype){
    for (jx in order_disease_tissue){
        order_clustersubtype<-c(order_clustersubtype,paste(ix,jx,sep='_'))
    }
}
df_frac<-readRDS(file_df_frac)
df_frac<-df_frac[df_frac$DiseaseTissue=='UC_R',]
df_frac<-df_frac[df_frac$min5pct==1,]
print('load enrichment analysis')
sig_sets_ImmDict <- readLines(file_sig_sets_ImmDict)
gsea_ImmDict<-readRDS(file_gsea_ImmDict)
sig_sets_PROGENy_CD64<-readLines(file_sig_sets_PROGENy_CD64)
gsea_PROGENy_CD64<-readRDS(file_gsea_PROGENy_CD64)
order_clustersubtype<-intersect(order_clustersubtype,names(gsea_ImmDict))
# pval enrichment heatmaps - pathway vs. time point
#all pathway not just sig pathway 
df_ImmDict<-NULL
for (ix in order_clustersubtype){
    curr_gsea<-gsea_ImmDict[[ix]]
    if (nrow(curr_gsea)>0){
      curr_df <- data.frame(pathway=curr_gsea[,'pathway'],
                          Padj=-log10(curr_gsea[,'padj'])*sign(curr_gsea[,'NES']), 
                          NES=curr_gsea[,'NES'], Context=ix)
      df_ImmDict<-rbind(df_ImmDict,curr_df)
    }
}
colnames(df_ImmDict)<-c('Cytokine','-log10Padj','NES','Context')
df_ImmDict$`-log10Padj`[which(is.na(df_ImmDict$`-log10Padj`))] <- 0
file_ImmDict_LR <- file.path(
  data_dir, "gene_sets", "ImmuneDictionary", "ImmDict_LR_pairs_by_cellchatDB.txt"
)
ImmDict_LR<-read.table(file_ImmDict_LR,header=T,sep='\t')
df_ImmDict<-merge(df_ImmDict,ImmDict_LR,by='Cytokine')
df_ImmDict$Context<-factor(df_ImmDict$Context,levels=order_clustersubtype)
df_ImmDict <- df_ImmDict %>%
  mutate(Status = ifelse(NES > 0, "I", ifelse(NES < 0, "NI", NA)))
df_frac$Context<-paste(df_frac$CellType,df_frac$DiseaseTissue,sep='_')
df_frac_ligand_disease_tissue_status<-df_frac %>%
  group_by(CellType,DiseaseTissue, Status) %>%
  summarise(Genes = list(unique(Gene)),.groups = "drop")
df_frac_ligand_disease_tissue_status<-as.data.frame(df_frac_ligand_disease_tissue_status)
for (ix in 1:nrow(df_ImmDict)){
  curr_context<-as.character(df_ImmDict[ix,'Context'])
  curr_status<-df_ImmDict[ix,'Status']
  curr_receptor<-unlist(strsplit(df_ImmDict[ix,'Receptor'],'_'))
  curr_min_5pct_gene_celltype_context<-df_frac[df_frac$Context%in%curr_context&df_frac$Status%in%curr_status,'Gene']
  df_ImmDict[ix,'Receptor_5pct_min']<-ifelse(all(curr_receptor %in% curr_min_5pct_gene_celltype_context), 1, 0)
  curr_ligand<-unlist(strsplit(df_ImmDict[ix,'Ligand'],'_'))
  curr_min_5pct_gene_context<-df_frac_ligand_disease_tissue_status[df_frac_ligand_disease_tissue_status$Status%in%curr_status,]
  df_ImmDict[ix,'Ligand_5pct_min']<-ifelse(any(sapply(curr_min_5pct_gene_context$Genes,function(g) all(curr_ligand %in% g))), 1, 0)
}
df_ImmDict<-as.data.frame(df_ImmDict)
file_out<-file.path(dir_out,paste0('ImmDict_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R.txt'))
write.table(df_ImmDict,file_out,sep='\t',col.names=T,row.names=F,quote=F)
df_pathway_ImmDict<-df_ImmDict %>%
  group_by(Cytokine, Pathway, Context,`-log10Padj`, NES, Status) %>%
  summarise(Receptor_pct = max(Receptor_5pct_min, na.rm = TRUE), Ligand_pct=max(Ligand_5pct_min,na.rm=TRUE),.groups = 'drop')
df_pathway_ImmDict<-as.data.frame(df_pathway_ImmDict)
file_out<-file.path(dir_out,paste0('ImmDict_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
write.table(df_pathway_ImmDict,file_out,sep='\t',col.names=T,row.names=F,quote=F)

#PROGENy results
df_PROGENy_CD64<-NULL
for (ix in order_clustersubtype){
    curr_gsea<-gsea_PROGENy_CD64[[ix]]
    if (nrow(curr_gsea)>0){
      curr_df <- data.frame(pathway=curr_gsea[,'pathway'],
                          Padj=-log10(curr_gsea[,'padj'])*sign(curr_gsea[,'NES']), 
                          NES=curr_gsea[,'NES'], Context=ix)
      df_PROGENy_CD64<-rbind(df_PROGENy_CD64,curr_df)
    }
}
colnames(df_PROGENy_CD64)<-c('Cytokine','-log10Padj','NES','Context')
df_PROGENy_CD64$`-log10Padj`[which(is.na(df_PROGENy_CD64$`-log10Padj`))] <- 0
df_PROGENy_CD64$Cytokine<-gsub('PROGENy_','',df_PROGENy_CD64$Cytokine)
df_PROGENy_CD64$Cytokine<-gsub('Gao2020_CD64_wtPadj_max200','CD64',df_PROGENy_CD64$Cytokine)
file_PROGENy_CD64_LR<-file.path(data_dir, "gene_sets", "PROGENy_CD64", "PROGENy_CD64_LR_pairs_by_cellchatDB.txt")
PROGENy_CD64_LR<-read.table(file_PROGENy_CD64_LR,header=T,sep='\t')
df_PROGENy_CD64<-merge(df_PROGENy_CD64,PROGENy_CD64_LR,by='Cytokine',all.x=T)
df_PROGENy_CD64$Context<-factor(df_PROGENy_CD64$Context,levels=order_clustersubtype)
df_PROGENy_CD64 <- df_PROGENy_CD64 %>%
  mutate(Status = ifelse(NES > 0, "I", ifelse(NES < 0, "NI", NA)))
for (ix in 1:nrow(df_PROGENy_CD64)){
  curr_context<-as.character(df_PROGENy_CD64[ix,'Context'])
  curr_status<-df_PROGENy_CD64[ix,'Status']
  curr_receptor<-unlist(strsplit(df_PROGENy_CD64[ix,'Receptor'],'_'))
  if (length(curr_receptor)>0){
    curr_min_5pct_gene_celltype_context<-df_frac[df_frac$Context%in%curr_context&df_frac$Status%in%curr_status,'Gene']
    df_PROGENy_CD64[ix,'Receptor_5pct_min']<-ifelse(all(curr_receptor %in% curr_min_5pct_gene_celltype_context), 1, 0)
  } else{
    df_PROGENy_CD64[ix,'Receptor_5pct_min']<-0
  }
  curr_ligand<-unlist(strsplit(df_PROGENy_CD64[ix,'Ligand'],'_'))
  if(length(curr_ligand)>0){
    curr_min_5pct_gene_context<-df_frac_ligand_disease_tissue_status[df_frac_ligand_disease_tissue_status$Status%in%curr_status,]
    df_PROGENy_CD64[ix,'Ligand_5pct_min']<-ifelse(any(sapply(curr_min_5pct_gene_context$Genes,function(g) all(curr_ligand %in% g))), 1, 0)
  }else{
    df_PROGENy_CD64[ix,'Ligand_5pct_min']<-0
  }
}
file_out<-file.path(dir_out,paste0('PROGENy_CD64_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R.txt'))
write.table(df_PROGENy_CD64,file_out,sep='\t',col.names=T,row.names=F,quote=F)
df_pathway_PROGENy_CD64<-df_PROGENy_CD64 %>%
  group_by(Cytokine,Pathway, Context,`-log10Padj`, NES, Status) %>%
  summarise(Receptor_pct = max(Receptor_5pct_min, na.rm = TRUE), Ligand_pct=max(Ligand_5pct_min,na.rm=TRUE),.groups = 'drop')
df_pathway_PROGENy_CD64<-as.data.frame(df_pathway_PROGENy_CD64)
file_out<-file.path(dir_out,paste0('PROGENy_CD64_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
write.table(df_pathway_PROGENy_CD64,file_out,sep='\t',col.names=T,row.names=F,quote=F)