rm(list=ls())
options(stringsAsFactors=FALSE)
set.seed(42)
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
  library(dplyr)
  library(dendextend)
  library(plotrix)
  library(RColorBrewer)
  library(tidyverse)
  library(scales)
  library(ggnewscale)
  library(ggtext)
})
project_dir <- "."
script_dir <- file.path(project_dir, "scripts")
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results","GSEA_RNA_Status")
# output directory
dir_out <- file.path(results_dir,"Downsampled_fgsea_RNAEnrich_co_viz","ImmDict_CytoSig_shared_cytokine_coviz")
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)
#params
order_disease_tissue<-c('CD_TI')
order_celltype<-c('CD4_T_cell','CD4_T_Eff','CD8_T_cell','Naive_T_cell','gdT_cell',
'B_cell','Plasma_cell','Macrophage','DC','Fibroblast','Endothelial',
'Enterocyte','Immature_Enterocyte','Colonocyte','Immature_Colonocyte',
'TA','Cycling_TA','Goblet','Immature_Goblet','Stem','Enteroendocrine')

#input results
file_df_results_ImmDict<-file.path(results_dir,paste('Downsample_statWald/PROGENy_CD64_ImmDict/CD_TI/ImmDict_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_df_results_ImmDict_RNAEnrich<-file.path(results_dir,paste('Downsample_RNAEnrich/PROGENy_CD64_ImmDict/CD_TI/ImmDict_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_df_results_CytoSig<-file.path(results_dir,paste('Downsample_statWald/CytoSig/CD_TI/CytoSig_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_df_results_CytoSig_RNAEnrich<-file.path(results_dir,paste('Downsample_RNAEnrich/CytoSig/CD_TI/CytoSig_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_cytokine_IBD_enhancer_targets_L<-file.path(data_dir,paste('gene_sets/ImmuneDictionary/ImmDict_pathway_Ligand_as_IBD_enhancer_targets.txt'))
file_cytokine_IBD_enhancer_targets_R<-file.path(data_dir,paste('gene_sets/ImmuneDictionary/ImmDict_pathway_Receptor_as_IBD_enhancer_targets.txt'))
file_df_order_cytokine_ImmDict<-file.path(results_dir,paste('Downsample_statWald/PROGENy_CD64_ImmDict/CD_R/sig_sets_ImmDict_cytokine_family_order_by_max_log10Padj.txt'))
file_keep_context<-file.path(data_dir,'keep_context/keep_context_CD_TI_Nsample_3_min_filter.txt')
file_df_cytokine_family<-file.path(data_dir,'gene_sets/CytoSig/CytoSig_cytokine_family_class.txt')

#results load
df_results_ImmDict<-read.table(file_df_results_ImmDict,sep='\t',header=T,check.names=F)
df_results_ImmDict_RNAEnrich<-read.table(file_df_results_ImmDict_RNAEnrich,sep='\t',header=T,check.names=F)
df_results_CytoSig<-read.table(file_df_results_CytoSig,sep='\t',header=T,check.names=F)
df_results_CytoSig_RNAEnrich<-read.table(file_df_results_CytoSig_RNAEnrich,sep='\t',header=T,check.names=F)
order_clustersubtype<-NULL
for (ix in order_celltype){
  for (jx in order_disease_tissue){
    order_clustersubtype<-c(order_clustersubtype,paste(ix,jx,sep='_'))
  }
}
keep_context<-readLines(file_keep_context)
keep_context<-intersect(order_clustersubtype,keep_context)
missing_context<-setdiff(order_clustersubtype,keep_context)
#filter ImmDict  results to inlcude only Nsample>=3 both I and NI results
df_results_ImmDict<-df_results_ImmDict[df_results_ImmDict$Context%in%keep_context,]
df_results_ImmDict_RNAEnrich<-df_results_ImmDict_RNAEnrich[df_results_ImmDict_RNAEnrich$Context%in%keep_context,]
df_results_CytoSig<-df_results_CytoSig[df_results_CytoSig$Context%in%keep_context,]
df_results_CytoSig_RNAEnrich<-df_results_CytoSig_RNAEnrich[df_results_CytoSig_RNAEnrich$Context%in%keep_context,]
overlap_pathway<-intersect(unique(df_results_ImmDict$Pathway),unique(df_results_CytoSig$Pathway))
df_results_ImmDict<-df_results_ImmDict[,-1]
df_missing <- data.frame(
    Pathway = rep(overlap_pathway, times = length(missing_context)),
    Context = rep(missing_context, each = length(overlap_pathway)))
df_missing$`-log10Padj`<-0
df_missing$NES<-0
df_missing$Status<-'NI'
df_missing$Receptor_pct<-0
df_missing$Ligand_pct<-0
df_results_ImmDict<-rbind(df_results_ImmDict,df_missing)
df_results_ImmDict$Database<-'ImmDict'
df_results_ImmDict$ContextDatabase<-paste(df_results_ImmDict$Context,df_results_ImmDict$Database,sep='_')
df_results_CytoSig<-rbind(df_results_CytoSig,df_missing)
df_results_CytoSig$Database<-'CytoSig'
df_results_CytoSig$ContextDatabase<-paste(df_results_CytoSig$Context,df_results_CytoSig$Database,sep='_')
df_results<-rbind(df_results_ImmDict,df_results_CytoSig)
df_results<-df_results[df_results$Pathway%in%overlap_pathway,]
#RNAEnrich results add the missing context 
df_missing_RNAEnrich <- data.frame(
    Pathway = rep(overlap_pathway, times = length(missing_context)),
    Context = rep(missing_context, each = length(overlap_pathway)))
df_missing_RNAEnrich$`-log10Padj`<-0
df_missing_RNAEnrich$NGene<-0
df_missing_RNAEnrich$coeff<-0
df_missing_RNAEnrich$odds.ratio<-0
df_missing_RNAEnrich$Status<-'NI'
df_missing_RNAEnrich$Receptor_pct<-0
df_missing_RNAEnrich$Ligand_pct<-0
df_results_ImmDict_RNAEnrich<-df_results_ImmDict_RNAEnrich[,-1]

df_results_ImmDict_RNAEnrich<-rbind(df_results_ImmDict_RNAEnrich,df_missing_RNAEnrich)
df_results_ImmDict_RNAEnrich$Database<-'ImmDict'
df_results_ImmDict_RNAEnrich$ContextDatabase<-paste(df_results_ImmDict_RNAEnrich$Context,df_results_ImmDict_RNAEnrich$Database,sep='_')
df_results_CytoSig_RNAEnrich<-rbind(df_results_CytoSig_RNAEnrich,df_missing_RNAEnrich)
df_results_CytoSig_RNAEnrich$Database<-'CytoSig'
df_results_CytoSig_RNAEnrich$ContextDatabase<-paste(df_results_CytoSig_RNAEnrich$Context,df_results_CytoSig_RNAEnrich$Database,sep='_')
df_results_RNAEnrich<-rbind(df_results_ImmDict_RNAEnrich,df_results_CytoSig_RNAEnrich)
#order column by cell type, database
order_clustersubtype_database <- as.vector(rbind(
  paste0(order_clustersubtype, "_ImmDict"),
  paste0(order_clustersubtype, "_CytoSig")))
# order cytokines
df_order_cytokine_ImmDict<-read.table(file_df_order_cytokine_ImmDict,sep='\t',header=T)
df_order_cytokine_ImmDict<-df_order_cytokine_ImmDict[,c('Pathway','Family')]
df_order_cytokine_ImmDict<-df_order_cytokine_ImmDict[df_order_cytokine_ImmDict$Pathway%in%overlap_pathway,]
extra_pathway<-setdiff(overlap_pathway,df_order_cytokine_ImmDict$Pathway)
# add missing cytokine to the corresponding family group
df_cytokine_family<-read.table(file_df_cytokine_family,sep='\t',header=T)
df_cytokine_family<-df_cytokine_family[df_cytokine_family$Pathway%in%extra_pathway,]
for (fam in unique(df_cytokine_family$Family)) {

  add_rows <- df_cytokine_family %>% filter(Family == fam)

  last_idx <- max(which(df_order_cytokine_ImmDict$Family == fam))

  df_order_cytokine_ImmDict <- bind_rows(
    df_order_cytokine_ImmDict[1:last_idx, ],
    add_rows,
    df_order_cytokine_ImmDict[(last_idx + 1):nrow(df_order_cytokine_ImmDict), ]
  )
} 
order_cytokine<-df_order_cytokine_ImmDict$Pathway
#IBD enhancer targets
IBD_enhancer_targets_L_pathway<-readLines(file_cytokine_IBD_enhancer_targets_L)
IBD_enhancer_targets_R_pathway<-readLines(file_cytokine_IBD_enhancer_targets_R)
IBD_enhancer_targets_pathway<-unique(c(IBD_enhancer_targets_L_pathway,IBD_enhancer_targets_R_pathway))
IBD_enhancer_targets_pathway<-intersect(IBD_enhancer_targets_pathway,overlap_pathway)
# pathway ligand detected
pathway_ligand_detected<-unique(df_results[df_results$Ligand_pct==1,'Pathway'])
# visualize the -log10Padj 
df_viz_padj <- acast(df_results, Pathway ~ ContextDatabase, value.var = "-log10Padj")
df_viz_padj<-df_viz_padj[order_cytokine,order_clustersubtype_database]
df_viz_padj[is.na(df_viz_padj)]<-0
#set NA for missing context results 
missing_context_database<-as.vector(rbind(
  paste0(missing_context, "_ImmDict"),
  paste0(missing_context, "_CytoSig")))
df_viz_padj[,missing_context_database]<-NA
df_viz_receptor <- acast(df_results, Pathway ~ ContextDatabase, value.var = "Receptor_pct")
df_viz_receptor<-df_viz_receptor[order_cytokine,order_clustersubtype_database]
df_viz_receptor[is.na(df_viz_receptor)]<-0
df_viz_receptor[,missing_context_database]<-NA
#RNA Enrich results 
df_viz_padj_RNAEnrich <- acast(df_results_RNAEnrich, Pathway ~ ContextDatabase, value.var = "-log10Padj")
df_viz_padj_RNAEnrich<-df_viz_padj_RNAEnrich[order_cytokine,order_clustersubtype_database]
df_viz_padj_RNAEnrich[is.na(df_viz_padj_RNAEnrich)]<-0
df_viz_padj_RNAEnrich[,missing_context_database]<-NA
#Visualize the top cell type annotation
meta<-data.frame(order_clustersubtype_database)
colnames(meta)<-'CellTypeDiseaseTissueDatabase'
meta$CellType<-sapply(strsplit(meta$CellTypeDiseaseTissueDatabase,split='_'),function(x) paste(x[1:(length(x)-3)],collapse='_'))
meta$Database<-sapply(strsplit(meta$CellTypeDiseaseTissueDatabase,split='_'),function(x) x[length(x)])
meta<-meta[,c('CellType','Database')]
celltype_colors<-read.delim('/data/GastroAI/genomics/analysis/multiome/liver/IBD_analysis/figures/V2/color_palette/IBD_main_color_palette.txt',header=T,row.names=1,sep='\t')
celltype_colors<-celltype_colors[order_celltype,,drop=F]
heat_meta_col <- list()
heat_meta_col[['CellType']] <- setNames(celltype_colors$Color, order_celltype)
heat_meta_col[['Database']] <- c('ImmDict'='black','CytoSig'='grey')
heat_max_pval<-20
cols <- rev(brewer.pal(11, "PRGn"))  
heat_col<-colorRamp2(c(-heat_max_pval,-15,-10,-5,-3, 0,3,5,10,15,heat_max_pval),cols)
lgd <- Legend(col_fun=heat_col, title='-log10(Padj)*signed(NES)',at=c(-heat_max_pval,-15,-10,-5,-3, 0,3,5,10,15,heat_max_pval),labels=c('-20','-15','-10','-5','-3','0','3','5','10','15','20'),
              legend_height = unit(4, "cm"),
              labels_gp = gpar(fontsize=9, fontface = 1), title_gp = gpar(fontsize=10, fontface = 2))
celltype_annot <- HeatmapAnnotation(df=meta, col=heat_meta_col,show_legend=FALSE,simple_anno_size = unit(3, "mm"),
                                    annotation_name_gp = gpar(fontsize=0,fontface=2))
lgd_Database= Legend(labels = c('ImmDict','CytoSig'), title = "Database", legend_gp = gpar(fill = heat_meta_col[['Database']] ),
                   legend_height = unit(6, "cm"),legend_width = unit(4, "cm"),
                   labels_gp = gpar(fontsize=9, fontface = 1), title_gp = gpar(fontsize=10, fontface = 2))# Heatmap
print('Heatmap')
col_section<- colnames(df_viz_padj)
col_section <- gsub('_CD_TI','',col_section)
col_section <- gsub('_ImmDict','',col_section)
col_section <- gsub('_CytoSig','',col_section)
col_split = data.frame(col_section)
colnames(col_split)<-'celltype'
col_split$celltype<-factor(col_split$celltype,levels=order_celltype)
row_split= data.frame(df_order_cytokine_ImmDict$Family)
colnames(row_split)<-'Family'
row_split$Family<-factor(row_split$Family,levels=unique(df_order_cytokine_ImmDict$Family))
#add information of row labels of ligand expression and IBD risk enhancer targets
row_labels <- rownames(df_viz_padj)
is_bold <- row_labels %in% pathway_ligand_detected
is_dagger  <- row_labels %in% IBD_enhancer_targets_pathway
row_names_gp <- gpar(
  fontsize = 8,
  fontface = ifelse(is_bold, "bold", "plain")
)
row_labels <- ifelse(
  is_dagger, 
  paste0(row_labels, "<sup>i</sup>"), #press option+T on keyboard
  row_labels
)
# Label celltypes
label_celltype <- colnames(df_viz_padj)
label_celltype <- gsub('_CD_TI','',label_celltype)
label_celltype <- gsub('_ImmDict','',label_celltype)
label_celltype <- gsub('_CytoSig','',label_celltype)
label_celltype<-gsub('CD4_T_cell','Treg/Tfh',label_celltype)
label_celltype<-gsub('CD4_T_Eff','CD4+ T Eff',label_celltype)
label_celltype<-gsub('CD8_T_cell','CD8+ T cell',label_celltype)
label_celltype<-gsub('gdT_cell','gd T cell',label_celltype)
label_celltype<-gsub('Immature_','Imm.',label_celltype)
label_celltype<-gsub('_',' ',label_celltype)
label_celltype[which(duplicated(label_celltype))] <- ''
file_out <- file.path(dir_out, paste0('heat_ImmDict_vs_CytoSig_top200_co_viz_fgsea_RNAEnrich_filter_Nsample_3_min_CD_TI.pdf'))
ht <- Heatmap(df_viz_padj, col=heat_col, cluster_columns=FALSE, cluster_rows=FALSE,
               show_row_dend=FALSE, show_row_names=TRUE, show_column_names=TRUE, 
               column_names_rot=90, column_names_side='top',column_labels=gt_render(label_celltype),
               column_names_gp=gpar(fontsize=10, fontface=1),cluster_column_slices = FALSE,
               row_labels = gt_render(row_labels),row_names_gp=row_names_gp,show_heatmap_legend=FALSE,
               column_split = col_split,row_split=row_split, row_title_rot = 0,
               row_names_side='left',row_title_gp=gpar(fontsize=10, fontface=2),
               row_gap = unit(0, "mm"),column_title=NULL, column_title_gp=gpar(fontsize=12, fontface=2),
               width = ncol(df_viz_padj)*unit(3, "mm"), height = nrow(df_viz_padj)*unit(3, "mm"),
               top_annotation=celltype_annot, 
               column_gap = unit(0, "mm"), border = TRUE,use_raster = FALSE,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 val_1<- df_viz_padj[i, j]
                 val_2<- df_viz_padj_RNAEnrich[i, j]
                 val_3<-df_viz_receptor[i, j]
                 if (!is.na(val_3) &&val_3 == 1&&
                     !is.na(val_1) && val_1 > 1&&
                     !is.na(val_2) && val_2 > 1) {
                  grid.rect(x = x, y = y,width = width,height = height,
                  gp = gpar(fill = NA,col = "black",lwd = 1.5))}
                 if (!is.na(val_3) &&val_3 == 1&&
                     !is.na(val_1) && val_1 <(-1)&&
                     !is.na(val_2) && val_2 <(-1)) {
                  grid.rect(x = x, y = y,width = width,height = height,
                  gp = gpar(fill = NA,col = "black",lwd = 1.5))}
                 if (!is.na(val_1) && val_1 > 1&&
                     !is.na(val_2) && val_2 > 1) {
                   grid.text("*", x = x, y = y-height * 0.25,
                             gp = gpar(fontsize = 9, fontface = "bold"))}
                 if (!is.na(val_1) && val_1 < (-1)&&
                     !is.na(val_2) && val_2 < (-1)) {
                   grid.text("*", x = x, y = y-height * 0.25,
                             gp = gpar(fontsize = 9, fontface = "bold"))}
                             })
calc_ht_size = function(ht, unit = "inch") {
  pdf(NULL)
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()
  
  c(w, h)
}
size = calc_ht_size(ht)
pdf(file_out, width=size[1]+3, height=size[2])
ht <- draw(ht,padding = unit(c(2, 15, 2, 2), "mm"),annotation_legend_list = list(lgd_Database,lgd))
dev.off()

