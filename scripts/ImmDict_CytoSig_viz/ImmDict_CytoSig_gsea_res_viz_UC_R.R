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
#work directory
project_dir <- "."
script_dir <- file.path(project_dir, "scripts")
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results","GSEA_RNA_Status")
#input files
file_gsea_ImmDict<-file.path(results_dir,paste('Downsample_statWald/PROGENy_CD64_ImmDict/UC_R/ImmDict_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_gsea_ImmDict_RNAEnrich<-file.path(results_dir,paste('Downsample_RNAEnrich/PROGENy_CD64_ImmDict/UC_R/ImmDict_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_df_order_cytokine_ImmDict<-file.path(results_dir,paste('Downsampled_fgsea_RNAEnrich_co_viz/UC_R/sig_sets_ImmDict_cytokine_family_order.txt'))

file_gsea_CytoSig<-file.path(results_dir,paste('Downsample_statWald/CytoSig/UC_R/CytoSig_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_gsea_CytoSig_RNAEnrich<-file.path(results_dir,paste('Downsample_RNAEnrich/CytoSig/UC_R/CytoSig_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_df_order_cytokine_CytoSig<-file.path(results_dir,paste('Downsampled_fgsea_RNAEnrich_co_viz/UC_R/sig_sets_CytoSig_cytokine_family_order.txt'))

file_cytokine_IBD_enhancer_targets_L_ImmDict<-file.path(data_dir,paste('gene_sets/ImmuneDictionary/ImmDict_pathway_Ligand_as_IBD_enhancer_targets.txt'))
file_cytokine_IBD_enhancer_targets_R_ImmDict<-file.path(data_dir,paste('gene_sets/ImmuneDictionary/ImmDict_pathway_Receptor_as_IBD_enhancer_targets.txt'))
file_cytokine_IBD_enhancer_targets_L_CytoSig<-file.path(data_dir,paste('gene_sets/CytoSig/CytoSig_pathway_Ligand_as_IBD_enhancer_targets.txt'))
file_cytokine_IBD_enhancer_targets_R_CytoSig<-file.path(data_dir,paste('gene_sets/CytoSig/CytoSig_pathway_Receptor_as_IBD_enhancer_targets.txt'))

file_keep_context<-file.path(data_dir,'keep_context/keep_context_UC_R_Nsample_3_min_filter.txt')
file_df_Nsample<-file.path(data_dir,'keep_context/df_sample_gene_expression_Nreads.txt')

#####
# Params
file_save <- 'IBD_multiome' # file save base
data_assay <- 'RNA' # data assay
order_disease_tissue<-c('UC_R')
order_celltype<-c('CD4_T_cell','CD4_T_Eff','CD8_T_cell','Naive_T_cell','gdT_cell','NK_cell',
                  'B_cell','Plasma_cell','Macrophage','DC','Fibroblast','Endothelial','Enterocyte','Immature_Enterocyte',
                  'BEST4_Colonocyte','Colonocyte','Immature_Colonocyte','TA','Cycling_TA','Goblet','Immature_Goblet','Stem',
                  'Enteroendocrine','Paneth','M_cell','Tuft')
# output directories
dir_out <- file.path(results_dir,'Downsampled_fgsea_RNAEnrich_co_viz/UC_R')
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
#ImmDict results 
gsea_ImmDict<-read.table(file_gsea_ImmDict,sep='\t',header=T,check.names=F)
gsea_ImmDict<-gsea_ImmDict[gsea_ImmDict$Context%in%order_clustersubtype,]
gsea_ImmDict_RNAEnrich<-read.table(file_gsea_ImmDict_RNAEnrich,sep='\t',header=T,check.names=F)
gsea_ImmDict_RNAEnrich<-gsea_ImmDict_RNAEnrich[gsea_ImmDict_RNAEnrich$Context%in%order_clustersubtype,]
df_order_cytokine_ImmDict<-read.table(file_df_order_cytokine_ImmDict,sep='\t',header=T)
df_order_cytokine_ImmDict<-df_order_cytokine_ImmDict[,c('Pathway','Family')]
order_cytokine_ImmDict<-df_order_cytokine_ImmDict$Pathway
gsea_ImmDict<-gsea_ImmDict[gsea_ImmDict$Pathway%in%order_cytokine_ImmDict,]
gsea_ImmDict_RNAEnrich<-gsea_ImmDict_RNAEnrich[gsea_ImmDict_RNAEnrich$Pathway%in%order_cytokine_ImmDict,]
#CytoSig results
gsea_CytoSig<-read.table(file_gsea_CytoSig,sep='\t',header=T,check.names=F)
gsea_CytoSig<-gsea_CytoSig[gsea_CytoSig$Context%in%order_clustersubtype,]
gsea_CytoSig_RNAEnrich<-read.table(file_gsea_CytoSig_RNAEnrich,sep='\t',header=T,check.names=F)
gsea_CytoSig_RNAEnrich<-gsea_CytoSig_RNAEnrich[gsea_CytoSig_RNAEnrich$Context%in%order_clustersubtype,]
df_order_cytokine_CytoSig<-read.table(file_df_order_cytokine_CytoSig,sep='\t',header=T)
df_order_cytokine_CytoSig<-df_order_cytokine_CytoSig[,c('Pathway','Family')]
order_cytokine_CytoSig<-df_order_cytokine_CytoSig$Pathway
gsea_CytoSig<-gsea_CytoSig[gsea_CytoSig$Pathway%in%order_cytokine_CytoSig,]
gsea_CytoSig_RNAEnrich<-gsea_CytoSig_RNAEnrich[gsea_CytoSig_RNAEnrich$Pathway%in%order_cytokine_CytoSig,]
#IBD enhancer targets
ImmDict_pathway_IBD_enhancer_targets_L<-readLines(file_cytokine_IBD_enhancer_targets_L_ImmDict)
ImmDict_pathway_IBD_enhancer_targets_R<-readLines(file_cytokine_IBD_enhancer_targets_R_ImmDict)
ImmDict_pathway_IBD_enhancer_targets<-unique(c(ImmDict_pathway_IBD_enhancer_targets_L,ImmDict_pathway_IBD_enhancer_targets_R))
CytoSig_pathway_IBD_enhancer_targets_L<-readLines(file_cytokine_IBD_enhancer_targets_L_CytoSig)
CytoSig_pathway_IBD_enhancer_targets_R<-readLines(file_cytokine_IBD_enhancer_targets_R_CytoSig)
CytoSig_pathway_IBD_enhancer_targets<-unique(c(CytoSig_pathway_IBD_enhancer_targets_L,CytoSig_pathway_IBD_enhancer_targets_R))

#information for bold ligand if detected
pathway_ligand_detected_ImmDict<-unique(gsea_ImmDict[gsea_ImmDict$Ligand_pct==1,'Pathway'])
pathway_ligand_detected_CytoSig<-unique(gsea_CytoSig[gsea_CytoSig$Ligand_pct==1,'Pathway'])
# pval enrichment heatmaps - pathway vs. time point
#Visualize GSEA results
df_viz_padj_ImmDict <- acast(gsea_ImmDict, Pathway ~ Context, value.var = "-log10Padj")
df_viz_padj_ImmDict<-df_viz_padj_ImmDict[order_cytokine_ImmDict,order_clustersubtype]
df_viz_padj_ImmDict[is.na(df_viz_padj_ImmDict)]<-0
df_viz_padj_ImmDict_RNAEnrich <- acast(gsea_ImmDict_RNAEnrich, Pathway ~ Context, value.var = "-log10Padj")
df_viz_padj_ImmDict_RNAEnrich<-df_viz_padj_ImmDict_RNAEnrich[order_cytokine_ImmDict,order_clustersubtype]
df_viz_padj_ImmDict_RNAEnrich[is.na(df_viz_padj_ImmDict_RNAEnrich)]<-0
df_viz_receptor_ImmDict <- acast(gsea_ImmDict, Pathway ~ Context, value.var = "Receptor_pct")
df_viz_receptor_ImmDict<-df_viz_receptor_ImmDict[order_cytokine_ImmDict,order_clustersubtype]
df_viz_receptor_ImmDict[is.na(df_viz_receptor_ImmDict)]<-0

df_viz_padj_CytoSig <- acast(gsea_CytoSig, Pathway ~ Context, value.var = "-log10Padj")
df_viz_padj_CytoSig<-df_viz_padj_CytoSig[order_cytokine_CytoSig,order_clustersubtype]
df_viz_padj_CytoSig[is.na(df_viz_padj_CytoSig)]<-0
df_viz_padj_CytoSig_RNAEnrich <- acast(gsea_CytoSig_RNAEnrich, Pathway ~ Context, value.var = "-log10Padj")
df_viz_padj_CytoSig_RNAEnrich<-df_viz_padj_CytoSig_RNAEnrich[order_cytokine_CytoSig,order_clustersubtype]
df_viz_padj_CytoSig_RNAEnrich[is.na(df_viz_padj_CytoSig_RNAEnrich)]<-0
df_viz_receptor_CytoSig <- acast(gsea_CytoSig, Pathway ~ Context, value.var = "Receptor_pct")
df_viz_receptor_CytoSig<-df_viz_receptor_CytoSig[order_cytokine_CytoSig,order_clustersubtype]
df_viz_receptor_CytoSig[is.na(df_viz_receptor_CytoSig)]<-0

meta <- data.frame(order_clustersubtype)
colnames(meta)<-c('CellType')
rownames(meta)<-meta$CellType
meta <- meta[order_clustersubtype,,drop=F]
meta[,'CellType']<-gsub('_UC_R','',meta[,'CellType'])
#Visualize the gene expression profile
name_annot <- c('CellType'='CellType')
celltype_colors<-read.delim(file.path(dir_data,"color_palette","IBD_main_color_palette.txt"),header=T,row.names=1,sep='\t')
celltype_colors<-celltype_colors[meta$CellType,,drop=F]
heat_meta_col <- list()
heat_meta_col[['CellType']] <- setNames(celltype_colors$Color, meta$CellType)

# Heatmap
print('Heatmap')
celltype_annot <- HeatmapAnnotation(df=meta, col=heat_meta_col,show_legend=FALSE,
                                    annotation_name_gp = gpar(fontsize=0,fontface=2),
                                    simple_anno_size = unit(3, "mm"))
heat_max_pval<-20
cols <- rev(brewer.pal(11, "PRGn"))  
heat_col<-colorRamp2(c(-heat_max_pval,-15,-10,-5,-3,0,3,5,10,15,heat_max_pval),cols)
lgd <- Legend(col_fun=heat_col, title='-log10(Padj)*signed(NES)',at=c(-heat_max_pval,-15,-10,-5,-3,0,3,5,10,15,heat_max_pval),labels=c('-20','-15','-10','-5','-3','0','3','5','10','15','20'),
              legend_height = unit(4, "cm"),
              labels_gp = gpar(fontsize=9, fontface = 1), title_gp = gpar(fontsize=10, fontface = 2))
row_split_ImmDict= data.frame(df_order_cytokine_ImmDict$Family)
colnames(row_split_ImmDict)<-'Family'
row_split_ImmDict$Family<-factor(row_split_ImmDict$Family,levels=unique(df_order_cytokine_ImmDict$Family))
row_split_CytoSig= data.frame(df_order_cytokine_CytoSig$Family)
colnames(row_split_CytoSig)<-'Family'
row_split_CytoSig$Family<-factor(row_split_CytoSig$Family,levels=unique(df_order_cytokine_CytoSig$Family))

label_celltype <- order_clustersubtype
label_celltype <- gsub('_UC_R','',label_celltype)
label_celltype<-gsub('CD4_T_cell','Treg/Tfh',label_celltype)
label_celltype<-gsub('CD4_T_Eff','CD4+ T Eff',label_celltype)
label_celltype<-gsub('CD8_T_cell','CD8+ T cell',label_celltype)
label_celltype<-gsub('gdT_cell','gd T cell',label_celltype)
label_celltype<-gsub('Immature_','Imm.',label_celltype)
label_celltype<-gsub('_',' ',label_celltype)
#add information for IBD enhancer targets
row_labels_CytoSig <- rownames(df_viz_padj_CytoSig)
is_bold_CytoSig <- row_labels_CytoSig %in% pathway_ligand_detected_CytoSig
is_dagger_CytoSig  <- row_labels_CytoSig %in% CytoSig_pathway_IBD_enhancer_targets
row_names_gp_CytoSig <- gpar(
  fontsize = 8,
  fontface = ifelse(is_bold_CytoSig, "bold", "plain")
)
row_labels_CytoSig <- ifelse(
  is_dagger_CytoSig, 
  paste0(row_labels_CytoSig, "<sup>i</sup>"), #press option+T on keyboard
  row_labels_CytoSig
)

row_labels_ImmDict <- rownames(df_viz_padj_ImmDict)
is_bold_ImmDict <- row_labels_ImmDict %in% pathway_ligand_detected_ImmDict
is_dagger_ImmDict  <- row_labels_ImmDict %in% ImmDict_pathway_IBD_enhancer_targets
row_names_gp_ImmDict <- gpar(
  fontsize = 8,
  fontface = ifelse(is_bold_ImmDict, "bold", "plain")
)
row_labels_ImmDict <- ifelse(
  is_dagger_ImmDict, 
  paste0(row_labels_ImmDict, "<sup>i</sup>"), #press option+T on keyboard
  row_labels_ImmDict
)
#Nsample information
df_Nsample<-read.table(file_df_Nsample,sep='\t',header=T)
df_Nsample<-df_Nsample[df_Nsample$DiseaseTissue=='UC_R',]
df_Nsample$CellTypeDiseaseTissue<-paste(df_Nsample$CellType,df_Nsample$DiseaseTissue,sep='_')
df_Nsample<-df_Nsample[df_Nsample$CellTypeDiseaseTissue%in%order_clustersubtype,]
df_Nsample_I<-df_Nsample[df_Nsample$Status=='I',]
df_Nsample_NI<-df_Nsample[df_Nsample$Status=='NI',]
df_viz_I <- df_Nsample_I %>%
  group_by(CellTypeDiseaseTissue) %>%
  summarize(
    NSample = n_distinct(Sample)
  )
df_viz_I$Status<-'NR'
df_viz_NI <- df_Nsample_NI %>%
  group_by(CellTypeDiseaseTissue) %>%
  summarize(
    NSample = n_distinct(Sample)
  )
df_viz_NI$Status<-'R'
df_viz_Nsample<-rbind(df_viz_I,df_viz_NI)
df_viz_Nsample <- acast(df_viz_Nsample, Status ~ CellTypeDiseaseTissue, value.var = "NSample")
df_viz_Nsample <-df_viz_Nsample[c('R','NR'),order_clustersubtype]
file_out <- file.path(dir_out, paste0('heat_',file_save,'_sig_cytokine_FDR10_ImmDict_CytoSig_viz_sig_fgsea_RNAEnrich_Nsample_3_filter_UC_R.pdf'))
ht1 <- Heatmap(df_viz_padj_ImmDict, col=heat_col, cluster_columns=FALSE, cluster_rows=FALSE,
               show_row_dend=FALSE, show_row_names=TRUE, show_column_names=TRUE, 
               column_names_rot=90, column_names_side='top',column_labels=gt_render(label_celltype),
               column_names_gp=gpar(fontsize=10, fontface=1),cluster_column_slices = FALSE,
               row_labels = gt_render(row_labels_ImmDict),row_names_gp=row_names_gp_ImmDict,show_heatmap_legend=FALSE,
               column_split = NULL,row_split=row_split_ImmDict, row_title_rot = 0,
               top_annotation=celltype_annot,
               row_names_side='left',row_title_gp=gpar(fontsize=10, fontface=2),
               row_gap = unit(0, "mm"),column_title=NULL, column_title_gp=gpar(fontsize=12, fontface=2),
               width = ncol(df_viz_padj_ImmDict)*unit(3, "mm"), height = nrow(df_viz_padj_ImmDict)*unit(3, "mm"),
               column_gap = unit(0, "mm"), border = TRUE,use_raster = FALSE,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 val_1<- df_viz_padj_ImmDict[i, j]
                 val_2<- df_viz_padj_ImmDict_RNAEnrich[i, j]
                 val_3<- df_viz_receptor_ImmDict[i,j]
                 if (!is.na(val_3) &&val_3 == 1&&
                     !is.na(val_1) && (val_1) > 1&& 
                     !is.na(val_2) && (val_2) > 1){
                   grid.rect(x = x, y = y,width = width,height = height,
                             gp = gpar(fill = NA,col = "black",lwd = 1.5))}
                 if (!is.na(val_3) &&val_3 == 1&&
                     !is.na(val_1) && (val_1) < (-1)&& 
                     !is.na(val_2) && (val_2) < (-1)){
                   grid.rect(x = x, y = y,width = width,height = height,
                             gp = gpar(fill = NA,col = "black",lwd = 1.5))}
                 if (!is.na(val_1) && abs(val_1) > 1&&
                     !is.na(val_2) && abs(val_2) > 1) {
                   grid.text("*", x = x, y = y-height * 0.25,
                             gp = gpar(fontsize = 9, fontface = "bold"))}
                 if (!is.na(val_1) && abs(val_1) < (-1)&&
                     !is.na(val_2) && abs(val_2) < (-1)) {
                   grid.text("*", x = x, y = y-height * 0.25,
                             gp = gpar(fontsize = 9, fontface = "bold"))}})
ht2 <- Heatmap(df_viz_padj_CytoSig, col=heat_col, cluster_columns=FALSE, cluster_rows=FALSE,
               show_row_dend=FALSE, show_row_names=TRUE, show_column_names=TRUE, 
               column_names_rot=90, column_names_side='top',column_labels=gt_render(label_celltype),
               column_names_gp=gpar(fontsize=0, fontface=1),cluster_column_slices = FALSE,
               row_labels = gt_render(row_labels_CytoSig),row_names_gp=row_names_gp_CytoSig,show_heatmap_legend=FALSE,
               column_split = NULL,row_split=row_split_CytoSig, row_title_rot = 0,
               row_names_side='left',row_title_gp=gpar(fontsize=10, fontface=2),
               row_gap = unit(0, "mm"),column_title=NULL, column_title_gp=gpar(fontsize=12, fontface=2),
               width = ncol(df_viz_padj_CytoSig)*unit(3, "mm"), height = nrow(df_viz_padj_CytoSig)*unit(3, "mm"),
               column_gap = unit(0, "mm"), border = TRUE,use_raster = FALSE,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 val_1<- df_viz_padj_CytoSig[i, j]
                 val_2<- df_viz_padj_CytoSig_RNAEnrich[i, j]
                 val_3<- df_viz_receptor_CytoSig[i,j]
                 if (!is.na(val_3) &&val_3 == 1&&
                     !is.na(val_1) && (val_1) > 1&& 
                     !is.na(val_2) && (val_2) > 1){
                   grid.rect(x = x, y = y,width = width,height = height,
                             gp = gpar(fill = NA,col = "black",lwd = 1.5))}
                 if (!is.na(val_3) &&val_3 == 1&&
                     !is.na(val_1) && (val_1) < (-1)&& 
                     !is.na(val_2) && (val_2) < (-1)){
                   grid.rect(x = x, y = y,width = width,height = height,
                             gp = gpar(fill = NA,col = "black",lwd = 1.5))}
                 if (!is.na(val_1) && abs(val_1) > 1&&
                     !is.na(val_2) && abs(val_2) > 1) {
                   grid.text("*", x = x, y = y-height * 0.25,
                             gp = gpar(fontsize = 9, fontface = "bold"))}
                 if (!is.na(val_1) && abs(val_1) < (-1)&&
                     !is.na(val_2) && abs(val_2) < (-1)) {
                   grid.text("*", x = x, y = y-height * 0.25,
                             gp = gpar(fontsize = 9, fontface = "bold"))}})
ht3 <- Heatmap(df_viz_Nsample, name = "N donors",col = colorRamp2(c(0, max(df_viz_Nsample)), c("white", "white")),  # white background
               cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE,
               show_column_names = TRUE, row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 0), row_names_side='left',
               width = ncol(df_viz_Nsample)*unit(3, "mm"), height = nrow(df_viz_Nsample)*unit(3, "mm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(df_viz_Nsample[i, j],x = x,y = y,
                gp = gpar(col = "black", fontsize = 8))},
               rect_gp = gpar(col = "grey80", lwd = 0.5),
               show_heatmap_legend=FALSE,
               border = TRUE)  # gives grid lines between cells
                                                  
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
ht_list<-ht1 %v% ht2  %v% ht3
size = calc_ht_size(ht_list)
pdf(file_out, width=size[1]+3, height=size[2])
ht_list <- draw(ht_list,ht_gap = unit(3, "mm"),padding = unit(c(2, 15, 2, 2), "mm"),annotation_legend_list = list(lgd))
dev.off()

