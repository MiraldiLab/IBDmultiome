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
file_gsea_PROGENy_CD64<-file.path(results_dir,paste('Downsample_statWald/PROGENy_CD64_ImmDict/UC_R/PROGENy_CD64_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_gsea_PROGENy_CD64_RNAEnrich<-file.path(results_dir,paste('Downsample_RNAEnrich/PROGENy_CD64_ImmDict/UC_R/PROGENy_CD64_I_vs_NI_-log10_Padj_all_cytokine_df_frac_L_R_by_pathway.txt'))
file_order_cytokine_PROGENy<-file.path(results_dir,paste('Downsampled_fgsea_RNAEnrich_co_viz/UC_R/sig_sets_PROGENy_cytokine_order_by_max_log10Padj.txt'))
file_keep_context<-file.path(data_dir,'keep_context/keep_context_UC_R_Nsample_3_min_filter.txt')
data_assay <- 'RNA' # data assay
order_disease_tissue<-c('UC_R')
order_celltype<-c('CD4_T_cell','CD4_T_Eff','CD8_T_cell','Naive_T_cell','gdT_cell','NK_cell',
'B_cell','Plasma_cell','Macrophage','DC','Fibroblast','Endothelial','Enterocyte','Immature_Enterocyte',
'BEST4_Colonocyte','Colonocyte','Immature_Colonocyte','TA','Cycling_TA','Goblet','Immature_Goblet','Stem',
'Enteroendocrine','Paneth','M_cell','Tuft')
#####
# Params
file_save <- 'IBD_multiome' # file save base
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
gsea_PROGENy_CD64<-read.table(file_gsea_PROGENy_CD64,sep='\t',header=T,check.names=F)
gsea_PROGENy_CD64<-gsea_PROGENy_CD64[gsea_PROGENy_CD64$Context%in%order_clustersubtype,]
gsea_PROGENy_CD64_RNAEnrich<-read.table(file_gsea_PROGENy_CD64_RNAEnrich,sep='\t',header=T,check.names=F)
gsea_PROGENy_CD64_RNAEnrich<-gsea_PROGENy_CD64_RNAEnrich[gsea_PROGENy_CD64_RNAEnrich$Context%in%order_clustersubtype,]
order_cytokine_PROGENy<-readLines(file_order_cytokine_PROGENy)
gsea_PROGENy_CD64<-gsea_PROGENy_CD64[gsea_PROGENy_CD64$Pathway%in%order_cytokine_PROGENy,]
gsea_PROGENy_CD64_RNAEnrich<-gsea_PROGENy_CD64_RNAEnrich[gsea_PROGENy_CD64_RNAEnrich$Pathway%in%order_cytokine_PROGENy,]
#information for bold ligand if detected
pathway_ligand_detected_PROGENy_CD64<-unique(gsea_PROGENy_CD64[gsea_PROGENy_CD64$Ligand_pct==1,'Pathway'])
# pval enrichment heatmaps - pathway vs. time point
#Visualize GSEA results
df_viz_padj_PROGENy_CD64 <- acast(gsea_PROGENy_CD64, Pathway ~ Context, value.var = "-log10Padj")
df_viz_padj_PROGENy_CD64<-df_viz_padj_PROGENy_CD64[order_cytokine_PROGENy,order_clustersubtype]
df_viz_padj_PROGENy_CD64[is.na(df_viz_padj_PROGENy_CD64)]<-0
df_viz_padj_PROGENy_CD64_RNAEnrich <- acast(gsea_PROGENy_CD64_RNAEnrich, Pathway ~ Context, value.var = "-log10Padj")
df_viz_padj_PROGENy_CD64_RNAEnrich<-df_viz_padj_PROGENy_CD64_RNAEnrich[order_cytokine_PROGENy,order_clustersubtype]
df_viz_padj_PROGENy_CD64_RNAEnrich[is.na(df_viz_padj_PROGENy_CD64_RNAEnrich)]<-0
df_viz_receptor_PROGENy_CD64 <- acast(gsea_PROGENy_CD64, Pathway ~ Context, value.var = "Receptor_pct")
df_viz_receptor_PROGENy_CD64<-df_viz_receptor_PROGENy_CD64[order_cytokine_PROGENy,order_clustersubtype]
df_viz_receptor_PROGENy_CD64[is.na(df_viz_receptor_PROGENy_CD64)]<-0

meta <- data.frame(order_clustersubtype)
colnames(meta)<-c('CellType')
rownames(meta)<-meta$CellType
meta <- meta[order_clustersubtype,,drop=F]
meta[,'CellType']<-gsub('_UC_R','',meta[,'CellType'])
#Visualize the gene expression profile
name_annot <- c('CellType'='CellType')
celltype_colors<-read.delim('/data/GastroAI/genomics/analysis/multiome/liver/IBD_analysis/figures/V2/color_palette/IBD_main_color_palette.txt',header=T,row.names=1,sep='\t')
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
label_celltype <- order_clustersubtype
label_celltype <- gsub('_UC_R','',label_celltype)
label_celltype<-gsub('CD4_T_cell','Treg/Tfh',label_celltype)
label_celltype<-gsub('CD4_T_Eff','CD4+ T Eff',label_celltype)
label_celltype<-gsub('gdT_cell','gd T cell',label_celltype)
label_celltype<-gsub('Immature_','Imm.',label_celltype)
label_celltype<-gsub('Cycling_TA','Cycling TA',label_celltype)
#add information for IBD enhancer targets
row_labels_PRGOENy_CD64 <- rownames(df_viz_padj_PROGENy_CD64)
is_bold_PROGENy_CD64 <- row_labels_PRGOENy_CD64 %in% pathway_ligand_detected_PROGENy_CD64
row_names_gp_PROGENy_CD64 <- gpar(
  fontsize = 8,
  fontface = ifelse(is_bold_PROGENy_CD64, "bold", "plain")
)

file_out <- file.path(dir_out, paste0('heat_',file_save,'_sig_cytokine_FDR10_PROGENy_CD64_viz_sig_fgsea_RNAEnrich_Nsample_3_filter_UC_R.pdf'))
ht<- Heatmap(df_viz_padj_PROGENy_CD64, col=heat_col, cluster_columns=FALSE, cluster_rows=FALSE,
               show_row_dend=FALSE, show_row_names=TRUE, show_column_names=TRUE, 
               column_names_rot=90, column_names_side='top',column_labels=gt_render(label_celltype),
               column_names_gp=gpar(fontsize=10, fontface=1),cluster_column_slices = FALSE,
               row_names_gp=row_names_gp_PROGENy_CD64,show_heatmap_legend=FALSE,
               column_split = NULL,row_split=NULL,
               row_names_side='left',row_title_gp=gpar(fontsize=10, fontface=2),
               row_gap = unit(0, "mm"),column_title=NULL, column_title_gp=gpar(fontsize=12, fontface=2),
               width = ncol(df_viz_padj_PROGENy_CD64)*unit(3, "mm"), height = nrow(df_viz_padj_PROGENy_CD64)*unit(3, "mm"),
               top_annotation=celltype_annot, 
               column_gap = unit(0, "mm"), border = TRUE,use_raster = FALSE,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 val_1<- df_viz_padj_PROGENy_CD64[i, j]
                 val_2<- df_viz_padj_PROGENy_CD64_RNAEnrich[i, j]
                 val_3<- df_viz_receptor_PROGENy_CD64[i,j]
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
ht_list <- draw(ht_list,ht_gap = unit(3, "mm"),padding = unit(c(2, 15, 2, 2), "mm"),annotation_legend_list = list(lgd))
dev.off()

