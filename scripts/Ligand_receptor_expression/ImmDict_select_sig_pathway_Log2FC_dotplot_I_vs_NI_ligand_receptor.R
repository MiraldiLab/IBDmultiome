rm(list=ls())
options(stringsAsFactors=FALSE)
set.seed(42)
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(magrittr)
  library(sva)
  library(RColorBrewer)
  library(scales)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(DESeq2)
  library(edgeR)
  library(progeny)
  library(ComplexHeatmap)
  library(circlize)
  library(scales)
  library(ggnewscale)
  library(ggtext)
})
project_dir <- "."
script_dir <- file.path(project_dir, "scripts")
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results","GSEA_RNA_Status")
####### INPUTS #######

# Output directory
dir_out <- file.path(results_dir,'Downsampled_fgsea_RNAEnrich_co_viz/Ligand_receptor_expression')
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

# scRNA data
file_data <- file.path(data_dir,paste('IBD_multiome_obj.rds')) #download from GEO
file_IBD_risk_gene<-file.path(data_dir,paste('/IBD_risk_loci_enhancer_gene_targets_intersect_eQTL_ABC_2kb_TSS_3D_chromatin_filter_5pct_min_expression_by_CellType.txt'))
file_CRC_risk_gene<-file.path(data_dir,paste('CRC_risk_loci_enhancer_gene_targets_intersect_eQTL_ABC_2kb_TSS_3D_chromatin_filter_5pct_min_expression_by_CellType.txt'))
data_assay <- 'RNA' # data assay
order_celltype<-c('CD4_T_cell','CD4_T_Eff','CD8_T_cell','Naive_T_cell','gdT_cell','NK_cell','ILC',
'B_cell','Plasma_cell','Macrophage','DC','Cycling_Immune','Fibroblast','Endothelial',
'Enterocyte','Immature_Enterocyte','BEST4_Colonocyte','Colonocyte','Immature_Colonocyte',
'TA','Cycling_TA','Goblet','Immature_Goblet','Stem','Enteroendocrine','Paneth','M_cell','Tuft')
order_disease_tissue<-c('CD_R','CD_TI','UC_R')
order_status<-c('NI','I')
# Params
file_save <- 'IBD' # file save base
meta_celltype <- 'CellType' # celltype metadata name
min_frac <- 0.05 # expressed in minimum fraction of cells per cluster
min_cells <- 50 # min cells to keep gene scrna
min_feat_umi<-20 # minimum counts per gene
min_group_umi <- 5E4 # Min UMIs to keep pseudobulk celltype/condition
cutoff_log2fc <- 0.58 # logFC cutoff for celltype signature genes
cutoff_fdr <- 0.1# FDR cutoff for celltype signature gene
type_de <- 'DEseq2' # type of DE test, options: 'DESeq2','wilcox', 'MAST', 'LR'
######################
# Load scRNA
print('Load scRNA')
obj <- readRDS(file_in_data)
DefaultAssay(obj) <- data_assay
Idents(obj) <- meta_celltype
obj<-subset(obj,idents=order_celltype)
order_clustersubtype<-NULL
for (ix in order_celltype){
  for (jx in order_disease_tissue){
    for (kx in order_status){
      order_clustersubtype<-c(order_clustersubtype,paste(ix,jx,kx,sep='_'))
    }
  }
}
obj@meta.data[,'ClusterSubtype']<-paste(obj@meta.data[,meta_celltype],obj@meta.data[,'DiseaseTissueStatus'],sep='_')
Idents(obj)<-'ClusterSubtype'
keep_idents<- names(which(table(Idents(obj)) >=10))
obj<-subset(obj,idents=keep_idents)
order_clustersubtype<-intersect(order_clustersubtype,keep_idents)
#ImmDict and CytoSig L-R pairs, order the same as in the heatmap
file_ImmDict_LR <- file.path(data_dir, "gene_sets", "ImmuneDictionary", "ImmDict_LR_pairs_by_cellchatDB.txt")
df_ImmDict<-read.table(file_ImmDict_LR,sep='\t',header=T,check.names=F)
order_cytokine<-c('TNFA','TL1A','CD40L','IL18','IL36','IL1A','IL1B','IFNG',
                  'IL12','IL23','IL27','IL15','IL7','IL21','IL13','IL4','IL2',
                  'TGFB1')
df_ImmDict<-df_ImmDict[df_ImmDict$Pathway%in%order_cytokine,]
order_L_R<-NULL
for (ix in order_cytokine){
  curr_L<-unique(unlist(strsplit(df_ImmDict[df_ImmDict$Pathway==ix,'Ligand'],split='_')))
  curr_R<-unique(unlist(strsplit(df_ImmDict[df_ImmDict$Pathway==ix,'Receptor'],split='_')))
  order_L_R<-unique(c(order_L_R,curr_L,curr_R))
}
#Log2FC between I vs NI
df_exp_Log2FC<-read.table(file.path(dir_out,paste('results/pseudobulk_RNA_DEseq2_Status/res_DESeq2_IBD_multiome_RNA_Disease_Tisssue_Status_shrink_ashr_unfiltered.txt')),sep='\t',header=T,check.names=F)
df_exp_Log2FC$CellTypeDiseaseTissue<-paste(df_exp_Log2FC$CellType,df_exp_Log2FC$DiseaseTissue,sep='_')
rownames(df_exp_Log2FC)<-paste(df_exp_Log2FC$CellTypeDiseaseTissue,df_exp_Log2FC$Gene,sep='_')
#load IBD/CRC risk gene
IBD_risk_gene<-readLines(file_IBD_risk_gene)
CRC_risk_gene<-readLines(file_CRC_risk_gene)
IBD_enhancer_targets_L<-intersect(IBD_risk_gene,unique(unlist(strsplit(df_ImmDict$Ligand,split='_'))))
IBD_enhancer_targets_R<-intersect(IBD_risk_gene,unique(unlist(strsplit(df_ImmDict$Receptor,split='_'))))
CRC_enhancer_targets_L<-intersect(CRC_risk_gene,unique(unlist(strsplit(df_ImmDict$Ligand,split='_'))))
CRC_enhancer_targets_R<-intersect(CRC_risk_gene,unique(unlist(strsplit(df_ImmDict$Receptor,split='_'))))
#visualize the z-scored 
counts <- LayerData(obj, layer='data')
order_markers <- order_L_R
order_markers<-intersect(order_markers,rownames(counts))
df_sc <- as.data.frame(as.matrix(counts)[order_markers,])
df_sc <- melt(as.matrix(df_sc))
colnames(df_sc) <- c('Gene','Cell','Count')
df_sc$ClusterSubtype <- as.character(Idents(obj)[as.character(df_sc$Cell)])
df_sc$PctExp <- ifelse(df_sc$Count > 0, 1, 0)
df_pct <- aggregate(PctExp ~ ClusterSubtype + Gene, df_sc, sum)
df_pct <- dcast(df_pct, ClusterSubtype ~ Gene)
rownames(df_pct) <- df_pct$ClusterSubtype
df_pct <- as.matrix(df_pct[,-1])
df_pct <- t(df_pct)
df_ncell <- table(Idents(obj))
df_ncell <- df_ncell[colnames(df_pct)]
df_pct <- t(sweep(df_pct, 2, df_ncell, FUN = '/')*100)
df_pct <- melt(df_pct)
colnames(df_pct) <- c('ClusterSubtype','Gene','PctExp')
df_pct$ClusterSubtype<-as.character(df_pct$ClusterSubtype)
df_pct$CellTypeDiseaseTissue<-sapply(strsplit(df_pct$ClusterSubtype,'_'), function(x) paste(x[1:(length(x)-1)],collapse='_'))
df_pct$DiseaseTissue<-sapply(strsplit(df_pct$CellTypeDiseaseTissue,'_'), function(x) paste(x[(length(x)-1):length(x)],collapse='_'))
df_pct$Status<-sapply(strsplit(df_pct$ClusterSubtype,'_'), function(x) x[length(x)])
df_pct<-df_pct[df_pct$DiseaseTissue%in%order_disease_tissue,]
df_pct <- df_pct %>%
  pivot_wider(
    id_cols = c(Gene, CellTypeDiseaseTissue, DiseaseTissue),
    names_from = Status,
    values_from = PctExp,
    names_prefix = "PctExp_"
  )
df_pct<-as.data.frame(df_pct)
rownames(df_pct) <- paste(as.character(df_pct$CellTypeDiseaseTissue), as.character(df_pct$Gene), sep='_')
df_dot <- df_exp_Log2FC
df_dot[rownames(df_pct),'PctExp_I'] <- df_pct$PctExp_I
df_dot[rownames(df_pct),'PctExp_NI'] <- df_pct$PctExp_NI
df_dot <- subset(df_dot, Gene %in% order_markers)
df_dot <- unique(df_dot[,c('Gene','Log2FC','pval','padj','CellTypeDiseaseTissue','PctExp_I','PctExp_NI')])
order_celltype_disease_tissue<-NULL
for (ix in order_celltype){
  for (jx in order_disease_tissue){
      order_celltype_disease_tissue<-c(order_celltype_disease_tissue,paste(ix,jx,sep='_'))
  }
}
order_celltype_disease_tissue<-intersect(order_celltype_disease_tissue,unique(df_dot$CellTypeDiseaseTissue))
df_viz<-df_dot
df_viz$CellTypeDiseaseTissue <- factor(df_viz$CellTypeDiseaseTissue, levels=order_celltype_disease_tissue)
df_viz$Gene<-factor(df_viz$Gene,levels=rev(order_markers))
df_viz$PctExp <- pmax(df_viz$PctExp_I, df_viz$PctExp_NI)
#add lines between cell types
meta<-as.data.frame(order_celltype_disease_tissue)
colnames(meta)<-'CellTypeDiseaseTissue'
meta$CellType <- sapply(strsplit(meta$CellTypeDiseaseTissue, '_'), function(x) paste(x[1:(length(x)-2)], collapse = '_'))
meta$Disease<-sapply(strsplit(meta$CellTypeDiseaseTissue,'_'),function(x)x[length(x)-1])
meta$Tissue<-sapply(strsplit(meta$CellTypeDiseaseTissue,'_'),function(x)x[length(x)])
rownames(meta)<-meta$CellTypeDiseaseTissue
order_lines<-as.data.frame(meta %>% group_by(CellType) %>% summarize(n = n()))
order_lines<-cumsum(order_lines[order(match(order_lines$CellType,order_celltype)),'n'])
label_render_y_markdown <- function(pathways) {
  sapply(pathways, function(pw) {
    label <- pw
    
    # Add symbols
    if (pw %in% CRC_enhancer_targets_L) {
      label <- paste0(label,"<sup>c</sup>")
    }
    if (pw %in% IBD_enhancer_targets_L) {
      label <- paste0(label,"<sup>i</sup>")
    }
    if (pw %in% IBD_enhancer_targets_R) {
      label <- paste0(label, "<sup>i</sup>")
    }
    if (pw %in% CRC_enhancer_targets_R) {
      label <- paste0(label, "<sup>c</sup>")
    }
    label <- paste0("<i>", label, "</i>")
    return(label)
  })
}
colors<-rev(brewer.pal(11, "PRGn"))
ggplot(df_viz, aes(y=Gene, x=CellTypeDiseaseTissue, size=PctExp, color=Log2FC)) + 
  geom_point(alpha=0.8) + theme_classic() + 
  scale_color_gradientn(
    colours = c(colors[1],colors[3],colors[5],colors[6],colors[7],colors[9],colors[11]),
    values = scales::rescale(c(-2.5,-1.5, -0.26,0, 0.26,1.5,2.5)),
    limits = c(-2.5, 2.5),
    breaks=c(-2.5,-1.5, -0.26,0, 0.26,1.5,2.5),
    oob = scales::squish
  )+
  scale_x_discrete(position='top') +
  scale_y_discrete(breaks = rev(order_markers), labels=label_render_y_markdown(rev(order_markers)))+    
  theme(axis.text.x.top=element_text(angle=45, hjust=0, vjust=0, size=11)) +
  theme(axis.text.y = ggtext::element_markdown(size = 10, color = "black"))+
  labs(x='',y='',color='Log2FC(NR vs R)',size='% Cells\nDetected') +
  theme(legend.position='bottom') +
  guides(color=guide_colourbar(title.position='top', title.hjust=-0.5))+
  theme(plot.margin = margin(10,80,10,40))+
  geom_vline(xintercept = order_lines+0.5, 
             color = "black", linetype = "solid", size = 0.5)+
  geom_point(data = subset(df_viz, PctExp_I >= 5& PctExp_NI >=5), 
           aes(y = Gene, x = CellTypeDiseaseTissue, size = PctExp), 
           color = "black", 
           shape=1,
           inherit.aes = FALSE)+
  geom_point(data = subset(df_viz, PctExp_I >= 5& PctExp_NI<5), 
           aes(y = Gene, x = CellTypeDiseaseTissue, size = PctExp), 
           color = "#762A83", 
           shape=1,
           inherit.aes = FALSE)+
    geom_point(data = subset(df_viz, PctExp_I <5& PctExp_NI>=5), 
           aes(y = Gene, x = CellTypeDiseaseTissue, size = PctExp), 
           color = "#1B7837", 
           shape=1,
           inherit.aes = FALSE)+
  geom_text(data = subset(df_viz, PctExp >= 5 & padj <= 0.1 & Log2FC > 0.26),
          aes(y = Gene, x = CellTypeDiseaseTissue, label = "*"),
          color = "black",
          position = position_nudge(y = -0.3),
          size = 6)+ # Adjust for visual appearance; `size = 6` roughly matches dot size
  geom_text(data = subset(df_viz, PctExp >= 5 & padj <= 0.1 &Log2FC<(-0.26)), 
           aes(y = Gene, x = CellTypeDiseaseTissue, label = "*"), 
           color = "black",
           position = position_nudge(y = -0.3),
           size = 6)
file_out <- file.path(dir_out, paste0('dot_markers_ImmuneDict_select_sig_cytokine_FDR10_CellTypeDiseaseTissue_Log2FC.pdf'))
ggsave(file_out, height=10.5, width=11)

meta<-meta[,c('CellType','Disease','Tissue')]
meta$Tissue <- sub("TI", "Ileum", meta$Tissue)
meta$Tissue <- sub("R", "Rectum", meta$Tissue)
df_viz<-acast(df_viz,Gene~CellTypeDiseaseTissue,value.var='Log2FC')
df_viz<-df_viz[intersect(order_markers,rownames(df_viz)),order_celltype_disease_tissue]
heat_max_padj<-0.58
heat_col = colorRamp2(c(-heat_max_padj, 0, heat_max_padj), 
                         c('dodgerblue2', 'white', 'red'))
name_annot <- c('CellType'='Cell Type', 'Disease'='Disease','Tissue'='Tissue')
celltype_colors<-read.delim(file.path(dir_data,"color_palette","IBD_main_color_palette.txt"),header=T,row.names=1,sep='\t')
celltype_colors<-celltype_colors[order_celltype,,drop=F]
heat_meta_col <- list()
heat_meta_col[['CellType']] <- setNames(celltype_colors$Color, order_celltype)
heat_meta_col[['Disease']] <- c('CD'='#9d8189',"UC"='#3d405b')
heat_meta_col[['Tissue']] <- c('Rectum'='#7F6545','Ileum'='#B59469')

lgd <- Legend(col_fun=heat_col, title='Log2FC(I vs NI)',at=c(-heat_max_padj,0,heat_max_padj),labels=c(paste0(-heat_max_padj),'0',paste0(heat_max_padj)),
              legend_height = unit(4, "cm"),
              labels_gp = gpar(fontsize=16, fontface = 1), title_gp = gpar(fontsize=18, fontface = 2))
celltype_annot <- HeatmapAnnotation(df=meta, col=heat_meta_col,show_legend=FALSE,
                                      annotation_name_gp = gpar(fontsize=0,fontface=2))

lgd_Disease= Legend(labels = c("CD",'UC'), title = "Disease", legend_gp = gpar(fill = heat_meta_col[['Disease']] ),
                   legend_height = unit(4, "cm"),legend_width = unit(4, "cm"),
                   labels_gp = gpar(fontsize=16, fontface = 1), title_gp = gpar(fontsize=18, fontface = 2))
lgd_Tissue= Legend(labels = c('Rectum','Ileum'), title = "Tissue", legend_gp = gpar(fill = heat_meta_col[['Tissue']] ),
                   legend_height = unit(4, "cm"),legend_width = unit(4, "cm"),
                   labels_gp = gpar(fontsize=16, fontface = 1), title_gp = gpar(fontsize=18, fontface = 2))

col_section<- colnames(df_viz)
for (ix in order_disease_tissue){
    col_section <- gsub(paste0('_',ix),'',col_section)
  }
col_split = data.frame(col_section)
colnames(col_split)<-'celltype'
col_split$celltype<-factor(col_split$celltype,levels=order_celltype)
# Heatmap
print('Heatmap')
# Label celltypes
label_celltype <- colnames(df_viz)
for (ix in order_disease_tissue){
    label_celltype <- gsub(paste0('_',ix),'',label_celltype)
  }
label_celltype[which(duplicated(label_celltype))] <- ''
label_celltype<-gsub('CD4_T_cell','Treg/Tfh',label_celltype)
label_celltype<-gsub('CD4_T_Eff','CD4+ T Eff',label_celltype)
label_celltype<-gsub('_',' ',label_celltype)
label_celltype<-gsub('Immature','Imm.',label_celltype)
file_out <- file.path(dir_out, paste0('heat_ImmDict_Log2FC_ligand_receptor_top_annotation.pdf'))
pdf(file_out, width=12, height=9)
ht <- Heatmap(df_viz, col=heat_col, cluster_columns=FALSE, cluster_rows=FALSE,
              show_row_dend=FALSE, show_row_names=FALSE, show_column_names=TRUE, 
              column_names_rot=90, column_names_side='top',column_labels=gt_render(label_celltype),
              column_names_gp=gpar(fontsize=12, fontface=1),cluster_column_slices = FALSE,
              top_annotation=celltype_annot, show_heatmap_legend=FALSE,
              column_split = col_split,row_names_side='left',
              row_gap = unit(0, "mm"),column_title=NULL, 
              column_gap = unit(0, "mm"), border = TRUE,use_raster = FALSE)
ht <- draw(ht,padding = unit(c(40, 2, 2, 2), "mm"),annotation_legend_list = list(lgd_Disease,lgd_Tissue,lgd))
dev.off()