
# Multiome IBD cell annotation - Rectum

We performed sample integration, dimension reduction and clustering on samples separately based on tissue. Below is cell-type annotation for cells from rectum with quality RNA and ATAC signal. We follow standard Seurat workflows and final clusters are based on weighted nearest neighbor graphs.

PMIDs: [35176508](https://pubmed.ncbi.nlm.nih.gov/35176508/), [34497389](https://pubmed.ncbi.nlm.nih.gov/34497389/)


## All cells

Clustering all cells based on WNN graph (res=0.8):

<img src="./figs/multiome_IBD_R_annotation_All_clusters.png" width=60% height=60%>


## T cells

T cell subclustering (res=0.5):

<img src="./figs/multiome_IBD_R_annotation_Tcell_clusters.png" width=60% height=60%>

Quality control metrics per cluster:

<img src="./figs/multiome_IBD_R_annotation_Tcell_qc_metrics.png" width=100% height=100%>

Cell cycle scores:

<img src="./figs/multiome_IBD_R_annotation_Tcell_markers_Cycling.png" width=60% height=60%>

Conventional T cell markers:

<img src="./figs/multiome_IBD_R_annotation_Tcell_markers_Tconv.png" width=100% height=100%>

Naive T cell markers:

<img src="./figs/multiome_IBD_R_annotation_Tcell_markers_naive.png" width=60% height=60%>

Gamma-delta T cell markers:

<img src="./figs/multiome_IBD_R_annotation_Tcell_markers_gdT.png" width=80% height=80%>

Treg markers:

<img src="./figs/multiome_IBD_R_annotation_Tcell_markers_Treg.png" width=100% height=100%>

MAIT cell markers:

<img src="./figs/multiome_IBD_R_annotation_Tcell_markers_MAIT.png" width=100% height=100%>

CD8+ cytotoxic cell markers:

<img src="./figs/multiome_IBD_R_annotation_Tcell_markers_CD8_CTL.png" width=80% height=80%>

NK cell markers:

<img src="./figs/multiome_IBD_R_annotation_Tcell_markers_NK.png" width=100% height=100%>

ILC3 cell markers:

<img src="./figs/multiome_IBD_R_annotation_Tcell_markers_ILC3.png" width=100% height=100%>

Cluster 0 has low UMIs and genes and express few unique markers. Cluster 11 express epithelial markers and has high UMIs/genes. Cluster 19 expresses B cell markers. Cluster 23 expressed endothelial markers. These cells are low quality and/or cell doublets, and I marked them for removal.


## B cells

B cell subclustering (res=0.8):

<img src="./figs/multiome_IBD_R_annotation_Bcell_clusters.png" width=60% height=60%>

Quality control metrics per cluster:

<img src="./figs/multiome_IBD_R_annotation_Bcell_qc_metrics.png" width=80% height=80%>

Cell cycle scores:

<img src="./figs/multiome_IBD_R_annotation_Bcell_markers_Cycling.png" width=60% height=60%>

Memory B cell markers:

<img src="./figs/multiome_IBD_R_annotation_Bcell_markers_Bmem.png" width=100% height=100%>

Follicular B cell markers:

<img src="./figs/multiome_IBD_R_annotation_Bcell_markers_FOB.png" width=100% height=100%>

Cluster 8 expresses T cell markers. Cluster 11 expresses epithelial markers. I've marked these cells for removal.


## Epithelial

Epithelial cell subclustering (res=0.5):

<img src="./figs/multiome_IBD_R_annotation_Epithelial_clusters.png" width=60% height=60%>

Quality control metrics per cluster:

<img src="./figs/multiome_IBD_R_annotation_Epithelial_qc_metrics.png" width=100% height=100%>

Cell cycle scores:

<img src="./figs/multiome_IBD_R_annotation_Epithelial_markers_Cycling.png" width=80% height=80%>

Colonocyte cell markers:

<img src="./figs/multiome_IBD_R_annotation_Epithelial_markers_Colonocyte.png" width=100% height=100%>

TA cell markers:

<img src="./figs/multiome_IBD_R_annotation_Epithelial_markers_TA.png" width=100% height=100%>

Crypt/stem cell markers:

<img src="./figs/multiome_IBD_R_annotation_Epithelial_markers_Crypt.png" width=80% height=80%>

Goblet cell markers:

<img src="./figs/multiome_IBD_R_annotation_Epithelial_markers_Goblet.png" width=100% height=100%>

M-cell markers:

<img src="./figs/multiome_IBD_R_annotation_Epithelial_markers_M_cell.png" width=100% height=100%>

Cluster 7 expresses B cell markers. Cluster 8 expresses T cell markers. I've marked these low quality cells for removal.


## Fibroblast

Fibroblast cell subclustering (res=0.8):

<img src="./figs/multiome_IBD_R_annotation_Fibroblast_clusters.png" width=60% height=60%>

Quality control metrics per cluster:

<img src="./figs/multiome_IBD_R_annotation_Fibroblast_qc_metrics.png" width=100% height=100%>

Cell cycle scores:

<img src="./figs/multiome_IBD_R_annotation_Fibroblast_markers_Cycling.png" width=60% height=60%>

Fibroblasts expressing follicular dendritic cell (FDC) markers:

<img src="./figs/multiome_IBD_R_annotation_Fibroblast_markers_FDC.png" width=100% height=100%>

Fibroblast reticular cell markers:

<img src="./figs/multiome_IBD_R_annotation_Fibroblast_markers_Tret.png" width=100% height=100%>

Myofibroblast markers:

<img src="./figs/multiome_IBD_R_annotation_Fibroblast_markers_Myofibroblast.png" width=100% height=100%>

Clusters 4 and 9 express epithelial markers. Cluster 3 expresses T cell markers. I've marked these low quality cells for removal.


## Myeloid

Myeloid cell subclustering (res=0.5):

<img src="./figs/multiome_IBD_R_annotation_Myeloid_clusters.png" width=60% height=60%>

Quality control metrics per cluster:

<img src="./figs/multiome_IBD_R_annotation_Myeloid_qc_metrics.png" width=100% height=100%>

Cell cycle scores:

<img src="./figs/multiome_IBD_R_annotation_Myeloid_markers_Cycling.png" width=60% height=60%>

Conventional DC1 (cDC1) markers:

<img src="./figs/multiome_IBD_R_annotation_Myeloid_markers_cDC1.png" width=80% height=80%>

cDC2 markers:

<img src="./figs/multiome_IBD_R_annotation_Myeloid_markers_cDC2.png" width=80% height=80%>

Plasmacytoid DC (pDC) markers:

<img src="./figs/multiome_IBD_R_annotation_Myeloid_markers_pDC.png" width=80% height=80%>

Activated DC markers:

<img src="./figs/multiome_IBD_R_annotation_Myeloid_markers_DC_Act.png" width=80% height=80%>

Macrophage markers:

<img src="./figs/multiome_IBD_R_annotation_Myeloid_markers_MP.png" width=100% height=100%>

Cluster 3 expresses epithelial markers. Cluster 6 expresses markers of B and T cells. I've marked these low quality cells for removal.
