# Differential Expression & Accessibility Analyses

## Signal-response gene set enrichment analysis

[`IBD_multiome_analysis_GSEA_RNA_Status.R`](../scripts/IBD_multiome_analysis_GSEA_RNA_Status.md) performs gene set enrichment analysis (GSEA) of gene expression differences between endoscopic healing groups (I vs. NI) within each cell type-disease-tissue context. The analysis contributes to **Figure 3**.

Cell type-resolved pseudobulk differential-expression statistics were generated using **DESeq2**. Genes are ranked by the DESeq2 Wald statistic, and enrichment of curated signal-response gene sets is evaluated using **fgsea**. The script can also run **RNA-Enrich**, which was used as a complementary enrichment method.

Signal-response gene sets were derived from **Immune Dictionary (ImmDict), CytoSig, PROGENy, and a CD64 response signature**.

Analyses are performed separately for CD rectum, CD terminal ileum, and UC rectum. Benjamini-Hochberg correction is applied within each signal-response resource, with **Padj < 0.1** used as the significance threshold. In the manuscript, significant signal responses were required to be independently enriched by both fgsea and RNA-Enrich with concordant directionality.
