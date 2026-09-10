# Differential Expression & Accessibility Analyses

## Signal-response gene set enrichment analysis

[`IBD_multiome_analysis_GSEA_RNA_Status.R`](../scripts/IBD_multiome_analysis_GSEA_RNA_Status.R) performs gene set enrichment analysis (GSEA) of gene expression differences between endoscopic healing groups (I vs. NI) within each cell type-disease-tissue context. The analysis contributes to **Figure 3**.

Cell type-resolved pseudobulk differential-expression statistics were generated using **DESeq2**. Genes are ranked by the DESeq2 Wald statistic, and enrichment of curated signal-response gene sets is evaluated using **fgsea**. The script can also run **RNA-Enrich**, which was used as a complementary enrichment method.

Signal-response gene sets were derived from **Immune Dictionary (ImmDict), CytoSig, PROGENy, and a CD64 response signature**.

Analyses are performed separately for CD rectum, CD terminal ileum, and UC rectum. Benjamini-Hochberg correction is applied within each signal-response resource, with **Padj < 0.1** used as the significance threshold. In the manuscript, significant signal responses were required to be independently enriched by both fgsea and RNA-Enrich with concordant directionality.


### Signal-response gene set processing and curation

[`process_ImmuneDictionary_Cui2024_combine_pval.R`](../scripts/process_ImmuneDictionary_Cui2024_combine_pval.R)
Immune Dictionary cytokine-response signatures were processed from cell type-resolved differential expression results. For each cytokine–gene pair, signed Z-scores were calculated from the adjusted *P* values and direction of effect and combined across immune cell populations using an unweighted Stouffer Z-score. Positively regulated genes were ranked by the resulting meta-analysis *P* value, and the top genes were retained as cytokine-response gene sets.


