# Differential Expression & Accessibility Analyses

## Signal-response gene set enrichment analysis

[`IBD_multiome_analysis_GSEA_RNA_Status.R`](../scripts/IBD_multiome_analysis_GSEA_RNA_Status.R) performs gene set enrichment analysis (GSEA) of gene expression differences between endoscopic healing groups (I vs. NI) within each cell type-disease-tissue context. The analysis contributes to **Figure 3**.

Cell type-resolved pseudobulk differential-expression statistics were generated using **DESeq2**. Genes are ranked by the DESeq2 Wald statistic, and enrichment of curated signal-response gene sets is evaluated using **fgsea**. The script can also run **RNA-Enrich**, which was used as a complementary enrichment method.

Signal-response gene sets were derived from **Immune Dictionary (ImmDict), CytoSig, PROGENy, and a CD64 response signature**.

Analyses are performed separately for CD rectum, CD terminal ileum, and UC rectum. Benjamini-Hochberg correction is applied within each signal-response resource, with **Padj < 0.1** used as the significance threshold. In the manuscript, significant signal responses were required to be independently enriched by both fgsea and RNA-Enrich with concordant directionality.


### Signal-response gene set processing and curation

[`process_ImmuneDictionary_Cui2024_combine_pval.R`](../scripts/process_ImmuneDictionary_Cui2024_combine_pval.R)
Immune Dictionary cytokine-response signatures were processed from cell type-resolved differential expression results. For each cytokine–gene pair, signed Z-scores were calculated from the adjusted *P* values and direction of effect and combined across immune cell populations using an unweighted Stouffer Z-score. Positively regulated genes were ranked by the resulting meta-analysis *P* value, and the top genes were retained as cytokine-response gene sets.

### Signal-response enrichment results process and visualization

Following the GSEA and RNA-Enrich analysis of signal response, process the results to incorporate nominal expression of corresponding ligand and receptor complexes (detected in >=5% nuclei in any disease, tissue, response group for a given cell type) information. For **fgsea** results following the scripts here: [`Downsampled_fgsea_Wald_statistics`](../scripts/Downsampled_fgsea_Wald_statistics); **RNA-Enrich** results following the scripts here: [`Downsampled_RNAEnrich`](../scripts/Downsampled_RNAEnrich). 

To reproduce the main and supplemental figure decks, firstly order the signals based on cytokine family and maximum enrichment scores from fgsea for significantly enriched signals in both fgsea and RNA-Enrich results, followed by the orders in CD, rectum, add the extra signals to the corresponding family for CD, TI or UC,  rectum visualization. Following the scripts here: 

For heatmap visualizations, enrichments were shown as signed -log10Padj (upregulated in TNFi NR >0, upregulated in TNFi response <0) from fgsea in the heatmaps, only significantly enriched signals from both RNA-Enrich and fgsea results were included and marked with asterisk, whereas significantly enriched signals exclusively in fsgea but not RNA-Enrich were adjusted to 0 for -log10Padj. 
Ligands complexes with nominal expression of all subunits in any cell type, response groups in the given disease, tissue context were labelled in bold, similarly receptor complexes with nominal expression of subunits in the corresponding cell type, disease, tissue, and upregulated response group were highlighted in black box. 
Signals with either ligand or receptor as IBD risk loci enhancer gene targets based on gene proximal(+/- 2kb of TSS) or distal(three-dimensional experimental data such as promoter-capture HiC, microC, TRAC-loop) from matched cell type predictions were labelled with superscript i.

To reproduce the visualization of 20 shared signals between ImmDict and CytoSig database in Fig 3.A, following the scripts here:

To reproduce the visualization of PROGENy results in Fig 3. C, following the scripts here:

To reproduce the scatter plot visualization of comparing same or closely related cell type enrichment results in Fig S4. D-G, following the scripts here: [`ScatterPlot`](../scripts/ScatterPlot)

To reproduce the visualization of all the significantly enriched (**Padj < 0.1**) ImmDict and CytoSig signal responses 


