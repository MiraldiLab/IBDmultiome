# Genetic Risk Integration

## RELI GWAS enrichment analyses

Cell type, tissue, therapy response-resolved IBD patients and cell type-resolved PBMC control pseudobulk ATAC-seq peaks were generated using **MACS2**, only pesudobulk library with at least 5M reads were kept for downstream genetic analysis based on GM12878 downsampled benchmark results as shown in Fig S10.

Phenotype GWAS data were derived from [**GWAS Catalog**](https://www.ebi.ac.uk/gwas/), a set of 465 human diseases and phenotypes (excluding traits primarily reflecting socioeconomic status, medication usage or
pleiotropy) that had ≥20 independent risk loci (i.e., loci for which tag variants are in linkage disequilibrium (LD) r2<0.2), each with robust genome-wide genetic association (p<5x10-8).

[**RELI**](https://github.com/WeirauchLab/RELI) enrichment analysis are performed separately for East Asian (EAS) and European (EUR) ancestry. Benjamini-Hochberg correction is applied within each pseudobulk ATAC-seq library, with **Padj < 0.1** and >= 5 independent risk loci overlaps used as the significance threshold.

Relative enrichment across contexts were visualized as z-scored enrichment scores across contexts as shown in Fig 2B for 171 significantly enriched phenotypes in EUR anceestry and selected IBD, Colorectal cancer (CRC) phenotypes across EUR and EAS ancestry.
The full set of significantly enriched results were shown in Fig. S3 with rows labelled

## Enhancer-to-gene mapping

Gene targets of IBD, CRC genetic risk loci overlaps with IBD ATAC-seq peaks were predicted through the following ways, including:
1. Activity-by-Contact (ABC) model
2. eQTL Catalogue integration
3. Curated 3D chromatin interaction data in related cell types
4. Promoter proximity mapping

## Robustness cytokine response enrichments of IBD, CRC enhancer gene targets

