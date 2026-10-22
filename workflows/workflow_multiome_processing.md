# Multiome-seq Processing

We processed multiome samples primarily in R v4.2.3 with Seurat v5 and Signac v1.9.0.

Here we describe sample processing steps including: mapping, feature counts, quality control, peak calling, and sample integration.

---

## Raw data processing

For each multiome sample, we performed alignment, initial joint cell barcode calling/filtering, RNA UMI counting, ATAC transposase cut site identification, and initial peak calling using [Cell Ranger ARC v2.0.0](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/algorithms/overview) and the prebuilt [GRCh38-2020-A-2.0.0 reference](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build).

The Cell Ranger genome reference includes gene annotations from [GENCODE V28](https://www.gencodegenes.org/human/release_28.html).

Example ```cellranger-arc count``` call on the HPC cluster:

```
module load cellranger-arc/2.0.0

cellranger-arc count --id=Denson-ADL0019-TI-20220525-Multiv1 \
--reference=<path_to_Cell_Ranger_reference>/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--libraries=<path_to_libraries_csv>/lib_cellranger_arc_multiome_<Sample>.csv \
--jobmode=lsf --maxjobs=12 --disable-ui
```

Separately, we processed the RNA portion of each multiome sample using [cellranger count v7.0.1](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) to generate the input for ambient transcript removal package CellBender (described below).

Example ```cellranger count``` call on the HPC cluster:

```
module load cellranger/7.0.1

cellranger count --id=<Sample> \
--transcriptome=<path_to_Cell_Ranger_reference>/refdata-gex-GRCh38-2020-A \
--fastqs=<path_to_Sample_fastq_directories>/<Sample_RNA_fastq_directory> \
--sample=<Sample> \
--chemistry=ARC-v1 --nosecondary --jobmode=lsf --maxjobs=12 --disable-ui
```

---

## Ambient RNA removal

We removed background RNA contamination using [CellBender](https://github.com/broadinstitute/CellBender) generating an adjusted raw counts matrix for each sample.

Example ```cellbender remove-background``` call running on a GPU on the HPC cluster:

```
#BSUB -q gpu-v100
#BSUB -gpu "num=1"

module load apptainer
apptainer run --nv \
	--bind <path_to_cellranger_count_Sample_directory>/outs \
	--bind <path_to_CellBender_output_directory>/cellbender_030 \
	docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.0 \
	cellbender remove-background \
	--input <path_to_cellranger_count_Sample_directory>/outs/raw_feature_bc_matrix.h5 \
	--output <path_to_CellBender_output_directory>/cellbender_030/<Sample>/<Sample>.h5 \
	--cuda \
	--epochs 150
```

---

## [Gene filtering](./scripts/IBD_multiome_data_gene_filter.R)

We filter genes for downstream analysis based on the following criteria:
* Keep genes encoded on autosomal, chrX, and chrM chromosomes.
* Keep genes annotated as [GENCODE/Ensembl biotype](https://useast.ensembl.org/info/genome/genebuild/biotypes.html) protein coding, lncRNAs, and TCR and immunoglobulin constant chain genes.
* Remove genes classified as pseudogenes.
* Remove hemoglobin genes.

We extracted these genes from the Cell Ranger gene annotation gtf file. We downloaded [human pseudogenes](./data/gene_pseudogene_hg38_Ensembl_v109_GRCh38_p13_processed.txt) from [Ensembl v109 GRCh30.p13](http://useast.ensembl.org/Homo_sapiens/Info/Index).

---

## Quality control

### Cell doublet identification

We identified cell doublets based on RNA or ATAC, separately. For cells with both quality RNA and ATAC signal, we removed cells determined to be doublets by both methods. For cells with only quality ATAC, we removed cell doublets identified by the ATAC-based method.

#### [ATAC cell doublets](./scripts/IBD_multiome_qc_atac_doublets_ArchR.R)

We identified cell doublets from ATAC using [ArchR](https://github.com/GreenleafLab/ArchR) setting the doublet filter ratio to 0.5.

#### [RNA cell doublets](./scripts/IBD_multiome_qc_rna_doublets.R)

We identified cell doublets from RNA using [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder/tree/master), assuming a doublet rate of [0.8% per 1000 cells](https://kb.10xgenomics.com/hc/en-us/articles/360054599512-What-is-the-cell-multiplet-rate-when-using-the-3-CellPlex-Kit-for-Cell-Multiplexing-) recovered and identifying pK by running ```paramSweep_v3```.

### [Filter cell barcodes](./scripts/IBD_multiome_qc_metrics.R)

For each sample, we filtered cell barcodes called by Cell Ranger ARC based on the following metrics indicating quality transcript and accessibility signal:
* RNA: total UMIs, total genes detected, percent mitochondrial gene expression
* ATAC: total fragments, TSS enrichment, percent mitochondrial fragments, percent blacklist fragments, nucleosome signal, percent fragments in peaks

We identified singlet cell barcodes with quality RNA, quality ATAC, and those with both quality RNA and ATAC.

### [Filtered Seurat objects](./scripts/IBD_multiome_qc_filter.R)

For each sample, we generated Seurat objects containing cells with either quality ATAC or quality RNA and ATAC for downstream processing and analysis.

---

## Peak calling

### [Fragment file modification](./scripts/IBD_multiome_process_atac_fragments.R)

The [Cell Ranger fragment file](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments) specifies genomic coordinates for each ATAC fragment. We modified the fragment file of each sample to specify the transposase cut sites at each end of the fragment. We used these modified fragment files for downstream peaking calling with MACS2 (described below).

We kept only fragments from cells passing QC, and we removed fragments overlapping [blacklist regions](./data/blacklist_hg38_ENCODE_maxATAC_extended.bed). We used the extended hg38 blacklist curated by [Cazares et al.](https://pubmed.ncbi.nlm.nih.gov/36719906/) for [maxATAC](https://github.com/MiraldiLab/maxATAC).

### [Call peaks per sample](./scripts/IBD_multiome_process_atac_peaks_sample.R)

We called peaks per sample using modified fragment files and [MACS2](https://github.com/macs3-project/MACS). We broadly clustered cells with quality ATAC using the standard Seurat workflow. We ensured a minimum of 3 million fragments per cluster by merging small clusters with nearest neighbors based on correlation of pseudobulk accessibility. Example MACS2 call with parameter values:

```
macs2 callpeak -t <fragments.bed> -f BED -g 2.9e9 --keep-dup all --nomodel --llocal 10000 --shift -100 --extsize 200 -q 0.001
```

### [Reference peak mapping](./scripts/IBD_multiome_process_atac_peaks_ref.R)

We used the union of cluster-specific peaks across all samples as reference peaks for downstream processing of ATAC. For each sample, we mapped fragments to these reference peaks creating new sample-specific Seurat objects for downstream integration.

---

## Sample integration

We integrated samples using [Harmony](https://github.com/immunogenomics/harmony).

### [ATAC integration](./scripts/IBD_multiome_process_integrate_atac.R)

For cells with quality ATAC, we applied Harmony on the LSI dimensional reduction and constructed a UMAP from the Harmony dimensional reduction.

### [Multimodal integration](./scripts/IBD_multiome_process_integrate_multiome.R)

For cells with quality RNA and ATAC, we applied Harmony on each modality separately. For RNA, we ran Harmony on the PCA dimensional reduction. We then construct a [Seurat weighted nearest neighbor (WNN) graph](https://pubmed.ncbi.nlm.nih.gov/34062119/) from both Harmony dimensional reductions for downstream clustering.

### [Sub-population multimodal integration](./scripts/IBD_multiome_process_integrate_subset.R)

After broadly annotating cells, we further subset populations and re-integrated samples. Clustering on these sub-populations was used for annotation of finer cell types/states.

---

## Cell clustering

We clustered cells with quality ATAC and RNA using Leiden clustering on the Seurat WNN graph. The clustering resolution parameter was selected based on visual inspection and marker gene expression.

### [Clustering all cells](./scripts/IBD_multiome_analysis_cluster.R)

### [Clustering sub-populations](./scripts/IBD_multiome_analysis_cluster_subset.R)

---

## Cell type annotation

### Label transfer

To guide cell-type annotation, we performed Seurat anchor-based label transfer from previously published gut scRNA-seq atlases from [Smillie et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31348891/), [Elmentaite et al. 2021](https://pubmed.ncbi.nlm.nih.gov/34497389/), and [Hickey et al. 2023](https://pubmed.ncbi.nlm.nih.gov/37468586/).

For each sample, we transfered annotations to cells with quality RNA and ATAC using the RNA modality.

#### [Pan-cell-type label transfer](./scripts/IBD_multiome_analysis_transfer.R)

#### [Sub-population label transfer](./scripts/IBD_multiome_analysis_transfer_subset.R)

We annotated cell types/states based on marker gene expression and label transfer from publicly available scRNA-seq datasets from human gut.

### [Apply Cluster -> Celltype annotation](./scripts/IBD_multiome_analysis_annotation.R)
