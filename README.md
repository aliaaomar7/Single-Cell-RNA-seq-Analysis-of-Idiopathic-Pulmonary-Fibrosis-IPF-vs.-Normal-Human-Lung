# Single-Cell-RNA-seq-Analysis-of-Idiopathic-Pulmonary-Fibrosis-IPF-vs.-Normal-Human-Lung

Dataset: GSE132771 · Platform: 10X Genomics Chromium · Tool: Seurat v4 (R)

Overview
Idiopathic Pulmonary Fibrosis (IPF) is a progressive, fatal lung disease characterized by excessive fibrotic remodeling of the lung parenchyma. Despite its severity, the cell-type-specific transcriptomic changes that drive IPF pathology remain incompletely understood at single-cell resolution.
This project performs an end-to-end single-cell RNA-seq (scRNA-seq) analysis comparing lung tissue from 3 IPF patients and 3 normal (control) donors. The goal was to integrate multi-donor data, resolve the major lung cell populations, and pinpoint the transcriptomic signature of pathological fibroblast activation in IPF — with a focus on identifying actionable disease markers.

Data
GroupDonorsGEO AccessionsNormal3GSM3891621, GSM3891623, GSM3891625IPF3GSM3891627, GSM3891629, GSM3891631
Raw data downloaded directly from NCBI GEO as 10X-format files (matrix.mtx, barcodes.tsv, features.tsv).

Pipeline
1. Data Loading & Quality Control

All 6 samples loaded via Read10X() and converted to Seurat objects
Minimum thresholds: 3 cells per gene, 200 genes per cell
QC filters applied per sample:

Genes per cell: 500 – 4,000 (IPF1: up to 6,000 to accommodate disease-related transcriptional complexity)
UMI counts: < 20,000 (IPF1: < 40,000)
Mitochondrial content: < 10%


Each sample subsampled to a maximum of 3,000 cells to balance computational load and donor contribution

2. Normalization & Feature Selection

Log-normalization applied per sample (NormalizeData)
2,000 highly variable genes selected per sample using variance-stabilizing transformation (VST)

3. Integration

Cross-sample integration performed using Canonical Correlation Analysis (CCA) via Seurat's anchor-based framework (FindIntegrationAnchors → IntegrateData)
Corrects for donor-level batch effects while preserving biological variation

4. Dimensionality Reduction & Clustering

Data scaled and PCA run on 50 components
Elbow plot used to select top 20 PCs for downstream analysis
UMAP embedding computed for visualization
Graph-based clustering at resolution 0.1 (FindNeighbors + FindClusters)

5. Cell Type Annotation
Clusters annotated using canonical lineage markers:
Cell TypeMarkers UsedMonocyte-derived MacrophagesS100A8, FCN1EndothelialEMCN, PECAM1, CLDN5EpithelialEPCAM, SFTPC, AGER, SCGB1A1NK CellsKLRD1, PRF1FibroblastsCOL1A2, DCN, LUM, PDGFRADendritic CellsCD1C, FCER1APlasma B CellsIGHG1, IGKCMast CellsTPSB2
8 major lung cell populations resolved and annotated in the final UMAP.
6. Differential Expression — Fibroblasts (IPF vs. Normal)

DE analysis restricted to the Fibroblast cluster
Method: Wilcoxon rank-sum test via FindMarkers
Thresholds: |log2FC| > 0.5, adj. p-value < 0.05, min.pct = 0.25
Significant genes classified as Up in IPF or Up in Normal
Top hits visualized on a volcano plot

7. Pathway Enrichment

Upregulated IPF fibroblast genes submitted to EnrichR
Databases queried: GO Biological Process 2023, KEGG 2021 Human
Immunoglobulin genes (IGH/IGK/IGL) removed prior to enrichment to reduce noise


Key Findings
CTHRC1 emerged as the top upregulated gene in IPF fibroblasts — a well-established pathological myofibroblast marker that is virtually absent in normal lung fibroblasts. Its exclusive expression in IPF on the UMAP split plot confirms disease-specific fibroblast activation.
GO enrichment of IPF-upregulated fibroblast genes revealed significant enrichment in:

Extracellular Matrix Organization
Collagen Fibril Organization
Supramolecular Fiber Organization
Regulation of Apoptotic Process

These pathways collectively reflect the hallmark fibrotic program: excessive ECM deposition, collagen crosslinking, and resistance to apoptosis in activated myofibroblasts.

Output Files
FileDescriptionCONIPF_integrated.rdsIntegrated Seurat object (pre-annotation)CONIPF_annotated.rdsFinal annotated Seurat objectcell_proportion_CON_vs_IPF.pdfCell type composition bar chartfeatureplot_markers.pdfLineage marker feature plotsCTHRC1_featureplot.pdfCTHRC1 expression split by conditionGO_enrichment.pdfGO Biological Process enrichment chart

Dependencies
rlibrary(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(enrichR)
R version 4.x recommended. Set seed 42 used throughout for reproducibility.



Subsampling to 3,000 cells/donor was applied before integration to prevent dominant donors from skewing cluster structure.
IPF1 QC thresholds were widened (nFeature < 6,000; nCount < 40,000) to retain biologically relevant high-complexity cells characteristic of fibrotic tissue.
The genes.tsv file for IPF1 and NML1/NML2 was renamed to features.tsv to comply with Seurat's Read10X requirements.
