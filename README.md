# Single-Cell-RNA-seq-Analysis-of-Idiopathic-Pulmonary-Fibrosis-IPF-vs.-Normal-Human-Lung

![Dataset](https://img.shields.io/badge/Dataset-GSE132771-0C447C?style=flat-square) ![Platform](https://img.shields.io/badge/Platform-10X%20Genomics-085041?style=flat-square) ![Tool](https://img.shields.io/badge/Tool-Seurat%20v4%20%7C%20R-3C3489?style=flat-square) ![Donors](https://img.shields.io/badge/Donors-6%20%7C%20up%20to%2018%2C000%20cells-633806?style=flat-square) ![Integration](https://img.shields.io/badge/Integration-CCA-444441?style=flat-square)

Idiopathic Pulmonary Fibrosis (IPF) is a progressive, fatal lung disease driven by excessive fibrotic remodeling of the lung parenchyma. This project performs an end-to-end scRNA-seq analysis comparing lung tissue from 3 IPF patients and 3 normal donors — resolving major lung cell populations and pinpointing the transcriptomic signature of pathological fibroblast activation in IPF.

---

## Data

| Group | Donors | GEO Accessions |
|-------|--------|----------------|
| Normal | 3 | `GSM3891621` · `GSM3891623` · `GSM3891625` |
| IPF | 3 | `GSM3891627` · `GSM3891629` · `GSM3891631` |

Raw data downloaded from NCBI GEO as 10X-format files (`matrix.mtx`, `barcodes.tsv`, `features.tsv`).

---

## Pipeline

**1 · Quality Control**
Samples filtered for 500–4,000 genes/cell, <20,000 UMIs, and <10% mitochondrial content. IPF1 thresholds widened (up to 6,000 genes, 40,000 UMIs) to retain high-complexity fibrotic cells. Each donor subsampled to a maximum of `3,000 cells`.

**2 · Normalization & Feature Selection**
Log-normalization applied per sample via `NormalizeData()`. `2,000` highly variable genes selected using variance-stabilizing transformation (VST).

**3 · Integration**
Canonical Correlation Analysis (CCA) via Seurat's anchor-based framework (`FindIntegrationAnchors` → `IntegrateData`) to correct donor-level batch effects while preserving biological variation.

**4 · Dimensionality Reduction & Clustering**
PCA on 50 components → top `20 PCs` selected by elbow plot → UMAP embedding → graph-based clustering at resolution `0.1`.

**5 · Cell Type Annotation**
8 major populations annotated using canonical lineage markers:

| Cell Type | Markers |
|-----------|---------|
| Monocyte-derived Macrophages | `S100A8` · `FCN1` |
| Endothelial | `EMCN` · `PECAM1` · `CLDN5` |
| Epithelial | `EPCAM` · `SFTPC` · `AGER` |
| NK Cells | `KLRD1` · `PRF1` |
| Fibroblasts | `COL1A2` · `DCN` · `PDGFRA` |
| Dendritic Cells | `CD1C` · `FCER1A` |
| Plasma B Cells | `IGHG1` · `IGKC` |
| Mast Cells | `TPSB2` |

**6 · Differential Expression — Fibroblasts (IPF vs. Normal)**
Wilcoxon rank-sum test via `FindMarkers()`, restricted to the fibroblast cluster. Thresholds: `|log2FC| > 0.5`, `adj. p < 0.05`, `min.pct = 0.25`. Results visualized on a volcano plot.

**7 · Pathway Enrichment**
Upregulated IPF fibroblast genes submitted to EnrichR. Databases: GO Biological Process 2023, KEGG 2021 Human. Immunoglobulin genes (`IGH`/`IGK`/`IGL`) removed prior to enrichment to reduce noise.

---

## Key Findings

> **CTHRC1** — top upregulated gene in IPF fibroblasts. A pathological myofibroblast marker virtually absent in normal lung, exclusively expressed in IPF on UMAP split, confirming disease-specific fibroblast activation.

GO enrichment of IPF-upregulated fibroblast genes revealed:
- Extracellular Matrix Organization
- Collagen Fibril Organization
- Supramolecular Fiber Organization
- Regulation of Apoptotic Process

These pathways collectively reflect the hallmark fibrotic program: excessive ECM deposition, collagen crosslinking, and apoptosis resistance in activated myofibroblasts.

---

## Output Files

| File | Description |
|------|-------------|
| `CONIPF_integrated.rds` | Integrated Seurat object (pre-annotation) |
| `CONIPF_annotated.rds` | Final annotated Seurat object |
| `cell_proportion_CON_vs_IPF.pdf` | Cell type composition bar chart |
| `featureplot_markers.pdf` | Lineage marker feature plots |
| `CTHRC1_featureplot.pdf` | CTHRC1 expression split by condition |
| `GO_enrichment.pdf` | GO Biological Process enrichment chart |

---

## Dependencies

```r
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(enrichR)
```

R version 4.x recommended. `set.seed(42)` used throughout for reproducibility.

---

## Notes

> `genes.tsv` renamed to `features.tsv` for IPF1, NML1, NML2 to comply with Seurat's `Read10X()` requirements.

> Subsampling to 3,000 cells/donor applied before integration to prevent dominant donors from skewing cluster structure.

> IPF1 QC thresholds widened intentionally — fibrotic tissue carries higher transcriptional complexity and cells would be incorrectly dropped by standard cutoffs.
