<div align="center">

# 🧫 scRNA-Seq Pipeline

**Demultiplexing, alignment, and quality filtering of 10x Chromium single-cell RNA-seq data into a high-quality count matrix**

[![Galaxy](https://img.shields.io/badge/Workflow-Galaxy-1f6feb?style=flat-square&logo=galaxy&logoColor=white)](https://usegalaxy.org/)
[![STARsolo](https://img.shields.io/badge/Aligner-RNA%20STARsolo-16a34a?style=flat-square)](https://github.com/alexdobin/STAR)
[![DropletUtils](https://img.shields.io/badge/Filtering-DropletUtils-e11d48?style=flat-square)](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html)
[![10x Chromium](https://img.shields.io/badge/Platform-10x%20Chromium%20v3-f59e0b?style=flat-square)](https://www.10xgenomics.com/)
[![Genome](https://img.shields.io/badge/Reference-hg19%20GRCh37-7c3aed?style=flat-square)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/)

[Overview](#overview) • [Dataset](#dataset) • [Pipeline](#pipeline) • [Upstream Analysis](#upstream-analysis-galaxy) • [Downstream Analysis](#downstream-analysis-python--scanpy) • [AnnData](#the-anndata-structure) • [References](#references)

</div>

---

## Overview

This pipeline processes 10x Chromium v3 scRNA-seq FASTQ files into a filtered, high-quality gene × cell count matrix. It uses **RNA STARsolo** as a drop-in replacement for the Cell Ranger pipeline — significantly faster with equivalent output. Cell filtering is handled by **DropletUtils**.

- **Sample:** 1k PBMCs from a healthy donor (10x Genomics v3 chemistry)
- **Subsampled to:** ~300 cells for tutorial purposes

---

## Dataset

| File | Content |
| :--- | :--- |
| `subset_pbmc_1k_v3_S1_L001_R1_001.fastq.gz` | Cell barcodes + UMI (Lane 1) |
| `subset_pbmc_1k_v3_S1_L001_R2_001.fastq.gz` | cDNA reads (Lane 1) |
| `subset_pbmc_1k_v3_S1_L002_R1_001.fastq.gz` | Cell barcodes + UMI (Lane 2) |
| `subset_pbmc_1k_v3_S1_L002_R2_001.fastq.gz` | cDNA reads (Lane 2) |
| `Homo_sapiens.GRCh37.75.gtf` | Gene annotation |
| `3M-february-2018.txt.gz` | 10x cell barcode whitelist (~3.7M barcodes) |

---

## Pipeline

```text
R1 (barcodes) + R2 (cDNA) FASTQs — 2 lanes
            │
            ▼
[1] Spliced Alignment + Demultiplexing (RNA STARsolo)
      Uses barcode whitelist to assign reads to cells
            │
            ├── BAM (alignments)
            ├── Log + Feature Statistics
            └── Raw count matrix (genes × ~5200 barcodes)
                        │
                        ▼
            [2a] Cell Ranger method (DropletUtils DefaultDrops)
                  Expected cells: 3000 → yields ~272 cells
                        │
            [2b] Introspective method (DropletUtils EmptyDrops)
                  Barcode rank plot → knee/inflection thresholds
                  Custom lower bound: 200, FDR: 0.01 → yields ~279 cells
                        │
                        ▼
            Filtered count matrix (genes × high-quality cells)
                  Ready for downstream clustering & analysis
```

---

## Upstream Analysis (Galaxy)

### Stage 1: Alignment & Demultiplexing
**Tool:** `RNA STARsolo v2.7.11a`
Simultaneously aligns cDNA reads to the reference genome and demultiplexes them by cell barcode.

| Parameter | Value |
| :--- | :--- |
| Reference | hg19 GRCh37 |
| Chemistry | Chromium v3 |
| CB/UMI length | 16 bp / 12 bp |
| UMI Filtering | Remove N + homopolymers |
| Mapping Rate | 87.5% Unique |

### Stage 2: Cell Filtering
**Tool:** `DropletUtils v1.10.0`
We utilized the **EmptyDrops** method to distinguish real cells from empty droplets containing only ambient RNA.

| Threshold | Value | Interpretation |
| :--- | :--- | :--- |
| Knee Point | 4,861 | Upper bound: high-RNA cells |
| Inflection | 260 | Lower bound: potential cells |
| **Final Yield** | **279 cells** | High-quality filtered cells |

---

## Downstream Analysis (Python / Scanpy)

The filtered Galaxy matrix was imported into Python as an **AnnData** object to leverage its multi-layered data structure.

### Workflow
1. **Load & Format:** Assigned Human Gene Symbols to the matrix columns for biological interpretability.
2. **QC Filtering:** Removed cells with fewer than 100 genes to preserve data in this smaller dataset.
3. **Normalization:** Scaled total counts per cell and applied log transformation (`log1p`).
4. **Dimensionality Reduction:** Computed PCA to denoise data and UMAP to visualize cell relationships in 2D space.
5. **Annotation:** Clusters were identified via the Leiden algorithm (res 0.5) and annotated using **CellTypist** automated labeling alongside manual marker gene verification.

---

## The AnnData Structure

This project served as an exploration of the AnnData format, which keeps metadata and raw data synchronized:

* **`.X`**: The core gene expression count matrix.
* **`.obs`**: Cell metadata (cluster labels, QC metrics, sample IDs).
* **`.var`**: Gene metadata (symbols, highly variable flags).
* **`.obsm`**: Multi-dimensional coordinates (PCA, UMAP).

---

## Results & Key Findings

### Final Analysis Summary
| Metric | Value |
| :--- | :--- |
| Initial Cell Count | 252 |
| Final Cell Count (Post-QC) | 248 |
| Total Genes in Dataset | 2,392 |
| Median Genes per Cell | ~145 |
| Identified Populations | Lymphocytes (T-Cells), Myeloid (Monocytes) |
| Principal Components (PCs) | 50 |
| Software Stack | Galaxy, Scanpy, CellTypist, AnnData |

### Biological Insights
* **Population Diversity:** UMAP visualization revealed two distinct biological islands. The majority of the sample consists of **Lymphocytes**, with a smaller distinct cluster of **Myeloid cells**.
* **Marker Verification:** Cluster 1 was confirmed as Myeloid through high expression of `CYBB` and `SAT1`. Cluster 0 was confirmed as Lymphocyte-rich due to high ribosomal gene expression (`RPL` and `RPS`).
* **Mitochondrial Content:** The analysis showed 0% mitochondrial expression, indicating these were filtered upstream in Galaxy or not captured during sequencing.
* **Portability:** The entire analysis is preserved in a single `.h5ad` file, maintaining the relationship between raw counts and 2D visualization coordinates.

> [!NOTE]
> **On Scale:** While the reference tutorial used 8,785 cells, this exploration focused on a subsampled 252-cell set, requiring more permissive QC thresholds to maintain statistical power.

---

## Glossary

* **Cell Barcode:** A 16 bp sequence uniquely identifying each cell droplet.
* **UMI:** Unique Molecular Identifier used to tag individual molecules and remove PCR duplicates.
* **Empty Droplet:** A droplet containing no cell, only ambient RNA; filtered out by DropletUtils.
* **Knee Point:** A UMI rank threshold above which barcodes likely represent real cells.

---

## References

* Dobin et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*
* Lun et al. (2019). EmptyDrops: distinguishing cells from empty droplets. *Genome Biology*
* Wolf et al. (2018). SCANPY: Large-scale single-cell gene expression data analysis. *Genome Biology*
