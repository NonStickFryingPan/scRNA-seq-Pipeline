<div align="center">

# 🧫 scRNA-Seq Pipeline
**Demultiplexing, alignment, and quality filtering of 10x Chromium single-cell RNA-seq data into a high-quality count matrix**

[![Galaxy](https://img.shields.io/badge/Workflow-Galaxy-1f6feb?style=flat-square&logo=galaxy&logoColor=white)](https://usegalaxy.org/)
[![STARsolo](https://img.shields.io/badge/Aligner-RNA%20STARsolo-16a34a?style=flat-square)](https://github.com/alexdobin/STAR)
[![DropletUtils](https://img.shields.io/badge/Filtering-DropletUtils-e11d48?style=flat-square)](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html)
[![10x Chromium](https://img.shields.io/badge/Platform-10x%20Chromium%20v3-f59e0b?style=flat-square)](https://www.10xgenomics.com/)
[![Genome](https://img.shields.io/badge/Reference-hg19%20GRCh37-7c3aed?style=flat-square)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/)

[Overview](#overview) · [Dataset](#dataset) · [Pipeline](#pipeline) · [Tools & Outputs](#tools--outputs) · [Key Results](#key-results) · [References](#references)

</div>

---

## Overview

This pipeline processes 10x Chromium v3 scRNA-seq FASTQ files into a filtered, high-quality gene × cell count matrix. It uses **RNA STARsolo** as a drop-in replacement for the Cell Ranger pipeline — significantly faster with equivalent output. Cell filtering is handled by **DropletUtils**, using both a Cell Ranger-equivalent default method and a custom EmptyDrops approach for comparison.

**Sample:** 1k PBMCs from a healthy donor (10x Genomics v3 chemistry)  
**Subsampled to:** ~300 cells for tutorial purposes

---

## Dataset

| File | Content |
|------|---------|
| `subset_pbmc_1k_v3_S1_L001_R1_001.fastq.gz` | Cell barcodes + UMI (Lane 1) |
| `subset_pbmc_1k_v3_S1_L001_R2_001.fastq.gz` | cDNA reads (Lane 1) |
| `subset_pbmc_1k_v3_S1_L002_R1_001.fastq.gz` | Cell barcodes + UMI (Lane 2) |
| `subset_pbmc_1k_v3_S1_L002_R2_001.fastq.gz` | cDNA reads (Lane 2) |
| `Homo_sapiens.GRCh37.75.gtf` | Gene annotation |
| `3M-february-2018.txt.gz` | 10x cell barcode whitelist (~3.7M barcodes) |

---

## Pipeline

```
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

## Tools & Outputs

### Stage 1 — Alignment & Demultiplexing

**Tool:** `RNA STARsolo v2.7.11a`

Simultaneously aligns cDNA reads to the reference genome and demultiplexes them by cell barcode. R1 reads are parsed for barcodes/UMIs; R2 reads carry the actual cDNA sequence. Each read is matched against the 10x whitelist to assign it to a known cell barcode.

| Parameter | Value |
|-----------|-------|
| Reference | hg19 GRCh37 (built-in) |
| Annotation | GRCh37.75 GTF |
| Chemistry | Chromium v3 |
| CB length | 16 bp |
| UMI length | 12 bp |
| UMI deduplication | CellRanger2-4 algorithm |
| UMI filtering | Remove N + homopolymers |
| Barcode matching | 1MM_multi (CellRanger 2 style) |
| Cell filtering | Disabled (handled by DropletUtils) |
| Strandedness | Forward (same strand as RNA) |

**Mapping quality (MultiQC):**

| Metric | Value |
|--------|-------|
| Uniquely mapped reads | 87.5% |
| Total barcodes detected | 5,200 |
| Reads assigned to unique gene features | 3,841,347 |

> 87.5% unique mapping is expected and acceptable for 10x datasets.

**Outputs:** BAM file, log, feature statistics summary, raw count matrix (bundled: `matrix.mtx` + `barcodes.tsv` + `genes.tsv`)

---

### Stage 2a — Cell Ranger Method (DefaultDrops)

**Tool:** `DropletUtils v1.10.0`

Emulates the Cell Ranger filtering approach. Retains barcodes above a UMI count threshold derived from the top `upper_quantile` of the expected number of cells, then keeps barcodes with at least `lower_proportion` × that threshold.

| Parameter | Value |
|-----------|-------|
| Method | DefaultDrops |
| Expected cells | 3,000 |
| Upper quantile | 0.99 |
| Lower proportion | 0.1 |

**Result: 272 high-quality cells**

---

### Stage 2b — Introspective Method (EmptyDrops)

**Tool:** `DropletUtils v1.10.0` — run twice

**Step 1 — Rank Barcodes:** Generates a barcode rank plot (log UMI count vs. log rank) to identify the knee and inflection points separating real cells from empty droplets.

| Threshold | Value | Interpretation |
|-----------|-------|----------------|
| Knee | 4,861 | Upper bound — ~100–200 high-RNA cells |
| Inflection | 260 | Lower bound — up to ~400 cells |

**Step 2 — EmptyDrops filtering:** Uses a statistical model to distinguish real cells from ambient RNA in empty droplets, with a custom lower-bound UMI threshold and FDR control.

| Parameter | Value |
|-----------|-------|
| Method | EmptyDrops |
| Lower-bound threshold | 200 |
| FDR threshold | 0.01 (1% false positives) |

**Result: 279 high-quality cells**

> EmptyDrops recovers 7 additional cells vs. DefaultDrops. On full-size datasets this difference can meaningfully reduce noise in downstream clustering.

---

## Key Results

| Metric | Value |
|--------|-------|
| Raw barcodes detected (STARsolo) | 5,200 |
| Cells after DefaultDrops filtering | 272 |
| Cells after EmptyDrops filtering | 279 |
| Unique mapping rate | 87.5% |
| Chemistry | 10x Chromium v3 (confirmed from R1 read length = 28 bp) |
| Final matrix dimensions | ~33,000 genes × ~279 cells |

---

## Glossary

| Term | Definition |
|------|------------|
| **Cell Barcode (CB)** | 16 bp sequence uniquely identifying each cell droplet |
| **UMI** | Unique Molecular Identifier — tags individual molecules to remove PCR duplicates |
| **Whitelist** | Known set of valid 10x cell barcodes used to correct sequencing errors in CB reads |
| **Empty droplet** | A droplet containing no cell, only ambient RNA — filtered out by DropletUtils |
| **Knee point** | UMI rank threshold above which barcodes likely represent real cells |
| **Count matrix** | Sparse genes × cells matrix of UMI counts — primary input for all downstream scRNA-seq analysis |

---

## References

- Dobin et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*
- Lun et al. (2019). EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. *Genome Biology*
- Tekman et al. (2020). A single-cell RNA-seq analysis pipeline. *GigaScience*
