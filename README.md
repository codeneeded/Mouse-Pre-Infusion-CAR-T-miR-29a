# Murine Pre-Infusion CAR-T scRNA-seq — Effect of miR-29a

> Single-cell RNA-seq of murine pre-infusion CAR-T cell products engineered to overexpress miR-29a, characterizing how miR-29a reshapes the transcriptional states and subset composition of the cell product prior to infusion. Code repository for a manuscript in preparation.

---

## Overview

This repository contains the R analysis pipeline for a single-cell RNA-sequencing (scRNA-seq) study of murine chimeric antigen receptor T cell (CAR-T) products engineered to overexpress **miR-29a**. The aim is to define, at single-cell resolution, how miR-29a modulates the differentiation state, effector and exhaustion programs, and proliferative behaviour of the CAR-T product at the pre-infusion stage.

miR-29a is a microRNA that represses the T-box transcription factors T-bet (*Tbx21*) and Eomes (*Eomes*) and the de novo DNA methyltransferases *Dnmt3a* / *Dnmt3b*, and has been implicated in restraining terminal effector differentiation and exhaustion in T cells. Profiling the pre-infusion product resolves how miR-29a engineering alters the starting cell state that is subsequently infused.

**Manuscript in preparation.** No publication data are available yet.

---

## Experimental Design

Murine CAR-T products were generated under three conditions and profiled by scRNA-seq across two independent biological replicates (six libraries total).

| Condition | Description |
|---|---|
| **miR-29a** | CAR-T engineered to overexpress miR-29a |
| **Scr** | Scramble control — non-targeting miRNA, matched construct backbone |
| **EV** | Empty-vector control |

- **Replicates:** Two independent biological replicates per condition (separate production runs)
- **Model:** Murine CAR-T, pre-infusion product
- **Primary question:** How does miR-29a overexpression alter the composition and transcriptional state of the CAR-T product relative to scramble and empty-vector controls?

---

## Repository Structure

```
Mouse-Pre-Infusion-CAR-T-miR-29a/
│
├── Scripts/                           # Seurat v5 analysis pipeline (R)
├── QC/                                # Quality-control outputs
│                                      # (filtering metrics, doublet removal,
│                                      #  cell-cycle scoring)
├── Integration/                       # Batch correction and integration
│                                      # (method comparison, composition checks)
├── Annotation/                        # Cluster annotation outputs
│                                      # (UMAPs, marker dot plots, heatmaps,
│                                      #  lineage assignment, cluster labels)
├── Differential_Abundance/            # Subset proportion testing across
│                                      # conditions
├── Differential_Expression/           # Pseudobulk differential expression
│                                      # across conditions
├── saved_R_data/                      # Serialized Seurat objects (.qs2)
│                                      # — not tracked (see .gitignore)
│
├── Mouse-Pre-Infusion-CAR-T-miR-29a.Rproj
├── .gitignore
└── README.md
```

---

## Analysis Pipeline

### 1. Alignment
- Reads aligned and counted with **Cell Ranger** (`cellranger count`) against the mouse reference genome (GRCm39)
- One gene-expression library per condition per replicate

### 2. Quality Control (`QC/`)
- Per-cell filtering on minimum UMI count, minimum genes detected, gene complexity, and maximum mitochondrial and haemoglobin RNA fractions
- Doublet detection and removal per library
- Cell-cycle phase scoring

### 3. Integration (`Integration/`)
- Log-normalization, variable-feature selection, scaling, and PCA
- Batch correction across replicates using **Harmony** (primary), with **FastMNN** evaluated in parallel
- Integration assessed for over-correction by confirming condition structure is preserved across clusters

### 4. Cell Type Annotation (`Annotation/`)
- Graph-based clustering on the integrated embedding; resolution selected using clustering-stability diagnostics
- Cluster identities assigned from canonical murine T-cell markers (naïve/memory, effector/cytotoxic, exhaustion, regulatory, proliferation, lineage transcription factors)
- Per-cell lineage assignment (CD4 / CD8 / γδ) from lineage markers
- Marker dot plots and z-scored expression heatmaps per cluster

### 5. Differential Abundance (`Differential_Abundance/`)
- Subset proportion testing across conditions at the replicate level using **propeller** (speckle)
- Per-replicate proportions reported alongside test statistics

### 6. Differential Expression (`Differential_Expression/`)
- **Pseudobulk** aggregation per condition × replicate within each lineage compartment
- Differential expression with **DESeq2**, modelling replicate as a covariate
- Primary contrast: miR-29a vs scramble; secondary: miR-29a vs empty vector; confident hits called by concordance across both controls
- Focused testing of canonical miR-29a target genes

### 7. Robustness
- Key abundance results reproduced on an alternative integration (FastMNN) to confirm conclusions are not dependent on integration method

---

## Scientific Questions

1. **Subset composition** — Does miR-29a overexpression alter the proportions of CAR-T subsets (naïve/stem-like, effector, exhausted, regulatory, proliferating)?
2. **Differentiation state** — Does miR-29a bias the product toward a less-differentiated, less-exhausted phenotype?
3. **Target repression** — Are canonical miR-29a targets (*Tbx21*, *Eomes*, *Dnmt3a* / *Dnmt3b*) and downstream effector programs repressed at single-cell resolution?
4. **Proliferation** — Does miR-29a change the proliferative composition of the pre-infusion product?
5. **Control specificity** — Are miR-29a–associated effects specific relative to both scramble and empty-vector controls?

---

## Dependencies

All scripts are written in **R**. Key packages:

| Package | Purpose |
|---|---|
| `Seurat` (v5) | Clustering, differential expression, visualization |
| `harmony` | Batch integration (primary) |
| `SeuratWrappers` / `batchelor` | FastMNN integration |
| `scDblFinder` | Doublet detection |
| `DESeq2` | Pseudobulk differential expression |
| `speckle` | Differential abundance (propeller) |
| `qs2` | Fast serialization of Seurat objects |
| `SeuratExtend` / `scCustomize` | Single-cell visualization |
| `dplyr` / `ggplot2` / `patchwork` | Data wrangling and visualization |
| `pheatmap` / `ggrepel` | Heatmaps and labeled plots |

---

## Status

**Manuscript in preparation.** This repository is under active development. Code and figures will be updated as the project progresses.

> ⚠️ Raw sequencing data and processed Seurat objects are not included in this repository. Raw FASTQ files and processed count matrices will be deposited at **NCBI GEO** upon publication (accession: _TBD_).

---

## Citation

A full citation will be added upon publication.

---

## Contact

For questions, please open an issue in this repository.
