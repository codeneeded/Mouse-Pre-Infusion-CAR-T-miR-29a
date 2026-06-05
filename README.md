# Murine Pre-Infusion CAR-T scRNA-seq — Effect of miR-29a

> Single-cell RNA-seq of murine pre-infusion CAR-T cell products engineered to overexpress miR-29a, characterizing how miR-29a reshapes the transcriptional states, metabolic programs, and subset composition of the cell product prior to infusion. Code repository for a manuscript in preparation.

---

## Overview

This repository contains the R analysis pipeline for a single-cell RNA-sequencing (scRNA-seq) study of murine chimeric antigen receptor T cell (CAR-T) products engineered to overexpress **miR-29a**. The aim is to define, at single-cell resolution, how miR-29a modulates the differentiation state, effector and exhaustion programs, metabolic state, and proliferative behaviour of the CAR-T product at the pre-infusion stage.

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
- **Primary question:** How does miR-29a overexpression alter the composition, transcriptional state, and metabolic programs of the CAR-T product relative to scramble and empty-vector controls?

---

## Repository Structure

```
Mouse-Pre-Infusion-CAR-T-miR-29a/
│
├── Scripts/                           # Seurat v5 analysis pipeline (R)
├── Resources/                         # Curated input resources
│                                      # (module gene lists; miR-29a
│                                      #  TargetScan target list)
├── QC/                                # Quality-control outputs
│                                      # (filtering metrics, doublet removal,
│                                      #  cell-cycle scoring)
├── Integration/                       # Batch correction and integration
│                                      # (method comparison, composition checks)
├── Annotation/                        # Cluster annotation outputs
│                                      # (UMAPs, marker dot plots, heatmaps,
│                                      #  lineage assignment, cluster labels)
├── Differential_Expression/           # Per-lineage and per-cluster DE
│                                      # (pseudobulk DESeq2 and MAST in
│                                      #  parallel, with method overlap)
├── Pathway_Analysis_EnrichR/          # Transcription factor, pathway, and
│                                      # miRNA-target enrichment (EnrichR)
├── Module_Scores/                     # Per-cell module scores for metabolic
│                                      # programs, FOXO axis, T cell state,
│                                      # and miR-29a target repression
├── miR29a_Target_Enrichment/          # Direct test of miR-29a target
│                                      # repression against the conserved
│                                      # TargetScan list
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
- Cell-cycle phase scoring with mouse-converted Tirosh gene sets

### 3. Integration (`Integration/`)
- Log-normalization, variable-feature selection, scaling, and PCA
- Cell-cycle regression of the S − G2M difference to preserve cycling-vs-resting structure while removing phase fragmentation
- Batch correction across replicates using **Harmony** (primary), with **FastMNN** retained as an alternative integration for robustness checks
- Integration assessed for over-correction by confirming condition structure is preserved across clusters

### 4. Cell Type Annotation (`Annotation/`)
- Graph-based clustering on the integrated embedding; resolution selected using clustering-stability diagnostics
- Cluster identities assigned from canonical murine T-cell markers (naïve/memory, effector/cytotoxic, exhaustion, regulatory, proliferation, lineage transcription factors)
- Per-cell lineage assignment (CD4 / CD8 / γδ) from lineage markers
- Marker dot plots, FeaturePlots, and z-scored expression heatmaps per cluster

### 5. Differential Expression (`Differential_Expression/`)
- **Pseudobulk DESeq2** (primary): aggregation per condition × replicate within each lineage compartment, design modelling replicate as a covariate, log fold-change shrinkage. Calibrated at low replicate count.
- **MAST** (parallel exploratory): single-cell-level differential expression with sequencing depth as a latent variable, run on the same compartments and contrasts for comparability.
- Two compartment scopes: **per lineage** (CD4, CD8, and minor non-T compartments analysed separately) and **per cluster** (every cluster in the integrated atlas).
- Contrasts: miR-29a vs scramble (primary, specificity); miR-29a vs empty vector (secondary); empty vector vs scramble (QC).
- Method comparison output summarises gene overlap between pseudobulk and MAST per compartment × contrast.

### 6. Pathway Analysis (`Pathway_Analysis_EnrichR/`)
- Enrichment of DE gene lists via **EnrichR** across three database groups:
  - **Transcription factor regulators** (TRRUST, ChEA, JASPAR PWMs)
  - **Pathways and ontologies** (KEGG / WikiPathways Mouse, GO Biological Process, MSigDB Hallmark, Reactome, BioPlanet, Panther)
  - **miRNA-target databases** (miRTarBase, TargetScan)
- miRNA enrichment is run separately on up- and down-regulated gene lists; significant miR-29 hits are flagged in dedicated output files.

### 7. Module Scores (`Module_Scores/`)
- Per-cell module scores via `AddModuleScore` for curated gene sets from `Resources/Module_Gene_Lists.xlsx`:
  - **Metabolism modules** — OXPHOS / ETC, TCA, FAO, glycolysis, mitochondrial biogenesis, mitochondrial stress / ROS, mitochondrial dynamics, mTOR / MYC anabolic state, proliferative metabolic state
  - **FOXO axis modules** — core FOXO TFs, memory/quiescence-linked targets, stress/autophagy targets, cell-cycle restraint
  - **T cell state modules** — cell cycle / G2M, exhaustion / inhibitory receptors, stem-like / memory
  - **miR-29a target module** — derived from the conserved TargetScan list (top-N strongest predicted targets)
- The actual gene set used for each module (after filtering to genes present in the data) is exported as a long-format table for reproducibility.
- Visualizations include UMAP feature plots, per-cluster violin plots, and condition-stratified violins within each cluster — organized into per-view subfolders to support cross-module comparison.
- One-sided Wilcoxon tests compare the miR-29a target module score between miR-29a and each control, overall and per cluster.

### 8. miR-29a Target Enrichment (`miR29a_Target_Enrichment/`)
- Direct test that miR-29a is functionally repressing its annotated target set, against the **conserved mouse miR-29-3p TargetScan list** (`Resources/miR29a_targetscan_conserved.csv`).
- Two complementary analyses per DE contrast / compartment / method:
  - **Spearman correlation** between TargetScan cumulative weighted context++ score and experimental log2 fold change. A positive correlation indicates that genes predicted to be the strongest targets are the most strongly repressed in the data.
  - **Hypergeometric enrichment** of significantly down-regulated genes against the top-N TargetScan targets and the full conserved list, with up-regulated genes tested as a negative control.
- Outputs per-contrast scatter plots (TargetScan score vs experimental log2 fold change, points coloured by differential-expression class — sig-down, sig-up, n.s. — with the strongest sig-down predicted targets labelled), a summary table with BH-adjusted p-values, and a headline summary plot for the primary lineage-level contrasts.

---

## Preliminary Findings

Initial analysis indicates that **miR-29a is functionally active** in the CAR-T pre-infusion product: TargetScan-predicted miR-29a targets are coherently repressed in both CD4 and CD8 compartments versus scramble control (Spearman ρ = 0.18–0.20, *p*<sub>BH</sub> < 10⁻⁶; 5–7× hypergeometric enrichment of significantly down-regulated genes in the top-200 predicted targets). The empty-vector vs scramble comparison shows no correlation or enrichment, confirming that the scramble construct controls appropriately for backbone effects. Canonical miR-29 targets including *Tet2* and *Tet3* are among the most strongly repressed genes.

The miR-29a transcriptome also shows the expected **miRNA-repressor direction asymmetry** — significantly more down-regulated than up-regulated genes (CD4: 1.9×, CD8: 2.4×) — while the scramble vs empty-vector control is balanced, as expected for a null.

### Summary figures (miR29a_Target_Enrichment/Plots/Summary_<method>/)

Four manuscript-grade summary figures are generated per DE method (pseudobulk and MAST). The pseudobulk versions are the inferentially-honest manuscript figures; MAST is shown alongside as a sensitivity analysis and gives concordant patterns at larger sample sizes.

1. **`1_spearman_forest.png`** — Spearman ρ ± 95% CI for each of 18 compartments × 3 contrasts. Every miR-29a contrast (navy / mid-blue) sits to the right of ρ = 0; every EV-vs-Scr control (grey) has its CI crossing zero. Read this as: "miR-29a coherently represses TargetScan-predicted targets across every compartment we can power-test; the scramble control behaves as a clean null."

2. **`2_hyper_enrichment_dotplot.png`** — Hypergeometric enrichment of sig-down genes in the top-200 strongest predicted miR-29-3p targets. Larger and warmer = stronger enrichment. The two miR-29a columns are populated with warm orange-to-red dots; the EV-vs-Scr column is uniformly pale and small. Smaller compartments (Non-T, Innate-like, the two review clusters) show the highest odds ratios — these are real, but driven by small denominators; lineage-level CD4 and CD8 effect sizes (5–7×) are the manuscript numbers.

3. **`3_direction_asymmetry.png`** — Log-log scatter of sig-down vs sig-up gene counts per (compartment, contrast). The lineage-level points (large circles) for miR-29a contrasts sit clearly below the y=x diagonal — more genes go down than up — while the EV-vs-Scr control points sit on or slightly above the diagonal. This is the simplest, target-list-independent demonstration that miR-29a behaves as a transcriptome-wide repressor.

4. **`4_canonical_targets_heatmap.png`** — log2FC of canonical miR-29 targets and strong data-driven hits across all 18 compartments, primary contrast only (miR-29a vs Scr). The DNA-methylation axis (*Tet2, Tet3, Dnmt3a, Dnmt3b, Tdg*) is coordinately repressed across nearly every compartment. *Eomes* is the single most consistently repressed T-cell-relevant gene. *Bach2* and *Foxo3* are state-dependent — repressed in major lineages but *upregulated* in the stem-like CD4 cluster (5) and the intermediate-activated cluster (12), suggesting either secondary network effects or cluster-level cell-state heterogeneity that would benefit from per-lineage re-clustering.

See [`RESULTS_EXPLANATION.txt`](RESULTS_EXPLANATION.txt) for a full walkthrough with effect sizes, robustness checks across methods, caveats, and a list of analyses still pending before manuscript submission.

> ⚠️ Cluster composition shifts (relative abundance of CAR-T subsets between conditions) are visible in the integration composition plots but have not yet been formally tested with per-replicate differential abundance methods. These observations are not reported as findings in the manuscript until that validation is complete. Cluster 12 in particular shows several flags (large apparent enrichment in miR-29a, *Mycn* and *Bach2* anomalously upregulated in the target heatmap) and is the highest priority for per-replicate validation.

---

## Scientific Questions

1. **Subset composition** — Does miR-29a overexpression alter the proportions of CAR-T subsets (naïve/stem-like, effector, exhausted, regulatory, proliferating)?
2. **Differentiation state** — Does miR-29a bias the product toward a less-differentiated, less-exhausted phenotype?
3. **Target repression** — Are canonical miR-29a targets (*Tbx21*, *Eomes*, *Dnmt3a* / *Dnmt3b*) and the broader conserved TargetScan target set repressed at single-cell resolution?
4. **Metabolic and FOXO programs** — Does miR-29a remodel metabolic state (OXPHOS, glycolysis, FAO, mitochondrial dynamics) or FOXO-axis programs (memory / quiescence versus effector)?
5. **Proliferation** — Does miR-29a change the proliferative composition of the pre-infusion product?
6. **Control specificity** — Are miR-29a–associated effects specific relative to both scramble and empty-vector controls?

---

## Dependencies

All scripts are written in **R**. Key packages:

| Package | Purpose |
|---|---|
| `Seurat` (v5) | Clustering, single-cell DE, visualization, module scoring |
| `harmony` | Batch integration (primary) |
| `SeuratWrappers` / `batchelor` | FastMNN integration |
| `scDblFinder` | Doublet detection |
| `DESeq2` / `ashr` | Pseudobulk differential expression with LFC shrinkage |
| `MAST` | Single-cell differential expression |
| `enrichR` | Pathway, TF, and miRNA-target enrichment |
| `qs2` | Fast serialization of Seurat objects |
| `SeuratExtend` / `scCustomize` | Single-cell visualization |
| `readxl` / `openxlsx` | Curated gene-list and enrichment workbook I/O |
| `dplyr` / `tidyr` / `ggplot2` / `patchwork` | Data wrangling and visualization |
| `pheatmap` / `ggrepel` | Heatmaps and labelled plots |

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
