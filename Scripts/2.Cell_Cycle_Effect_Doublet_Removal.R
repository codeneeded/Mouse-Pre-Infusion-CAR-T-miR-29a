###############################################################################
## 02_cellcycle_doublets_mouse_CARTmiR29a.R
## Murine pre-infusion CAR-T (miR-29a) — single-cell, GEX only.
## Runs on the QC-filtered object from script 01.
##
##   1) preliminary (unintegrated) embedding for QC visualization
##   2) cell-cycle scoring  (MOUSE gene lists; SCORED + VISUALIZED, not regressed)
##   3) per-sample doublet detection with scDblFinder -> remove doublets
##   4) save the cell-cycle-scored, doublet-clean object
##
## ADAPTED FROM the human CITE-seq version — note the differences:
##   - Azimuth `pbmcref` REMOVED: it is a HUMAN reference and cannot annotate
##     mouse data. Cell-type annotation belongs in a later step with mouse
##     markers / a mouse reference.
##   - cc.genes.updated.2019 are HUMAN symbols; converted to mouse casing below
##     (Pcna, Mki67, ...). The match count is printed so you can confirm it.
##   - Cell cycle is SCORED here but NOT regressed. Decide on regression in the
##     integration script after seeing whether Phase drives clustering.
##   - This embedding is UNINTEGRATED and exists only for QC visualization;
##     batch integration across replicates happens in script 03.
###############################################################################

# ---- Libraries ----
library(Seurat)
library(scDblFinder)
library(scCustomize)
library(tidyverse)     # dplyr + ggplot2 + stringr
library(patchwork)
library(qs2)

set.seed(123)          # scDblFinder + UMAP are stochastic

# ============================ Paths ==========================================
project_root   <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a"
saved_dir      <- file.path(project_root, "saved_R_data")
cell_cycle_dir <- file.path(project_root, "QC", "Cell_Cycle")
doublet_dir    <- file.path(project_root, "QC", "Doublets")
for (d in c(cell_cycle_dir, doublet_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Load filtered object from script 01 ----
filtered_seurat <- qs_read(file.path(saved_dir,
                                     "Mouse_CARTmiR29a_filtered_seurat.qs2"))
DefaultAssay(filtered_seurat) <- "RNA"
seurat_phase <- filtered_seurat

# ===================== 1) Preliminary embedding (QC view) ====================
# Unintegrated — only for visualizing Phase + doublets. Integration is script 03.
seurat_phase <- NormalizeData(seurat_phase)
seurat_phase <- FindVariableFeatures(seurat_phase, selection.method = "vst",
                                     nfeatures = 2000)
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase, npcs = 30, verbose = FALSE)
seurat_phase <- FindNeighbors(seurat_phase, dims = 1:30, reduction = "pca")
seurat_phase <- FindClusters(seurat_phase, resolution = 0.8,
                             cluster.name = "unintegrated_clusters")
seurat_phase <- RunUMAP(seurat_phase, dims = 1:30, reduction = "pca",
                        reduction.name = "umap.unintegrated")

# ===================== 2) Cell-cycle scoring (MOUSE) =========================
# Human Tirosh lists -> mouse casing. Title-casing maps the large majority of
# these symbols (PCNA->Pcna, MKI67->Mki67); keep only those present in the data.
s.genes.mm   <- intersect(str_to_title(cc.genes.updated.2019$s.genes),
                          rownames(seurat_phase))
g2m.genes.mm <- intersect(str_to_title(cc.genes.updated.2019$g2m.genes),
                          rownames(seurat_phase))
message("Cell-cycle genes matched -> S: ", length(s.genes.mm), "/",
        length(cc.genes.updated.2019$s.genes), " | G2M: ", length(g2m.genes.mm),
        "/", length(cc.genes.updated.2019$g2m.genes))

seurat_phase <- CellCycleScoring(seurat_phase,
                                 s.features   = s.genes.mm,
                                 g2m.features = g2m.genes.mm,
                                 set.ident    = FALSE)
seurat_phase$Phase <- factor(seurat_phase$Phase, levels = c("G1", "S", "G2M"))

# ---- visualize cell cycle ----
ridge_markers <- intersect(c("Pcna", "Top2a", "Mcm6", "Mki67"),
                           rownames(seurat_phase))
png(file.path(cell_cycle_dir, "RidgePlot_CellCycleMarkers.png"),
    width = 1800, height = 1200)
RidgePlot(seurat_phase, features = ridge_markers, ncol = 2)
dev.off()

png(file.path(cell_cycle_dir, "CellCycle_Scores_bySample.png"),
    width = 1800, height = 1200)
VlnPlot(seurat_phase, features = c("S.Score", "G2M.Score"),
        group.by = "orig.ident", pt.size = 0.1)
dev.off()

png(file.path(cell_cycle_dir, "CellCycle_Scores_byCondition.png"),
    width = 1800, height = 1200)
VlnPlot(seurat_phase, features = c("S.Score", "G2M.Score"),
        group.by = "condition", pt.size = 0.1)
dev.off()

# PCA + UMAP colored by Phase — the key question: does proliferation drive
# the structure? If clusters separate by Phase, consider regressing in script 03.
ggsave(file.path(cell_cycle_dir, "PCA_by_Phase.png"),
       DimPlot_scCustom(seurat_phase, reduction = "pca", group.by = "Phase"),
       width = 11, height = 8, dpi = 300)
ggsave(file.path(cell_cycle_dir, "UMAP_by_Phase.png"),
       DimPlot_scCustom(seurat_phase, reduction = "umap.unintegrated",
                        group.by = "Phase"),
       width = 11, height = 8, dpi = 300)
ggsave(file.path(cell_cycle_dir, "UMAP_Phase_byCondition.png"),
       DimPlot_scCustom(seurat_phase, reduction = "umap.unintegrated",
                        group.by = "Phase", split.by = "condition"),
       width = 16, height = 6, dpi = 300)

# proportion of cells per phase, per sample
phase_tab <- seurat_phase@meta.data %>%
  dplyr::count(orig.ident, Phase) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::mutate(pct = round(100 * n / sum(n), 1)) %>%
  dplyr::ungroup()
write.csv(phase_tab,
          file.path(cell_cycle_dir, "CellCycle_phase_proportions.csv"),
          row.names = FALSE)
print(phase_tab)

# ===================== 3) Doublet detection (scDblFinder) ====================
# Per-sample in a single call via `samples=`. Operates on raw counts; each
# replicate is processed independently (correct — they're separate captures).
counts_mat <- GetAssayData(seurat_phase, assay = "RNA", layer = "counts")
sce <- scDblFinder(counts_mat, samples = seurat_phase$orig.ident)

# scDblFinder preserves column order, so assign straight back
seurat_phase$scDblFinder.score <- sce$scDblFinder.score
seurat_phase$scDblFinder.class <- factor(sce$scDblFinder.class,
                                         levels = c("singlet", "doublet"))

# doublet rate per sample
dbl_tab <- seurat_phase@meta.data %>%
  dplyr::count(orig.ident, scDblFinder.class) %>%
  tidyr::pivot_wider(names_from = scDblFinder.class, values_from = n,
                     values_fill = 0) %>%
  dplyr::mutate(pct_doublet = round(100 * doublet / (doublet + singlet), 1))
write.csv(dbl_tab, file.path(doublet_dir, "Doublet_rates_per_sample.csv"),
          row.names = FALSE)
print(dbl_tab)

# visualize doublets on the unintegrated UMAP (singlet blue, doublet red)
ggsave(file.path(doublet_dir, "UMAP_doublet_class.png"),
       DimPlot_scCustom(seurat_phase, reduction = "umap.unintegrated",
                        group.by = "scDblFinder.class",
                        colors_use = c("#1f78b4", "#e31a1c")),
       width = 11, height = 8, dpi = 300)
ggsave(file.path(doublet_dir, "UMAP_doublet_class_bySample.png"),
       DimPlot_scCustom(seurat_phase, reduction = "umap.unintegrated",
                        group.by = "scDblFinder.class", split.by = "orig.ident",
                        colors_use = c("#1f78b4", "#e31a1c"), num_columns = 3),
       width = 16, height = 10, dpi = 300)

# ===================== 4) Remove doublets + save =============================
seurat_clean <- subset(seurat_phase, subset = scDblFinder.class == "singlet")
message("Cells before doublet removal: ", ncol(seurat_phase),
        " | after: ", ncol(seurat_clean))

qs_save(seurat_clean,
        file = file.path(saved_dir,
                         "Mouse_CARTmiR29a_CellCycle_DoubletClean.qs2"))
message("Saved -> Mouse_CARTmiR29a_CellCycle_DoubletClean.qs2")
###############################################################################