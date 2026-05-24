###############################################################################
## qc_pipeline_mouse_CARTmiR29a.R
## Murine pre-infusion CAR-T (miR-29a) — single-cell, GEX ONLY (no ADT/CITE).
##
## End-to-end basic QC:
##   1) read CellRanger `count` outputs (6 samples) -> per-sample Seurat objects
##   2) merge + attach mouse QC metadata -> save merged checkpoint
##   3) Pre-QC diagnostic plots
##   4) filter on (tunable) thresholds
##   5) Post-QC plots -> save filtered object
##
## Thresholds below are set from inspection of the Pre-QC plots for THIS run
## (clean, deeply-sequenced product: main UMI mass ~25k, ~5-6k genes/cell,
##  mito ~2%). Doublets are handled separately in script 03.
###############################################################################

# ---- Libraries ----
library(Seurat)
library(tidyverse)
library(Matrix)
library(patchwork)
library(scCustomize)   # QC_Plots_* helpers (install if missing)
library(qs2)           # fast save/load

# ============================ Paths ==========================================
# CellRanger output: one subfolder per sample, each containing /outs/
in.path <- "/media/akshay-iyer/Elements/data_from_hpc/Cell_Ranger_Out/"

# Project root + output folders (created if missing)
project_root <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a"
out.path     <- file.path(project_root, "saved_R_data")

qc_root    <- file.path(project_root, "QC")
preqc_dir  <- file.path(qc_root, "Pre-QC")
postqc_dir <- file.path(qc_root, "Post-QC")
for (d in c(out.path, qc_root, preqc_dir, postqc_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ===================== 1) Discover sample folders ============================
get_folder_names <- function(in.path) {
  fn <- list.dirs(in.path, full.names = FALSE, recursive = FALSE)
  fn[fn != ""]
}
f_names <- get_folder_names(in.path)

# keep only folders that actually contain a filtered matrix (drops stray dirs)
has_h5 <- file.exists(file.path(in.path, f_names,
                                "outs", "filtered_feature_bc_matrix.h5"))
f_names <- f_names[has_h5]

print(f_names)
# Expected: "EV_rep1" "Scr_rep1" "miR29a_rep1" "EV_rep2" "Scr_rep2" "miR29a_rep2"
stopifnot(length(f_names) == 6)

# ===================== Load + build per-sample objects =======================
seurat_list <- list()

for (i in f_names) {
  h5 <- file.path(in.path, i, "outs", "filtered_feature_bc_matrix.h5")
  message("Reading: ", i)
  
  mat <- Read10X_h5(h5)
  # A GEX-only run returns the matrix directly; if it's a list (multi feature
  # types), pull the Gene Expression block.
  if (is.list(mat)) mat <- mat[["Gene Expression"]]
  
  s <- CreateSeuratObject(counts = mat, project = i, min.cells = 3)
  
  # ---- sample-level metadata (parsed from the CellRanger --id) ----
  s$orig.ident <- i
  s$condition  <- sub("_rep[0-9]+$", "", i)           # EV / Scr / miR29a
  s$replicate  <- sub("^.*_(rep[0-9]+)$", "\\1", i)    # rep1 / rep2
  
  # ---- QC metrics (MOUSE gene nomenclature — note the case!) ----
  s <- PercentageFeatureSet(s, pattern = "^mt-",       col.name = "percent_mito")
  s <- PercentageFeatureSet(s, pattern = "^Rp[sl]",    col.name = "percent_ribo")
  s <- PercentageFeatureSet(s, pattern = "^Hb[ab]",    col.name = "percent_hb")
  s <- PercentageFeatureSet(s, pattern = "Pecam1|Pf4", col.name = "percent_plat")
  s$log10GenesPerUMI <- log10(s$nFeature_RNA) / log10(s$nCount_RNA)
  
  seurat_list[[i]] <- s
}

# sanity check: if these medians are all 0, the mouse mito pattern didn't match
sapply(seurat_list, function(x) round(median(x$percent_mito), 2))

# ===================== 2) Merge + checkpoint =================================
merged_seurat <- merge(
  x           = seurat_list[[1]],
  y           = seurat_list[-1],
  add.cell.id = names(seurat_list)
)

merged_seurat$condition <- factor(merged_seurat$condition,
                                  levels = c("EV", "Scr", "miR29a"))
merged_seurat$replicate <- factor(merged_seurat$replicate,
                                  levels = c("rep1", "rep2"))

print(table(merged_seurat$orig.ident))

# checkpoint save (lets you restart from the QC step without re-reading h5s)
qs_save(merged_seurat,
        file = file.path(out.path, "Mouse_CARTmiR29a_merged_seurat.qs2"))

# ===================== 3) PRE-QC PLOTS =======================================
setwd(preqc_dir)
metadata <- merged_seurat@meta.data

# GEX-only QC feature set (no ADT features)
feats.1 <- c("nCount_RNA", "nFeature_RNA",
             "percent_mito", "percent_ribo", "percent_hb", "percent_plat")

# Cells per sample
png("Cells_per_sample.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(x = orig.ident, fill = orig.ident)) +
  geom_bar() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title  = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("NCells per sample")
dev.off()

# Cells per condition
png("Cells_per_condition.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(x = condition, fill = condition)) +
  geom_bar() + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("NCells per condition")
dev.off()

# Grouped violins of QC features
png("Pre-QC_features_grouped.png", width = 1800, height = 1200)
VlnPlot(merged_seurat, group.by = "orig.ident", features = feats.1,
        pt.size = 0.1, ncol = 3) + NoLegend()
dev.off()

# UMI count density  (floor = 1000)
png("UMI_Count.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) +
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  ylab("Cell density") + geom_vline(xintercept = 1000)
dev.off()

# nGenes density  (floor = 500)
png("nGenes.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(color = orig.ident, x = nFeature_RNA, fill = orig.ident)) +
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  ylab("Cell density") + geom_vline(xintercept = 500)
dev.off()

# Complexity (log10 genes per UMI)  (floor = 0.75)
png("Complexity_Score.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(x = log10GenesPerUMI, color = orig.ident, fill = orig.ident)) +
  geom_density(alpha = 0.2) + theme_classic() + geom_vline(xintercept = 0.75)
dev.off()

# Mito ratio  (ceiling = 10)
png("Mito_Ratio.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(color = orig.ident, x = percent_mito, fill = orig.ident)) +
  geom_density(alpha = 0.2) + theme_classic() + geom_vline(xintercept = 10)
dev.off()

# Ribo ratio
png("Ribo_Ratio.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(color = orig.ident, x = percent_ribo, fill = orig.ident)) +
  geom_density(alpha = 0.2) + theme_classic()
dev.off()

# Hemoglobin ratio  (ceiling = 5)
png("Heme_Ratio.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(color = orig.ident, x = percent_hb, fill = orig.ident)) +
  geom_density(alpha = 0.2) + theme_classic() + geom_vline(xintercept = 5)
dev.off()

# Platelet ratio
png("Platelet_Ratio.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(color = orig.ident, x = percent_plat, fill = orig.ident)) +
  geom_density(alpha = 0.2) + theme_classic()
dev.off()

# scCustomize grouped cutoff + scatter views (annotations match filter floors)
p1 <- QC_Plots_Genes(merged_seurat, low_cutoff = 500,  high_cutoff = 6000)
p2 <- QC_Plots_UMIs(merged_seurat,  low_cutoff = 1000, high_cutoff = 50000)
p3 <- QC_Plots_Mito(merged_seurat,  high_cutoff = 10)
p4 <- QC_Plots_Complexity(merged_seurat, high_cutoff = 0.75)

png("Grouped_Cutoff.png", width = 1800, height = 1200)
wrap_plots(p1, p2, p3, p4, ncol = 4)
dev.off()

png("UMIvsGene.png", width = 1800, height = 1200)
QC_Plot_UMIvsGene(merged_seurat,
                  low_cutoff_gene = 500, high_cutoff_gene = 6000,
                  low_cutoff_UMI  = 1000, high_cutoff_UMI  = 50000,
                  group.by = "orig.ident")
dev.off()

png("MitovsGene_gradient.png", width = 1800, height = 1200)
QC_Plot_UMIvsGene(merged_seurat,
                  meta_gradient_name = "percent_mito",
                  low_cutoff_gene = 500, high_cutoff_gene = 6000,
                  high_cutoff_UMI = 50000)
dev.off()

# ===================== 4) FILTERING ==========================================
# Thresholds set from this run's Pre-QC plots. Re-inspect if you re-sequence.
# In Seurat v5 the merged object holds one counts layer per sample; JoinLayers
# collapses them into a single matrix so subsetting/gene-filtering behave.
merged_seurat <- JoinLayers(merged_seurat)
DefaultAssay(merged_seurat) <- "RNA"

filtered_seurat <- subset(
  merged_seurat,
  subset =
    nCount_RNA       >= 1000  &   # min UMIs  (well below the ~25k main mass)
    nFeature_RNA     >= 500   &   # min detected genes
    log10GenesPerUMI >  0.75  &   # complexity floor
    percent_mito     <  10    &   # dying/stressed cells (prep is clean, ~2%)
    percent_hb       <  5         # RBC contamination
  # Optional extras — uncomment / adjust if needed:
  # & nCount_RNA     >= 2000      # firmer floor; only trims a thin low tail
  # & percent_ribo   >  5
  # & percent_plat   <  2
  # NOTE: do NOT add an nFeature ceiling for doublets — handle those in
  # script 03 with scDblFinder so you don't delete deeply-sequenced real cells.
)

# Keep genes detected in >= 10 cells
counts     <- GetAssayData(filtered_seurat, assay = "RNA", layer = "counts")
keep_genes <- rownames(counts)[Matrix::rowSums(counts > 0) >= 10]
filtered_seurat <- subset(filtered_seurat, features = keep_genes)

message("Cells before: ", ncol(merged_seurat),
        " | after: ", ncol(filtered_seurat))

# ===================== 5) POST-QC PLOTS + save ===============================
setwd(postqc_dir)

png("Post-QC_Cells_per_sample.png", width = 1800, height = 1200)
filtered_seurat@meta.data %>%
  ggplot(aes(x = orig.ident, fill = orig.ident)) +
  geom_bar() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title  = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Post-QC NCells per sample")
dev.off()

png("Post-QC_Cells_per_condition.png", width = 1800, height = 1200)
filtered_seurat@meta.data %>%
  ggplot(aes(x = condition, fill = condition)) +
  geom_bar() + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Post-QC NCells per condition")
dev.off()

png("Post-QC_features_grouped.png", width = 1800, height = 1200)
VlnPlot(filtered_seurat, group.by = "orig.ident", features = feats.1,
        pt.size = 0.1, ncol = 3) + NoLegend()
dev.off()

qs_save(filtered_seurat,
        file = file.path(out.path, "Mouse_CARTmiR29a_filtered_seurat.qs2"))

message("Saved filtered object -> ",
        file.path(out.path, "Mouse_CARTmiR29a_filtered_seurat.qs2"))
###############################################################################