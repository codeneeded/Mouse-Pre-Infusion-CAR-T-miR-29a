###############################################################################
# 11_subcluster_CD8_CD4_mouse_CARTmiR29a.R
#
# Subset the CAR-T miR-29a object into CD8 and CD4 lineages and re-process each
# from scratch (HVG -> scale -> PCA -> Harmony -> UMAP -> clustering) so that
# within-lineage structure (naive / effector / memory / exhausted, etc.) is
# resolved instead of the global lineage axis.
#
# This script ONLY produces the re-embedded, re-clustered objects + QC and the
# markers for review. Annotation is script 12; trajectory is script 13.
#
# INTENDED WORKFLOW (per lineage): run sections 1-6 (subset -> embed -> harmony
# -> cluster SWEEP -> clustree -> QC UMAPs), then review clustree_stability.png
# and the sweep UMAPs, set default_res_cd8 / default_res_cd4 to the chosen
# resolution, and run section 7 to compute markers + average expression ONCE on
# that resolution. The sweep clustering is cheap; FindAllMarkers is not, so it
# is deliberately NOT run across every swept resolution.
#
# STRUCTURE: deliberately FLAT. Shared config at the top, then CD8 spelled out
# in numbered sections top-to-bottom, then CD4 the same way. No per-lineage
# wrapper function -- run/inspect one section at a time, and both objects
# (seu_cd8, seu_cd4) stay in the environment for debugging.
#
# Why re-embed (not reuse umap.harmony):
#   On the full object the HVGs are dominated by lineage-discriminating genes
#   (Cd4/Cd8, contaminants). On a single-lineage subset those are no longer
#   variable; recomputing HVGs surfaces the within-lineage differentiation
#   axis. Everything downstream of HVGs follows, and Harmony runs on the new
#   PCA, so it must re-run too.
#
# Normalization note: the object is RNA / LogNormalize. The `data` slot is
# per-cell, so subsetting does NOT require re-running NormalizeData. ScaleData
# onward DOES. Module scores already in the object survive subsetting.
#
# Inputs:
#   saved_R_data/Mouse_CARTmiR29a_WithModuleScores.qs2
#
# Outputs:
#   saved_R_data/Mouse_CARTmiR29a_CD8_subclustered.qs2
#   saved_R_data/Mouse_CARTmiR29a_CD4_subclustered.qs2
#   Subclustering/
#     lineage_assignment.csv              tentative_state -> CD8/CD4 map (review!)
#     CD8/  CD4/
#       QC/                               elbow, clustree (stability + gene
#                                         overlays), resolution-sweep UMAPs,
#                                         condition/sample UMAPs, lineage-check
#                                         feature plots
#       Markers/                          FindAllMarkers CSV + top-N per
#                                         cluster + avg expression, for the
#                                         CHOSEN resolution only
###############################################################################

library(Seurat)
library(SeuratExtend)    # DimPlot2, theme_umap_arrows
library(scCustomize)     # FeaturePlot_scCustom
library(harmony)         # RunHarmony
library(clustree)        # resolution-sweep stability tree
library(qs2)
library(dplyr); library(tidyr); library(tibble)
library(ggplot2); library(patchwork); library(viridis)

# ============================ Paths ==========================================
project_dir <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/10x_scRNAseq"
saved_dir   <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/saved_R_data"
sub_base    <- file.path(project_dir, "Subclustering")
dir.create(sub_base, recursive = TRUE, showWarnings = FALSE)

# ============================ Parameters (REVIEW) ============================
# Set these to match the ORIGINAL integration script (02/03). Harmony should
# correct the technical batch (replicate), NOT condition -- correcting on
# condition would remove the EV/Scr/miR29a biology.
harmony_group_var <- "replicate"          # <- confirm against original script
n_dims_cd8        <- 15                   # CD8 elbow flattens ~PC 10-15 (sparse, 2 parent states)
n_dims_cd4        <- 30                   # CD4 has more structure; check its elbow + adjust
n_hvg             <- 2000                 # FindVariableFeatures nfeatures
# Resolution sweeps are PER-LINEAGE. CD8 here is sparse (few parent states) so
# it gets a lower, finer range to avoid over-splitting; CD4 has many parent
# states so it gets a higher range. Tune after seeing each clustree.
res_sweep_cd8   <- c(0.1, 0.2, 0.3, 0.4, 0.5)
default_res_cd8 <- 0.5                    # resolution written to CD8 $subclusters
res_sweep_cd4   <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
default_res_cd4 <- 0.6                    # resolution written to CD4 $subclusters
# Section 7 computes markers + avg expression for EVERY resolution in these
# compare vectors (deliberate, bounded comparison -- NOT the whole sweep).
# Files are tagged with the resolution so candidates don't overwrite each
# other. $subclusters is still set from default_res_* (the committed choice).
compare_res_cd8 <- c(0.4, 0.5)            # CD8: compare these two side by side
compare_res_cd4 <- c(default_res_cd4)     # CD4: single (add more to compare)
scale_vars        <- NULL                 # ScaleData vars.to.regress, or NULL
# Genes overlaid on clustree nodes (median expression) so you can watch a
# state-defining marker resolve as resolution climbs. Present-only filtered.
# Per-lineage: CD8 differentiation axis vs CD4 helper/Treg axis.
clustree_overlay_genes_cd8 <- c("Sell", "Ccr7", "Tcf7", "Gzmb", "Mki67", "Havcr2")
clustree_overlay_genes_cd4 <- c("Foxp3", "Il2ra", "Tbx21", "Ifng",
                                "Cxcr5", "Bcl6", "Tcf7", "Mki67")
# Lineage-confirmation genes for QC feature plots (verify subset purity)
lineage_check_genes <- c("Cd3e", "Cd8a", "Cd8b1", "Cd4", "Foxp3", "Mki67")

cond_cols <- c(EV = "#E76F51", Scr = "#52B788", miR29a = "#5E60CE")

# ============================ Load object ====================================
message("Loading object...")
obj <- qs_read(file.path(saved_dir, "Mouse_CARTmiR29a_WithModuleScores.qs2"))
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")
if (!"data" %in% Layers(obj[["RNA"]]))
  obj <- NormalizeData(obj, verbose = FALSE)

obj$condition       <- factor(obj$condition, levels = c("EV", "Scr", "miR29a"))
obj$replicate       <- factor(obj$replicate)
obj$.sample         <- paste(obj$condition, obj$replicate, sep = "-")
obj$tentative_state <- droplevels(factor(obj$tentative_state))

stopifnot(harmony_group_var %in% colnames(obj@meta.data))

# ============================================================================
# LINEAGE ASSIGNMENT (shared) -- cluster-level, from tentative_state
# ============================================================================
# Auto-guess each state's lineage by string match, write to CSV for review,
# edit the manual override block for ambiguous states. Each state must end up
# as "CD8", "CD4", or "EXCLUDE"; a leftover NA (truly unassigned) halts the
# script so nothing is silently dropped or misassigned.
message("\n=== Lineage assignment ===")

state_levels <- levels(obj$tentative_state)
auto_lineage <- dplyr::case_when(
  grepl("CD8", state_levels, ignore.case = TRUE) ~ "CD8",
  grepl("CD4", state_levels, ignore.case = TRUE) ~ "CD4",
  TRUE                                            ~ NA_character_
)
lineage_map <- setNames(auto_lineage, state_levels)

# ---- MANUAL OVERRIDE: resolve ambiguous / lineage-less states here ----------
# Set each to "CD8", "CD4", or "EXCLUDE". "EXCLUDE" = intentionally dropped
# from both lineages (distinct from an unassigned NA, which halts the script).
lineage_map["Cycling non-T (review)"]         <- "EXCLUDE"  # non-T
lineage_map["Treg"]                           <- "CD4"      # Treg is a CD4 subset
lineage_map["NK-like / innate-like (review)"] <- "EXCLUDE"  # non-conventional
# (set the last to "CD8" instead if you consider these CD8-lineage innate-like)
# -----------------------------------------------------------------------------

lineage_df <- data.frame(tentative_state = names(lineage_map),
                         lineage         = unname(lineage_map),
                         stringsAsFactors = FALSE)
write.csv(lineage_df, file.path(sub_base, "lineage_assignment.csv"),
          row.names = FALSE)
message("Lineage assignment (review Subclustering/lineage_assignment.csv):")
print(lineage_df)

unassigned <- lineage_df$tentative_state[is.na(lineage_df$lineage)]
if (length(unassigned) > 0) {
  stop("Unassigned states -- set each to CD8, CD4, or EXCLUDE in the MANUAL ",
       "OVERRIDE block: ", paste(unassigned, collapse = ", "))
}

excluded <- lineage_df$tentative_state[lineage_df$lineage == "EXCLUDE"]
if (length(excluded) > 0)
  message(sprintf("Excluded from both lineages (%d): %s",
                  length(excluded), paste(excluded, collapse = ", ")))

states_cd8 <- lineage_df$tentative_state[lineage_df$lineage == "CD8"]
states_cd4 <- lineage_df$tentative_state[lineage_df$lineage == "CD4"]
message(sprintf("CD8 states (%d): %s", length(states_cd8),
                paste(states_cd8, collapse = ", ")))
message(sprintf("CD4 states (%d): %s", length(states_cd4),
                paste(states_cd4, collapse = ", ")))


###############################################################################
###############################################################################
##                                                                           ##
##                              C D 8                                         ##
##                                                                           ##
###############################################################################
###############################################################################

cd8_dir <- file.path(sub_base, "CD8")
cd8_qc  <- file.path(cd8_dir, "QC")
cd8_mk  <- file.path(cd8_dir, "Markers")
dir.create(cd8_qc, recursive = TRUE, showWarnings = FALSE)
dir.create(cd8_mk, recursive = TRUE, showWarnings = FALSE)

# ===== CD8 : 1. SUBSET =======================================================
message("\n=== CD8 : subset ===")
seu_cd8 <- subset(obj, subset = tentative_state %in% states_cd8)
seu_cd8$parent_state    <- droplevels(factor(seu_cd8$tentative_state))
seu_cd8$tentative_state <- NULL
# Drop clustering columns inherited from the global object. subset() carries
# over every RNA_snn_res.* / seurat_clusters column from obj, and FindClusters
# only overwrites the resolutions we sweep -- any others survive as STALE
# (global-clustering) columns that clustree would otherwise pick up. Strip them
# so the tree reflects only THIS subset's re-clustering.
old_clust_cd8 <- grep("^RNA_snn_res\\.|^seurat_clusters$",
                      colnames(seu_cd8@meta.data), value = TRUE)
for (cc in old_clust_cd8) seu_cd8[[cc]] <- NULL
message(sprintf("  CD8: %d cells across %d parent states",
                ncol(seu_cd8), nlevels(seu_cd8$parent_state)))

# ===== CD8 : 2. RE-EMBED (HVG -> scale -> PCA) ===============================
message("\n=== CD8 : re-embed ===")
seu_cd8 <- FindVariableFeatures(seu_cd8, selection.method = "vst",
                                nfeatures = n_hvg, verbose = FALSE)
seu_cd8 <- ScaleData(seu_cd8, vars.to.regress = scale_vars, verbose = FALSE)
seu_cd8 <- RunPCA(seu_cd8, npcs = max(50, n_dims_cd8), verbose = FALSE)
ggsave(file.path(cd8_qc, "elbow_PCA.png"),
       ElbowPlot(seu_cd8, ndims = max(50, n_dims_cd8)),
       width = 7, height = 5, dpi = 300, bg = "white")

# ===== CD8 : 3. HARMONY -> NEIGHBORS -> UMAP -> CLUSTER SWEEP =================
message("\n=== CD8 : harmony + cluster sweep ===")
seu_cd8 <- RunHarmony(seu_cd8, group.by.vars = harmony_group_var,
                      reduction = "pca", dims.use = 1:n_dims_cd8,
                      reduction.save = "harmony", verbose = FALSE)
seu_cd8 <- FindNeighbors(seu_cd8, reduction = "harmony", dims = 1:n_dims_cd8,
                         verbose = FALSE)
seu_cd8 <- RunUMAP(seu_cd8, reduction = "harmony", dims = 1:n_dims_cd8,
                   reduction.name = "umap.harmony", verbose = FALSE)
seu_cd8 <- FindClusters(seu_cd8, resolution = res_sweep_cd8, verbose = FALSE)

# Default-resolution clustering -> $subclusters
default_col <- paste0("RNA_snn_res.", default_res_cd8)
if (!default_col %in% colnames(seu_cd8@meta.data))
  stop("CD8 default_res_cd8 column not found: ", default_col)
seu_cd8$subclusters <- seu_cd8[[default_col]][, 1]
Idents(seu_cd8) <- "subclusters"

# ===== CD8 : 4. CLUSTREE =====================================================
message("\n=== CD8 : clustree ===")
ggsave(file.path(cd8_qc, "clustree_stability.png"),
       clustree(seu_cd8, prefix = "RNA_snn_res."),
       width = 10, height = 9, dpi = 300, bg = "white", limitsize = FALSE)
for (g in intersect(clustree_overlay_genes_cd8, rownames(seu_cd8))) {
  ggsave(file.path(cd8_qc, paste0("clustree_gene_", g, ".png")),
         clustree(seu_cd8, prefix = "RNA_snn_res.",
                  node_colour = g, node_colour_aggr = "median"),
         width = 10, height = 9, dpi = 300, bg = "white", limitsize = FALSE)
}

# ===== CD8 : 5. QC UMAPS =====================================================
message("\n=== CD8 : QC UMAPs ===")
for (r in res_sweep_cd8) {
  rc <- paste0("RNA_snn_res.", r)
  ggsave(file.path(cd8_qc, sprintf("UMAP_clusters_res%.1f.png", r)),
         DimPlot2(seu_cd8, features = rc, reduction = "umap.harmony",
                  label = TRUE, box = TRUE, theme = theme_umap_arrows()) +
           ggtitle(sprintf("CD8 -- clusters @ res %.1f", r)),
         width = 9, height = 8, dpi = 300, bg = "white")
}
ggsave(file.path(cd8_qc, "UMAP_by_condition.png"),
       DimPlot2(seu_cd8, features = "condition", reduction = "umap.harmony",
                cols = cond_cols, box = TRUE, theme = theme_umap_arrows()) +
         ggtitle("CD8 -- condition"),
       width = 9, height = 8, dpi = 300, bg = "white")
ggsave(file.path(cd8_qc, "UMAP_by_sample.png"),
       DimPlot2(seu_cd8, features = ".sample", reduction = "umap.harmony",
                box = TRUE, theme = theme_umap_arrows()) + ggtitle("CD8 -- sample"),
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(cd8_qc, "UMAP_by_parent_state.png"),
       DimPlot2(seu_cd8, features = "parent_state", reduction = "umap.harmony",
                label = TRUE, box = TRUE, theme = theme_umap_arrows()) +
         ggtitle("CD8 -- parent (global) state"),
       width = 11, height = 8, dpi = 300, bg = "white")

# ===== CD8 : 6. LINEAGE-CHECK FEATURE PLOTS ==================================
cd8_check <- intersect(lineage_check_genes, rownames(seu_cd8))
if (length(cd8_check) > 0) {
  cd8_panels <- lapply(cd8_check, function(g)
    FeaturePlot_scCustom(seu_cd8, reduction = "umap.harmony", features = g,
                         colors_use = viridis(10, option = "A"), order = TRUE) +
      theme(legend.position = "none", plot.title = element_text(face = "italic")))
  cd8_nc <- min(3, length(cd8_panels))
  ggsave(file.path(cd8_qc, "UMAP_lineage_check_genes.png"),
         wrap_plots(cd8_panels, ncol = cd8_nc),
         width = 5 * cd8_nc,
         height = 4.5 * ceiling(length(cd8_panels) / cd8_nc),
         dpi = 300, bg = "white", limitsize = FALSE)
}

# ===== CD8 : 7. MARKERS + AVERAGE EXPRESSION (compare_res_cd8) ===============
# Run AFTER reviewing clustree + the sweep UMAPs. Computes markers + avg
# expression for EACH resolution in compare_res_cd8 (here 0.4 and 0.5) so they
# can be compared; files are tagged with the resolution. $subclusters is set
# from default_res_cd8 (the committed choice for downstream / the saved object).
message("\n=== CD8 : markers @ compare_res_cd8 ===")
for (r in compare_res_cd8) {
  rc <- paste0("RNA_snn_res.", r)
  if (!rc %in% colnames(seu_cd8@meta.data)) {
    message("  skip res ", r, " (column ", rc, " not in object -- add to res_sweep_cd8)")
    next
  }
  Idents(seu_cd8) <- rc
  message(sprintf("  res %.2f -> %d clusters", r, length(levels(seu_cd8[[rc]][, 1]))))
  
  mk <- FindAllMarkers(seu_cd8, only.pos = TRUE, min.pct = 0.25,
                       logfc.threshold = 0.25, verbose = FALSE)
  write.csv(mk, file.path(cd8_mk, sprintf("FindAllMarkers_res%.2f.csv", r)),
            row.names = FALSE)
  if (nrow(mk) > 0) {
    top_mk <- mk %>% dplyr::group_by(cluster) %>%
      dplyr::slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
      dplyr::ungroup()
    write.csv(top_mk, file.path(cd8_mk, sprintf("top20_per_cluster_res%.2f.csv", r)),
              row.names = FALSE)
  }
  
  avg <- AverageExpression(seu_cd8, assays = "RNA", layer = "data",
                           group.by = rc)$RNA
  write.csv(as.data.frame(avg),
            file.path(cd8_mk, sprintf("avg_expression_per_cluster_res%.2f.csv", r)))
  message(sprintf("    %d marker rows; avg-expr %d genes x %d clusters",
                  nrow(mk), nrow(avg), ncol(avg)))
}

# Commit $subclusters to the chosen default resolution
seu_cd8$subclusters <- seu_cd8[[paste0("RNA_snn_res.", default_res_cd8)]][, 1]
Idents(seu_cd8) <- "subclusters"

# ===== CD8 : 8. README + SAVE ================================================
writeLines(c(
  "CD8 subclustering (script 11)",
  paste0("Cells: ", ncol(seu_cd8)),
  paste0("Parent states pulled in: ", paste(states_cd8, collapse = ", ")),
  paste0("Harmony batch var: ", harmony_group_var, "; n_dims_cd8: ", n_dims_cd8,
         "; HVGs: ", n_hvg),
  paste0("Resolutions swept: ", paste(res_sweep_cd8, collapse = ", ")),
  paste0("$subclusters set from res ", default_res_cd8,
         " (change default_res_cd8 + re-run, or repoint in script 12)."),
  "Annotation = script 12; trajectory = script 13."
), file.path(cd8_dir, "README.txt"))
qs_save(seu_cd8, file.path(saved_dir, "Mouse_CARTmiR29a_CD8_subclustered.qs2"))
message("  saved CD8 object.")


###############################################################################
###############################################################################
##                                                                           ##
##                              C D 4                                         ##
##                                                                           ##
###############################################################################
###############################################################################

cd4_dir <- file.path(sub_base, "CD4")
cd4_qc  <- file.path(cd4_dir, "QC")
cd4_mk  <- file.path(cd4_dir, "Markers")
dir.create(cd4_qc, recursive = TRUE, showWarnings = FALSE)
dir.create(cd4_mk, recursive = TRUE, showWarnings = FALSE)

# ===== CD4 : 1. SUBSET =======================================================
message("\n=== CD4 : subset ===")
seu_cd4 <- subset(obj, subset = tentative_state %in% states_cd4)
seu_cd4$parent_state    <- droplevels(factor(seu_cd4$tentative_state))
seu_cd4$tentative_state <- NULL
# Drop clustering columns inherited from the global object (see CD8 note).
old_clust_cd4 <- grep("^RNA_snn_res\\.|^seurat_clusters$",
                      colnames(seu_cd4@meta.data), value = TRUE)
for (cc in old_clust_cd4) seu_cd4[[cc]] <- NULL
message(sprintf("  CD4: %d cells across %d parent states",
                ncol(seu_cd4), nlevels(seu_cd4$parent_state)))

# ===== CD4 : 2. RE-EMBED (HVG -> scale -> PCA) ===============================
message("\n=== CD4 : re-embed ===")
seu_cd4 <- FindVariableFeatures(seu_cd4, selection.method = "vst",
                                nfeatures = n_hvg, verbose = FALSE)
seu_cd4 <- ScaleData(seu_cd4, vars.to.regress = scale_vars, verbose = FALSE)
seu_cd4 <- RunPCA(seu_cd4, npcs = max(50, n_dims_cd4), verbose = FALSE)
ggsave(file.path(cd4_qc, "elbow_PCA.png"),
       ElbowPlot(seu_cd4, ndims = max(50, n_dims_cd4)),
       width = 7, height = 5, dpi = 300, bg = "white")

# ===== CD4 : 3. HARMONY -> NEIGHBORS -> UMAP -> CLUSTER SWEEP =================
message("\n=== CD4 : harmony + cluster sweep ===")
seu_cd4 <- RunHarmony(seu_cd4, group.by.vars = harmony_group_var,
                      reduction = "pca", dims.use = 1:n_dims_cd4,
                      reduction.save = "harmony", verbose = FALSE)
seu_cd4 <- FindNeighbors(seu_cd4, reduction = "harmony", dims = 1:n_dims_cd4,
                         verbose = FALSE)
seu_cd4 <- RunUMAP(seu_cd4, reduction = "harmony", dims = 1:n_dims_cd4,
                   reduction.name = "umap.harmony", verbose = FALSE)
seu_cd4 <- FindClusters(seu_cd4, resolution = res_sweep_cd4, verbose = FALSE)

default_col <- paste0("RNA_snn_res.", default_res_cd4)
if (!default_col %in% colnames(seu_cd4@meta.data))
  stop("CD4 default_res_cd4 column not found: ", default_col)
seu_cd4$subclusters <- seu_cd4[[default_col]][, 1]
Idents(seu_cd4) <- "subclusters"

# ===== CD4 : 4. CLUSTREE =====================================================
message("\n=== CD4 : clustree ===")
ggsave(file.path(cd4_qc, "clustree_stability.png"),
       clustree(seu_cd4, prefix = "RNA_snn_res."),
       width = 10, height = 9, dpi = 300, bg = "white", limitsize = FALSE)
for (g in intersect(clustree_overlay_genes_cd4, rownames(seu_cd4))) {
  ggsave(file.path(cd4_qc, paste0("clustree_gene_", g, ".png")),
         clustree(seu_cd4, prefix = "RNA_snn_res.",
                  node_colour = g, node_colour_aggr = "median"),
         width = 10, height = 9, dpi = 300, bg = "white", limitsize = FALSE)
}

# ===== CD4 : 5. QC UMAPS =====================================================
message("\n=== CD4 : QC UMAPs ===")
for (r in res_sweep_cd4) {
  rc <- paste0("RNA_snn_res.", r)
  ggsave(file.path(cd4_qc, sprintf("UMAP_clusters_res%.1f.png", r)),
         DimPlot2(seu_cd4, features = rc, reduction = "umap.harmony",
                  label = TRUE, box = TRUE, theme = theme_umap_arrows()) +
           ggtitle(sprintf("CD4 -- clusters @ res %.1f", r)),
         width = 9, height = 8, dpi = 300, bg = "white")
}
ggsave(file.path(cd4_qc, "UMAP_by_condition.png"),
       DimPlot2(seu_cd4, features = "condition", reduction = "umap.harmony",
                cols = cond_cols, box = TRUE, theme = theme_umap_arrows()) +
         ggtitle("CD4 -- condition"),
       width = 9, height = 8, dpi = 300, bg = "white")
ggsave(file.path(cd4_qc, "UMAP_by_sample.png"),
       DimPlot2(seu_cd4, features = ".sample", reduction = "umap.harmony",
                box = TRUE, theme = theme_umap_arrows()) + ggtitle("CD4 -- sample"),
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(cd4_qc, "UMAP_by_parent_state.png"),
       DimPlot2(seu_cd4, features = "parent_state", reduction = "umap.harmony",
                label = TRUE, box = TRUE, theme = theme_umap_arrows()) +
         ggtitle("CD4 -- parent (global) state"),
       width = 11, height = 8, dpi = 300, bg = "white")

# ===== CD4 : 6. LINEAGE-CHECK FEATURE PLOTS ==================================
cd4_check <- intersect(lineage_check_genes, rownames(seu_cd4))
if (length(cd4_check) > 0) {
  cd4_panels <- lapply(cd4_check, function(g)
    FeaturePlot_scCustom(seu_cd4, reduction = "umap.harmony", features = g,
                         colors_use = viridis(10, option = "A"), order = TRUE) +
      theme(legend.position = "none", plot.title = element_text(face = "italic")))
  cd4_nc <- min(3, length(cd4_panels))
  ggsave(file.path(cd4_qc, "UMAP_lineage_check_genes.png"),
         wrap_plots(cd4_panels, ncol = cd4_nc),
         width = 5 * cd4_nc,
         height = 4.5 * ceiling(length(cd4_panels) / cd4_nc),
         dpi = 300, bg = "white", limitsize = FALSE)
}

# ===== CD4 : 7. MARKERS + AVERAGE EXPRESSION (compare_res_cd4) ===============
# Same structure as CD8: computes markers + avg expression for each resolution
# in compare_res_cd4 (single value by default). $subclusters set from
# default_res_cd4. Add resolutions to compare_res_cd4 to compare candidates.
message("\n=== CD4 : markers @ compare_res_cd4 ===")
for (r in compare_res_cd4) {
  rc <- paste0("RNA_snn_res.", r)
  if (!rc %in% colnames(seu_cd4@meta.data)) {
    message("  skip res ", r, " (column ", rc, " not in object -- add to res_sweep_cd4)")
    next
  }
  Idents(seu_cd4) <- rc
  message(sprintf("  res %.2f -> %d clusters", r, length(levels(seu_cd4[[rc]][, 1]))))
  
  mk <- FindAllMarkers(seu_cd4, only.pos = TRUE, min.pct = 0.25,
                       logfc.threshold = 0.25, verbose = FALSE)
  write.csv(mk, file.path(cd4_mk, sprintf("FindAllMarkers_res%.2f.csv", r)),
            row.names = FALSE)
  if (nrow(mk) > 0) {
    top_mk <- mk %>% dplyr::group_by(cluster) %>%
      dplyr::slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
      dplyr::ungroup()
    write.csv(top_mk, file.path(cd4_mk, sprintf("top20_per_cluster_res%.2f.csv", r)),
              row.names = FALSE)
  }
  
  avg <- AverageExpression(seu_cd4, assays = "RNA", layer = "data",
                           group.by = rc)$RNA
  write.csv(as.data.frame(avg),
            file.path(cd4_mk, sprintf("avg_expression_per_cluster_res%.2f.csv", r)))
  message(sprintf("    %d marker rows; avg-expr %d genes x %d clusters",
                  nrow(mk), nrow(avg), ncol(avg)))
}

# Commit $subclusters to the chosen default resolution
seu_cd4$subclusters <- seu_cd4[[paste0("RNA_snn_res.", default_res_cd4)]][, 1]
Idents(seu_cd4) <- "subclusters"

# ===== CD4 : 8. README + SAVE ================================================
writeLines(c(
  "CD4 subclustering (script 11)",
  paste0("Cells: ", ncol(seu_cd4)),
  paste0("Parent states pulled in: ", paste(states_cd4, collapse = ", ")),
  paste0("Harmony batch var: ", harmony_group_var, "; n_dims_cd4: ", n_dims_cd4,
         "; HVGs: ", n_hvg),
  paste0("Resolutions swept: ", paste(res_sweep_cd4, collapse = ", ")),
  paste0("$subclusters set from res ", default_res_cd4,
         " (change default_res_cd4 + re-run, or repoint in script 12)."),
  "Annotation = script 12; trajectory = script 13."
), file.path(cd4_dir, "README.txt"))
qs_save(seu_cd4, file.path(saved_dir, "Mouse_CARTmiR29a_CD4_subclustered.qs2"))
message("  saved CD4 object.")


# ============================================================================
# Done
# ============================================================================
message("\nDone.")
message("Review Subclustering/<lineage>/QC/clustree_*.png + the resolution-sweep")
message("UMAPs, then pick a resolution and either re-run with the matching")
message("default_res_cd8 / default_res_cd4 set, or point script 12 at the")
message("chosen RNA_snn_res.<r> column.")
message("Outputs: ", sub_base)
###############################################################################