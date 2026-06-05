###############################################################################
## 04_annotation_markerplots_mouse_CARTmiR29a.R
## Murine pre-infusion CAR-T (miR-29a) — single-cell, GEX only.
## Pre-annotation marker visualization on a chosen integrated object.
##
## Produces the views used to MANUALLY assign cluster identities:
##   - DimPlots (clusters, condition), cluster-distribution bars
##   - a summary DotPlot of the marker panel across clusters
##   - per-gene Violin (by cluster + split by condition) and Feature plots
##
## NB this is a CD4/CD8-sorted T-cell product, so the panel is T-cell focused
## (not the human PBMC panel from the CITE-seq workflow). GEX only -> no ADT.
###############################################################################

# ---- Libraries ----
library(Seurat)
library(SeuratExtend)   # DimPlot2, VlnPlot2, DotPlot2, ClusterDistrBar
library(scCustomize)    # FeaturePlot_scCustom
library(tidyverse)
library(patchwork)
library(viridis)
library(pheatmap)       # publication marker heatmap
library(qs2)

# ============================ Choices ========================================
# Which integrated object to annotate (decide via composition + referee first).
integration_choice <- "by_replicate"          # "by_replicate" or "by_origident"
umap_reduction     <- "umap.harmony"           # match the kept integration method
min_cluster_size   <- 20                        # drop clusters smaller than this

# ====================== RESOLUTION SELECTION (clustree) ======================
## clustree is generated in script 03 (it needs the integrated graph) and ALL
## sweep columns (RNA_snn_res.0.2 ... 1.2) are saved in the object. Do NOT
## re-cluster here: inspect Integration/<choice>/Clustree/clustree.png, then set
## the chosen column below. THIS is the single place the resolution is committed.
res_choice <- "RNA_snn_res.0.4"                 # <-- chosen from clustree (0.4 is stable;
#     0.6's extra clusters were unstable subdivisions).
#     Confirm it resolves CD4/CD8 + states in the DotPlot;
#     subcluster a lineage later if 0.4 lumps two states.

# ============================ Paths ==========================================
project_root <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a"
saved_dir    <- file.path(project_root, "saved_R_data")
annot_dir    <- file.path(project_root, "Annotation", "Pre-Annotation")
clust_plot_dir <- file.path(annot_dir, "Cluster_Plots")
rna_vln_dir       <- file.path(annot_dir, "Violin_Plot", "RNA")
rna_vln_cond_dir  <- file.path(annot_dir, "Violin_Plot", "RNA_splitBy_condition")
rna_feature_dir   <- file.path(annot_dir, "Feature_Plot", "RNA")
for (d in c(clust_plot_dir, rna_vln_dir, rna_vln_cond_dir, rna_feature_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Load chosen integrated object ----
obj <- qs_read(file.path(saved_dir,
                         paste0("Mouse_CARTmiR29a_integrated_", integration_choice, ".qs2")))
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])        # ensure single layer for plotting/DE

# ---- commit resolution + drop tiny clusters ----
stopifnot(res_choice %in% colnames(obj@meta.data))
Idents(obj)  <- res_choice
obj$clusters <- Idents(obj)

sizes <- table(obj$clusters)
keep  <- names(sizes[sizes >= min_cluster_size])
message("Clusters kept (>= ", min_cluster_size, " cells): ",
        length(keep), "/", length(sizes))
obj <- subset(obj, idents = keep)
obj$clusters <- droplevels(factor(obj$clusters))
# order cluster levels NUMERICALLY (0,1,2,...,10,11 — not the default 0,1,10,11,...,2)
lv <- levels(obj$clusters)
obj$clusters <- factor(obj$clusters, levels = lv[order(as.numeric(lv))])
Idents(obj)  <- "clusters"

# ====================== Per-cell lineage assignment (INTERIM) ================
# CD4/CD8-sorted product. The proper fix (re-cluster WITHIN each lineage) is
# deferred for compute reasons. As an interim, assign lineage PER CELL from
# marker expression so CD4 vs CD8 can be separated even inside lineage-mixed
# clusters (e.g. the proliferating ones) — with NO re-run of script 03.
# DN = dropout (no Cd4/Cd8 detected), concentrated in cycling cells; those get
# resolved properly when we re-cluster per lineage later.
cd8_genes <- intersect(c("Cd8a", "Cd8b1", "Cd8b"), rownames(obj[["RNA"]]))
cd4_genes <- intersect(c("Cd4"),                   rownames(obj[["RNA"]]))
gd_genes  <- intersect(c("Trdc", "Trgc1", "Trgc2"), rownames(obj[["RNA"]]))
lin <- FetchData(obj, vars = unique(c(cd8_genes, cd4_genes, gd_genes)),
                 layer = "data")
cd8_score <- if (length(cd8_genes)) do.call(pmax, lin[cd8_genes]) else 0
cd4_score <- if (length(cd4_genes)) lin[[cd4_genes[1]]]           else 0
gd_score  <- if (length(gd_genes))  do.call(pmax, lin[gd_genes])  else 0

obj$lineage <- dplyr::case_when(
  gd_score  > 0 & gd_score  >= cd8_score & gd_score >= cd4_score ~ "gdT",
  cd8_score > 0 & cd8_score >= cd4_score                         ~ "CD8",
  cd4_score > 0 & cd4_score >  cd8_score                         ~ "CD4",
  TRUE                                                            ~ "DN"
)
obj$lineage <- factor(obj$lineage, levels = c("CD4", "CD8", "gdT", "DN"))
message("Per-cell lineage assignment:")
print(table(obj$lineage))

# lineage x cluster cross-tab -> shows which clusters are pure vs mixed
lin_tab <- obj@meta.data %>%
  dplyr::count(clusters, lineage) %>%
  tidyr::pivot_wider(names_from = lineage, values_from = n, values_fill = 0)
write.csv(lin_tab, file.path(clust_plot_dir, "lineage_by_cluster.csv"),
          row.names = FALSE)

# per-cluster majority call (among lineage-assigned cells): CD8 / CD4 / Mixed
clust_call <- obj@meta.data %>%
  dplyr::filter(lineage %in% c("CD4", "CD8")) %>%
  dplyr::count(clusters, lineage) %>%
  dplyr::group_by(clusters) %>%
  dplyr::summarise(CD8_frac = sum(n[lineage == "CD8"]) / sum(n), .groups = "drop") %>%
  dplyr::mutate(call = dplyr::case_when(CD8_frac >= 0.70 ~ "CD8",
                                        CD8_frac <= 0.30 ~ "CD4",
                                        TRUE             ~ "Mixed"))
write.csv(clust_call, file.path(clust_plot_dir, "cluster_lineage_call.csv"),
          row.names = FALSE)
print(clust_call)

# plots: lineage on UMAP + lineage composition per cluster
ggsave(file.path(clust_plot_dir, "UMAP_lineage.png"),
       DimPlot2(obj, reduction = umap_reduction, group.by = "lineage", pt.size = 0.4),
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(clust_plot_dir, "Cluster_lineage_composition.png"),
       ClusterDistrBar(obj$lineage, obj$clusters, flip = FALSE, border = "black") +
         theme(axis.title.x = element_blank()),
       width = 14, height = 8, dpi = 300, bg = "white")

# ============================ Marker panel (MOUSE) ===========================
# Mouse symbols (title case). Note mouse-specific names: CD8b = Cd8b1.
marker_panel <- list(
  Pan_T          = c("Cd3d","Cd3e","Cd3g","Cd247"),
  Lineage        = c("Cd4","Cd8a","Cd8b1"),
  Naive_Memory   = c("Tcf7","Sell","Ccr7","Lef1","Il7r","Bach2","Slamf6","Cd27","Cd28"),
  Effector_Cyto  = c("Gzmb","Gzmk","Gzma","Prf1","Nkg7","Klrg1","Cx3cr1","Zeb2","Klrd1"),
  Activation     = c("Cd69","Il2ra","Tnfrsf9","Tnfrsf4","Icos","Cd44","Tnfrsf18"),
  Exhaustion     = c("Pdcd1","Havcr2","Lag3","Tigit","Ctla4","Entpd1","Tox","Nr4a1","Nr4a2"),
  Treg           = c("Foxp3","Ikzf2"),
  Lineage_TFs    = c("Tbx21","Gata3","Rorc","Bcl6","Prdm1","Maf","Id2","Id3","Eomes"),
  Cytokines      = c("Ifng","Tnf","Il2","Il10","Il17a","Il21","Csf2"),
  Proliferation  = c("Mki67","Top2a","Pcna","Birc5","Stmn1"),
  miR29a_targets = c("Tbx21","Eomes","Ifng","Dnmt3a","Dnmt3b"),  # canonical miR-29 targets
  Purity_check   = c("Cd19","Cd14","Itgam","Ncr1","Klrb1c")      # should be ~absent if sort is clean
)

# flatten, keep only genes present in the data
rna.features <- unique(unlist(marker_panel))
missing <- setdiff(rna.features, rownames(obj[["RNA"]]))
if (length(missing)) message("Markers not in data (skipped): ",
                             paste(missing, collapse = ", "))
rna.features <- intersect(rna.features, rownames(obj[["RNA"]]))

# ============================ Cluster overview plots =========================
ggsave(file.path(clust_plot_dir, paste0("UMAP_", res_choice, ".png")),
       DimPlot2(obj, reduction = umap_reduction, group.by = "clusters",
                label = TRUE, repel = TRUE, box = TRUE, label.size = 5) +
         ggtitle(res_choice),
       width = 10, height = 8, dpi = 300, bg = "white")

ggsave(file.path(clust_plot_dir, "UMAP_clusters_vs_condition.png"),
       wrap_plots(
         DimPlot2(obj, reduction = umap_reduction, group.by = "clusters",
                  label = TRUE, repel = TRUE, label.size = 4) + ggtitle("Clusters"),
         DimPlot2(obj, reduction = umap_reduction, group.by = "condition") +
           ggtitle("Condition"),
         ncol = 2),
       width = 18, height = 8, dpi = 300, bg = "white")

# cluster composition by sample and by condition (was OUD_status)
ggsave(file.path(clust_plot_dir, "Cluster_Distribution_by_sample.png"),
       ClusterDistrBar(obj$orig.ident, obj$clusters, flip = FALSE, border = "black") +
         theme(axis.title.x = element_blank()),
       width = 16, height = 9, dpi = 300, bg = "white")

ggsave(file.path(clust_plot_dir, "Cluster_Distribution_by_condition.png"),
       ClusterDistrBar(obj$condition, obj$clusters, flip = FALSE, border = "black") +
         theme(axis.title.x = element_blank()),
       width = 13, height = 9, dpi = 300, bg = "white")

# summary DotPlot — usually the fastest way to assign identities.
# ~70 genes on the y-axis: height must scale with gene count or it squishes.
ggsave(file.path(clust_plot_dir, "DotPlot_marker_panel.png"),
       DotPlot2(obj, features = marker_panel, group.by = "clusters") +
         theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)),
       width = 14, height = 22, dpi = 300, bg = "white", limitsize = FALSE)

# checkpoint
qs_save(obj, file = file.path(saved_dir, "Mouse_CARTmiR29a_PreAnnotation.qs2"))

# ===================== Publication heatmap (Fig-1C style) ====================
# RNA-only (no ADT). Mean expression per cluster, row z-scored, rows grouped by
# functional category with a coloured side annotation; columns clustered so
# similar clusters sit together (annotation aid).
heat_panel <- marker_panel[setdiff(names(marker_panel), "miR29a_targets")] # identity markers
gene2group <- unlist(lapply(names(heat_panel),
                            function(g) setNames(rep(g, length(heat_panel[[g]])),
                                                 heat_panel[[g]])))
gene2group <- gene2group[!duplicated(names(gene2group))]
gene2group <- gene2group[names(gene2group) %in% rownames(obj[["RNA"]])]
heat_feats <- names(gene2group)

avg_hm <- AverageExpression(obj, assays = "RNA", features = heat_feats,
                            group.by = "clusters", layer = "data")$RNA
avg_hm <- log1p(avg_hm)[heat_feats, , drop = FALSE]       # log space, panel row order
# AverageExpression prefixes numeric idents with "g" (g0, g1, ...). Strip it,
# then order columns NUMERICALLY (else they sort g0,g1,g10,...,g2).
colnames(avg_hm) <- sub("^g", "", colnames(avg_hm))
num <- suppressWarnings(as.numeric(colnames(avg_hm)))
if (!any(is.na(num))) avg_hm <- avg_hm[, order(num), drop = FALSE]

group_colors <- c(
  Pan_T         = "#9C6FD6", Lineage      = "#5E60CE", Naive_Memory = "#74C2E1",
  Effector_Cyto = "#EF5350", Activation   = "#FFD54F", Exhaustion   = "#A1887F",
  Treg          = "#F06292", Lineage_TFs  = "#66BB6A", Cytokines    = "#26A69A",
  Proliferation = "#FF8A65", Purity_check = "#BDBDBD"
)
row_annot <- data.frame(Group = factor(gene2group[heat_feats],
                                       levels = names(group_colors)),
                        row.names = heat_feats)
gaps <- which(head(as.character(row_annot$Group), -1) !=
                tail(as.character(row_annot$Group), -1))

png(file.path(clust_plot_dir, "Heatmap_marker_panel.png"),
    width = 12, height = 18, units = "in", res = 300, bg = "white")
pheatmap(as.matrix(avg_hm),
         scale = "row",                 # per-gene z-score
         cluster_rows = FALSE,          # fixed functional-group order
         cluster_cols = FALSE,          # keep numeric cluster order (TRUE = group similar clusters)
         annotation_row = row_annot,
         annotation_colors = list(Group = group_colors),
         gaps_row = gaps,
         color = colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(100),
         border_color = NA, fontsize_row = 9, fontsize_col = 12,
         main = "Mean expression (row z-score) by cluster - res 0.4")
dev.off()

# ===================== Tables to share for annotation =======================
# These CSVs are the inputs to hand back for annotation help.
annot_data_dir <- file.path(project_root, "Annotation", "Annotation_tables")
dir.create(annot_data_dir, recursive = TRUE, showWarnings = FALSE)

# (a) average expression of the marker panel per cluster (linear normalized)
avg_panel <- AverageExpression(obj, assays = "RNA", features = heat_feats,
                               group.by = "clusters", layer = "data")$RNA
colnames(avg_panel) <- sub("^g", "", colnames(avg_panel))   # g0->0, clean labels
np <- suppressWarnings(as.numeric(colnames(avg_panel)))
if (!any(is.na(np))) avg_panel <- avg_panel[, order(np), drop = FALSE]
write.csv(round(as.matrix(avg_panel), 4),
          file.path(annot_data_dir, "avg_expression_marker_panel.csv"))

# (b) data-driven top markers per cluster — the key annotation table
all_markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25,
                              max.cells.per.ident = 2000)   # subsample for speed
write.csv(all_markers,
          file.path(annot_data_dir, "FindAllMarkers_res0.4.csv"), row.names = FALSE)

top20 <- all_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
  dplyr::ungroup()
write.csv(top20, file.path(annot_data_dir, "FindAllMarkers_top20_res0.4.csv"),
          row.names = FALSE)

# (c) cluster sizes split by condition (context for annotation + abundance)
clust_counts <- obj@meta.data %>%
  dplyr::count(clusters, condition) %>%
  tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0)
write.csv(clust_counts,
          file.path(annot_data_dir, "cluster_counts_by_condition.csv"),
          row.names = FALSE)

message("Annotation tables written to: ", annot_data_dir)

# ============================ Per-gene plots =================================
pal <- viridis(n = 10, option = "A")   # magma

for (g in rna.features) {
  
  # 1) violin by cluster
  ggsave(file.path(rna_vln_dir, paste0(g, "_RNA_VlnPlot.png")),
         VlnPlot2(obj, features = g, cols = "default", show.mean = TRUE) +
           ggtitle(paste("RNA |", g)),
         dpi = 300, width = 14, height = 8, bg = "white")
  
  # 2) violin split by condition (EV / Scr / miR29a)
  ggsave(file.path(rna_vln_cond_dir, paste0(g, "_RNA_VlnPlot_byCondition.png")),
         VlnPlot2(obj, features = g, cols = "default", split.by = "condition",
                  stat.method = "wilcox.test") +
           ggtitle(paste("RNA |", g, "| split by condition")),
         dpi = 300, width = 14, height = 8, bg = "white")
  
  # 3) feature plot on the integrated UMAP
  ggsave(file.path(rna_feature_dir, paste0(g, "_RNA_FeaturePlot.png")),
         FeaturePlot_scCustom(obj, reduction = umap_reduction, features = g,
                              colors_use = pal, order = TRUE),
         dpi = 300, width = 8, height = 7, bg = "white")
}

# ===================== Tentative annotation =================================
# Working labels from wetlab scientist (NK_AI_06_04_26_scRNAseq_Updates.xlsx,
# tab 'Annotation') applied to the res 0.4 cluster set.
#
# Two clusters flagged for likely exclusion, pending further investigation:
#   - Cluster 8 : cycling cells lacking canonical T-cell markers (likely
#                 non-T). High DN fraction from the lineage gate (74%) is
#                 consistent. Cannot drive a CAR-T story from a non-T cell.
#   - Cluster 13: scientist labeled NK-like / innate-like (Ncr1/Klrb1c).
#                 NOTE: FindAllMarkers shows Trdc/Trgc TCR chains as top
#                 markers, which would support a gd-T identity. Resolve
#                 before final exclusion.
#
# Wetlab notes for downstream analysis planning (script 5+):
#   - CD4 island is over-split; consolidation may be appropriate after
#     miR-29a effects are established per cluster.
#   - Lab interest is CD8-primary; pseudobulk DE should run CD8 SEPARATELY
#     from CD4 to avoid CD4 dominance masking CD8-specific effects.
#   - Cluster 8 (non-T) must not drive a metabolic-pathway story in a
#     CAR-T manuscript -- pseudobulk per compartment will isolate it.
#   - Focus DE/abundance on clusters with strongest expected miR-29a
#     frequency shifts (see mir29a_freq_expected column below).
#
# Edit `state` / `lineage` / `exclude` as the annotation is refined.
annot <- data.frame(
  cluster = as.character(0:13),
  lineage = c("CD4", "CD4", "CD8", "CD4", "CD8", "CD4", "CD4", "CD4",
              "Non-T", "CD4", "CD4", "CD4", "CD4", "Innate-like"),
  state = c(
    "Activated intermediate CD4",                  # 0
    "Differentiated activated/effector CD4",       # 1
    "Proliferative CD8 effector",                  # 2
    "IFNg-associated early activated CD4",         # 3
    "Non-cycling exhaustion-like CD8",             # 4
    "CD4 TCF1hi naive/stem-like progenitor",       # 5
    "Transitional progenitor-like CD4",            # 6
    "Early activated CD4",                         # 7
    "Cycling non-T (review)",                      # 8
    "Treg",                                        # 9
    "Transitional activated CD4",                  # 10
    "CX3CR1+KLRG1+ differentiated effector CD4",   # 11
    "Intermediate activated CD4",                  # 12
    "NK-like / innate-like (review)"               # 13
  ),
  mir29a_freq_expected = c(
    "Increased", "Increased", "Reduced", "Similar/mild down", "Reduced",
    "Slight up (unclear)", "Similar", "Slight down", "Reduced", "Similar",
    "Similar", "Slight down", "Slight up", "Slight down (unclear)"
  ),
  exclude = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
              TRUE,  FALSE, FALSE, FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

# guard: cluster set in the object must match the mapping
lv_now <- levels(obj$clusters)
if (length(setdiff(lv_now, annot$cluster)))
  warning("No annotation row for cluster(s): ",
          paste(setdiff(lv_now, annot$cluster), collapse = ", "))
if (length(setdiff(annot$cluster, lv_now)))
  warning("Annotation has cluster(s) not in object: ",
          paste(setdiff(annot$cluster, lv_now), collapse = ", "))

# map onto cells
ord <- order(as.numeric(annot$cluster))
idx <- match(as.character(obj$clusters), annot$cluster)
obj$tentative_lineage    <- factor(annot$lineage[idx], levels = unique(annot$lineage[ord]))
obj$tentative_state      <- factor(annot$state[idx],   levels = unique(annot$state[ord]))
obj$tentative_annotation <- factor(
  paste0(as.character(obj$clusters), ": ", annot$state[idx]),
  levels = paste0(annot$cluster, ": ", annot$state)[ord]
)
obj$exclude_cluster      <- annot$exclude[idx]   # TRUE for clusters 8, 13

# (1) annotation table (the requested output) -- includes wetlab fields
write.csv(annot, file.path(annot_data_dir, "tentative_annotation.csv"),
          row.names = FALSE)

# (2) UMAP coloured by state (cluster identity), with cluster numbers on plot
ggsave(file.path(clust_plot_dir, "UMAP_tentative_annotation.png"),
       DimPlot2(obj, reduction = umap_reduction, group.by = "tentative_state",
                pt.size = 0.3) +
         ggtitle("Tentative annotation (wetlab labels \u2014 confirm before use)"),
       width = 14, height = 9, dpi = 300, bg = "white")

ggsave(file.path(clust_plot_dir, "UMAP_tentative_annotation_labled.png"),
       DimPlot2(obj, reduction = umap_reduction, group.by = "tentative_state",label  = T,box=T,
                pt.size = 0.3) +
         ggtitle("Tentative annotation (wetlab labels \u2014 confirm before use)"),
       width = 14, height = 9, dpi = 300, bg = "white")

# (3) UMAP by lineage compartment (CD4 / CD8 / Non-T / Innate-like)
ggsave(file.path(clust_plot_dir, "UMAP_tentative_lineage.png"),
       DimPlot2(obj, reduction = umap_reduction, group.by = "tentative_lineage",
                pt.size = 0.3) +
         ggtitle("Tentative lineage compartments"),
       width = 10, height = 8, dpi = 300, bg = "white")

# (4) re-save object WITH tentative annotation columns (overwrites checkpoint)
qs_save(obj, file = file.path(saved_dir, "Mouse_CARTmiR29a_PreAnnotation.qs2"))
message("Tentative annotation written to ",
        file.path(annot_data_dir, "tentative_annotation.csv"))

message("Done. Inspect Annotation/Pre-Annotation/ to assign cluster identities.")
###############################################################################