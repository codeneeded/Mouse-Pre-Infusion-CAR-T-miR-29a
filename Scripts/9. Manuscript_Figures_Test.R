###############################################################################
# 10_manuscript_figures_mouse_CARTmiR29a.R
#
# Publication-quality figures for the miR-29a CAR-T scRNA-seq manuscript.
#
# Plot types generated (each in its own subfolder under Manuscript_Figures/):
#   - Cluster abundance with formal significance testing (speckle::propeller)
#   - Condition-split UMAP with cluster annotations
#   - Per-cluster volcano plots highlighting miR-29a targets + curated panels
#   - Feature plots for cluster-specific gene panels (UMAP)
#   - Violin plots per gene split by condition, grouped by cluster
#   - Cluster distribution comparison (SeuratExtend::ClusterDistrPlot)
#
# ============================================================================
# IMPORTANT METHODOLOGICAL UPDATE (June 2026):
#
# EV is now the PRIMARY control. The in vivo anti-tumor experiment was
# performed against EV, so EV is the relevant biological reference. Scr is
# the SECONDARY control, used to verify mechanism specificity. This is the
# opposite of what scripts 05 and 08 assumed -- those treat Scr as primary.
# This script (10) honours the new convention. The older DGE/Pathway/Target
# outputs still exist and are valid; only the *interpretation* shifts.
# ============================================================================
#
# Inputs:
#   saved_R_data/Mouse_CARTmiR29a_WithModuleScores.qs2
#   Differential_Expression/DGE_MAST/by_cluster/*/  (per-cluster MAST CSVs)
#   Resources/miR29a_targetscan_conserved.csv
#
# Outputs:
#   Manuscript_Figures/
#     Cluster_Abundance/          propeller stats + per-replicate plot
#     Split_UMAP/                 condition-split UMAP
#     Cluster_Distribution/       ClusterDistrPlot from SeuratExtend
#     Volcano_Plots/<contrast>/   per-cluster volcanoes
#     Feature_Plots/<panel>/      UMAP feature plots per gene panel
#     Violin_Plots/<panel>/       violin plots per gene panel
#     cluster_legend.csv          number -> name mapping for figure captions
###############################################################################

library(Seurat)
library(SeuratExtend)    # DimPlot2, VlnPlot2, ClusterDistrPlot, theme_umap_arrows
library(scCustomize)     # FeaturePlot_scCustom
library(qs2)
library(speckle)         # propeller for cluster abundance significance
library(limma)           # propeller dependency
library(dplyr); library(tidyr); library(tibble)
library(ggplot2); library(ggrepel); library(patchwork); library(viridis)
library(cowplot); library(grid); library(scales)

# ============================ Paths ==========================================
project_dir <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/10x_scRNAseq"
saved_dir   <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/saved_R_data"
res_dir     <- file.path(project_dir, "Resources")
de_root     <- file.path(project_dir, "Differential_Expression")
fig_base    <- file.path(project_dir, "Manuscript_Figures")

dir_abundance <- file.path(fig_base, "Cluster_Abundance")
dir_umap      <- file.path(fig_base, "Split_UMAP")
dir_distr     <- file.path(fig_base, "Cluster_Distribution")
dir_volcano   <- file.path(fig_base, "Volcano_Plots")
dir_feature   <- file.path(fig_base, "Feature_Plots")
dir_violin    <- file.path(fig_base, "Violin_Plots")
for (d in c(dir_abundance, dir_umap, dir_distr, dir_volcano,
            dir_feature, dir_violin))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ============================ Load object ====================================
message("Loading object...")
obj <- qs_read(file.path(saved_dir, "Mouse_CARTmiR29a_WithModuleScores.qs2"))
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")
if (!"data" %in% Layers(obj[["RNA"]]))
  obj <- NormalizeData(obj, verbose = FALSE)

# Factors + sample identifier (use DASH to match script 05 / pseudobulk convention)
obj$condition <- factor(obj$condition, levels = c("EV", "Scr", "miR29a"))
obj$replicate <- factor(obj$replicate)
obj$.sample   <- paste(obj$condition, obj$replicate, sep = "-")
obj$clusters  <- droplevels(factor(obj$clusters))
obj$tentative_state <- droplevels(factor(obj$tentative_state))

umap_reduction <- if ("umap.harmony" %in% Reductions(obj)) "umap.harmony" else
  Reductions(obj)[grepl("umap", Reductions(obj), ignore.case = TRUE)][1]

# Cluster name mapping + safe filename stem
cluster_label_map <- tapply(as.character(obj$tentative_state),
                            obj$clusters, function(x) x[1])
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
cluster_stem <- function(cl_id) {
  lbl <- cluster_label_map[as.character(cl_id)]
  if (is.na(lbl) || nchar(lbl) == 0) return(paste0("cluster_", cl_id))
  safe_name(paste0(sprintf("%02d", as.integer(as.character(cl_id))), "_", lbl))
}

# Write the number -> name legend so figure captions can reference it
cluster_legend <- data.frame(
  cluster_id = as.integer(as.character(levels(obj$clusters))),
  name       = cluster_label_map[as.character(levels(obj$clusters))],
  stringsAsFactors = FALSE
)
cluster_legend <- cluster_legend[order(cluster_legend$cluster_id), ]
write.csv(cluster_legend, file.path(fig_base, "cluster_legend.csv"), row.names = FALSE)

# Condition color palette -- reused across plots
cond_cols <- c(EV = "#E76F51", Scr = "#52B788", miR29a = "#5E60CE")

# ============================ Gene panels ====================================
# Wetlab collaborator's curated cluster-specific signatures
panel_c2 <- list(
  down_in_miR29a = c("Epas1","Gpx8","Slc2a6","Cox6a2","Pim3","Hacd1","Atf5",
                     "Adora2a","Adora2b","Irs2","Agpat4","Bhlhe40","Bhlhe41","Nr4a3"),
  up_in_miR29a   = c("Slc25a23","Ldhb","Rxra","Camkk1","Sgms1","Abca3",
                     "Atp8a2","Atp10d")
)
panel_c4 <- list(
  down_in_miR29a = c("Epas1","Gpx8","Pim3","Slc2a6","Slc16a3","Agpat4","Hacd1",
                     "Atf5","Adora2a","Bhlhe40","Nr4a3"),
  up_in_miR29a   = c("Foxo3","Cpt1a","Ldhb","Mgll","Rxra","S1pr1")
)
panel_c5 <- list(
  up_in_miR29a = c("Tcf7","S1pr1","Ppargc1b","mt-Nd6","Pde4c","Btla")
)

# General-signature panel (collaborator-provided, mouse case)
panel_general <- c("Slc16a3","Rxra","Tfam","Nrf1","Opa1","Mfn1","Mfn2",
                   "Tcf7","Entpd1","Havcr2","Ccr7","P2rx7","Tigit")

# Canonical T cell category sets (always highlighted in volcanoes)
cat_methylation <- c("Tet2","Tet3","Dnmt3a","Dnmt3b","Tdg")
cat_effector_TF <- c("Tbx21","Eomes")
cat_memory_TF   <- c("Bach2","Foxo3","Tcf7","Lef1","Myb")
cat_exhaustion  <- c("Tox","Pdcd1","Havcr2","Lag3","Tigit","Entpd1","Ctla4")
cat_cytotox     <- c("Gzmb","Prf1","Nkg7","Gzma","Gzmk","Ifng","Gzmh")

# Per-cluster wetlab panel (flat vector for volcano highlighting)
cluster_specific_panel <- function(cl_id) {
  switch(as.character(cl_id),
         "2" = c(panel_c2$down_in_miR29a, panel_c2$up_in_miR29a),
         "4" = c(panel_c4$down_in_miR29a, panel_c4$up_in_miR29a),
         "5" = c(panel_c5$up_in_miR29a),
         character(0))
}

# Load TargetScan top-200 list for "miR-29a target" volcano coloring
targets_csv <- file.path(res_dir, "miR29a_targetscan_conserved.csv")
mir29a_targets_top200 <- character(0)
if (file.exists(targets_csv)) {
  ts <- read.csv(targets_csv, check.names = FALSE)
  ts <- ts[order(ts$cum_weighted_context_score), ]
  mir29a_targets_top200 <- head(ts$gene, 200)
  message("Loaded ", length(mir29a_targets_top200), " miR-29a top-200 targets.")
}


# ============================================================================
# SECTION 1: CLUSTER ABUNDANCE -- propeller + replicate-aware plot
# ============================================================================
message("\n=== Section 1: Cluster abundance ===")

run_propeller_pair <- function(seu, levels_vec) {
  sub <- subset(seu, subset = condition %in% levels_vec)
  sub$condition <- droplevels(factor(sub$condition, levels = levels_vec))
  res <- propeller(clusters = sub$clusters,
                   sample   = sub$.sample,
                   group    = sub$condition,
                   transform = "logit")
  res$cluster <- rownames(res)
  res$contrast <- paste(levels_vec[2], "vs", levels_vec[1], sep = "_")
  res
}

# Primary: miR29a vs EV. Secondary: miR29a vs Scr. QC: EV vs Scr.
prop_mir_ev  <- run_propeller_pair(obj, c("EV",  "miR29a"))
prop_mir_scr <- run_propeller_pair(obj, c("Scr", "miR29a"))
prop_ev_scr  <- run_propeller_pair(obj, c("Scr", "EV"))

prop_all <- bind_rows(prop_mir_ev, prop_mir_scr, prop_ev_scr)
write.csv(prop_all, file.path(dir_abundance, "propeller_results.csv"),
          row.names = FALSE)

# Per-replicate proportions (the data for the plot)
prop_df <- obj@meta.data %>%
  dplyr::count(condition, replicate, clusters, name = "n") %>%
  dplyr::group_by(condition, replicate) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cluster_id = as.integer(as.character(clusters)))

write.csv(prop_df, file.path(dir_abundance, "cluster_proportions_per_replicate.csv"),
          row.names = FALSE)

# Aggregate to mean + SD per condition
agg <- prop_df %>%
  dplyr::group_by(condition, clusters, cluster_id) %>%
  dplyr::summarise(mean_prop = mean(prop),
                   sd_prop   = sd(prop), .groups = "drop")

# Significance annotation rows -- primary contrast (miR29a vs EV) above bars
sig_primary <- prop_mir_ev %>%
  dplyr::transmute(cluster_id = as.integer(as.character(cluster)),
                   sig_lab_primary = dplyr::case_when(
                     is.na(FDR)   ~ "",
                     FDR < 0.001  ~ "***",
                     FDR < 0.01   ~ "**",
                     FDR < 0.05   ~ "*",
                     TRUE         ~ ""
                   ))
sig_secondary <- prop_mir_scr %>%
  dplyr::transmute(cluster_id = as.integer(as.character(cluster)),
                   sig_lab_secondary = dplyr::case_when(
                     is.na(FDR)   ~ "",
                     FDR < 0.001  ~ "***",
                     FDR < 0.01   ~ "**",
                     FDR < 0.05   ~ "*",
                     TRUE         ~ ""
                   ))

agg_y_max <- agg %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::summarise(y_top = max(mean_prop + sd_prop, na.rm = TRUE), .groups = "drop")

sig_annot <- agg_y_max %>%
  dplyr::left_join(sig_primary,   by = "cluster_id") %>%
  dplyr::left_join(sig_secondary, by = "cluster_id") %>%
  dplyr::mutate(
    # absolute offsets in proportion space prevent overlap at clusters with
    # very small y_top (where a percentage offset is essentially zero).
    y_primary   = y_top + pmax(y_top * 0.10, 0.008),
    y_secondary = y_top + pmax(y_top * 0.28, 0.022)
  )

p_abund <- ggplot(agg, aes(x = factor(cluster_id), y = mean_prop, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.85),
           color = "grey20", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(mean_prop - sd_prop, 0),
                    ymax = mean_prop + sd_prop),
                position = position_dodge(width = 0.85), width = 0.25,
                linewidth = 0.4) +
  geom_point(data = prop_df,
             aes(x = factor(cluster_id), y = prop, fill = condition),
             position = position_dodge(width = 0.85),
             shape = 21, size = 2, color = "grey20", stroke = 0.3) +
  geom_text(data = sig_annot,
            aes(x = factor(cluster_id), y = y_primary,
                label = sig_lab_primary),
            inherit.aes = FALSE, size = 5, fontface = "bold", color = "#5E60CE") +
  geom_text(data = sig_annot,
            aes(x = factor(cluster_id), y = y_secondary,
                label = sig_lab_secondary),
            inherit.aes = FALSE, size = 4, fontface = "bold", color = "grey50") +
  scale_fill_manual(values = cond_cols) +
  labs(title = "Cluster abundance by condition",
       subtitle = paste0("Bars: mean across replicates (n=2 each); points: replicates. ",
                         "Purple stars: miR29a vs EV (primary). Grey stars: miR29a vs Scr."),
       x = "Cluster",
       y = "Proportion of cells",
       fill = "Condition") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40", size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "right")

ggsave(file.path(dir_abundance, "cluster_proportions.png"),
       p_abund, width = 14, height = 7, dpi = 300, bg = "white")


# ============================================================================
# SECTION 2: SPLIT UMAP -- by condition, numbered cluster labels
# ============================================================================
message("\n=== Section 2: Split UMAP ===")

# DimPlot2 with split.by produces a broken layout (narrow strips, truncated
# facet labels, oversized legend). Workaround: make three separate DimPlot2
# objects with subset(obj, condition == X), then stitch with patchwork and
# collect legends. Each panel keeps its natural UMAP aspect ratio.

p_split_umap <- DimPlot2(
  obj,
  features    = "tentative_state",
  split.by    = "condition",
  reduction   = "umap.harmony",
  label       = TRUE,
  box         = TRUE,
  pt.size     = 0.6,
  index.title = "C",
  cols='default',
  theme       = theme_umap_arrows()
)
ggsave(file.path(dir_umap, "UMAP_split_by_condition.png"),
       p_split_umap, width = 22, height = 9, dpi = 300, bg = "white", limitsize = FALSE)

# Combined (un-split) version for reference, same labelling style
p_combined_umap <- DimPlot2(
  obj,
  features    = "tentative_state",
  reduction   = "umap.harmony",
  label       = TRUE,
  box         = TRUE,
  pt.size     = 0.6,
  index.title = "C",
  cols='default',
  theme       = theme_umap_arrows()
)
ggsave(file.path(dir_umap, "UMAP_combined.png"),
       p_combined_umap, width = 11, height = 9, dpi = 300, bg = "white")


# ============================================================================
# SECTION 3: CLUSTER DISTRIBUTION -- SeuratExtend ClusterDistrPlot
# ============================================================================
message("\n=== Section 3: Cluster distribution ===")

p_distr <- tryCatch(
  ClusterDistrPlot(
    origin    = obj$.sample,
    cluster   = obj$clusters,
    condition = obj$condition
  ),
  error = function(e) {
    message("ClusterDistrPlot failed: ", e$message,
            " -- falling back to ClusterDistrBar.")
    ClusterDistrBar(origin = obj$condition, cluster = obj$clusters)
  }
)
ggsave(file.path(dir_distr, "ClusterDistribution_by_condition.png"),
       p_distr, width = 14, height = 8, dpi = 300, bg = "white")


# ============================================================================
# SECTION 4: VOLCANO PLOTS -- per cluster, per contrast (MAST)
# ============================================================================
message("\n=== Section 4: Volcano plots ===")

# Category color palette for the GENERAL volcano (standard T cell categories).
# A separate highlight folder (Section 4b) emphasises the email gene lists for
# clusters 2, 4, 5 with their own colour scheme.
category_cols <- c(
  "miR-29a target"    = "#D97757",   # warm orange
  "Methylation axis"  = "#3B82F6",   # blue
  "Effector TF"       = "#A0303D",   # dark red
  "Memory / Stem TF"  = "#2D8659",   # green
  "Exhaustion"        = "#8B5CF6",   # purple
  "Cytotoxicity"      = "#EAB308",   # mustard yellow
  "Other"             = "grey80"
)

assign_category <- function(genes, mir29a_targets) {
  cat_vec <- rep("Other", length(genes))
  cat_vec[genes %in% cat_cytotox]     <- "Cytotoxicity"
  cat_vec[genes %in% cat_exhaustion]  <- "Exhaustion"
  cat_vec[genes %in% cat_memory_TF]   <- "Memory / Stem TF"
  cat_vec[genes %in% cat_effector_TF] <- "Effector TF"
  cat_vec[genes %in% cat_methylation] <- "Methylation axis"
  cat_vec[genes %in% mir29a_targets]  <- "miR-29a target"
  factor(cat_vec, levels = names(category_cols))
}

make_volcano <- function(de_df, title, out_png,
                         mir29a_targets,
                         dir_labels = c("Higher in miR-29a", "Higher in control")) {
  if (is.null(de_df) || nrow(de_df) == 0) return(invisible(NULL))
  
  d <- de_df
  d$gene           <- rownames(d)
  d$neg_log10_padj <- -log10(d$p_val_adj + 1e-300)
  d$category       <- assign_category(d$gene, mir29a_targets)
  
  # Label significant genes in highlighted categories
  d$label <- ""
  highlight_genes <- c(cat_methylation, cat_effector_TF, cat_memory_TF,
                       cat_exhaustion, cat_cytotox,
                       head(mir29a_targets, 25))
  to_label <- d$gene %in% highlight_genes & !is.na(d$p_val_adj) & d$p_val_adj < 0.05
  d$label[to_label] <- d$gene[to_label]
  
  # Drop non-finite rows BEFORE any axis-cap math -- these were the silent
  # killer for cluster 10 (huge -log10 p-values cascading into bad viewports).
  d <- d %>% dplyr::filter(is.finite(avg_log2FC), is.finite(neg_log10_padj))
  if (nrow(d) == 0) return(invisible(NULL))
  
  # Axis caps: 99.5th percentile + cushion, BUT hard-capped to prevent the
  # annotation_custom viewport from going to -y_cap*0.22 = -50+ when one
  # cluster has astronomically small p-values.
  x_vals <- abs(d$avg_log2FC)
  x_lim  <- if (length(x_vals) > 0) quantile(x_vals, 0.995, na.rm = TRUE) * 1.15 else 2
  if (!is.finite(x_lim) || x_lim < 1) x_lim <- 2
  x_lim  <- min(x_lim, 5)                                          # hard cap
  d$avg_log2FC_plot <- pmax(pmin(d$avg_log2FC, x_lim * 0.98), -x_lim * 0.98)
  
  y_vals <- d$neg_log10_padj
  y_cap  <- if (length(y_vals) > 0) quantile(y_vals, 0.995, na.rm = TRUE) * 1.10 else 30
  if (!is.finite(y_cap) || y_cap < 10) y_cap <- 20
  y_cap  <- min(y_cap, 100)                                        # hard cap
  d$neg_log10_padj_plot <- pmin(d$neg_log10_padj, y_cap)
  
  p <- ggplot(d, aes(x = avg_log2FC_plot, y = neg_log10_padj_plot,
                     color = category)) +
    geom_point(data = d %>% dplyr::filter(category == "Other"),
               size = 1.1, alpha = 0.3) +
    geom_point(data = d %>% dplyr::filter(category != "Other"),
               size = 3.5, alpha = 0.9) +
    ggrepel::geom_label_repel(
      data = d %>% dplyr::filter(label != ""),
      aes(label = label, fill = category),
      color = "white",                                       # white-on-coloured-fill -- matches the human template
      size = 5, fontface = "bold.italic", max.overlaps = 30,
      segment.size = 0.4, segment.color = "grey40",
      box.padding = 0.55, label.padding = 0.3,
      label.size = 0.4, show.legend = FALSE,
      alpha = 0.92, force = 2, force_pull = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = category_cols, name = "Gene category", drop = FALSE) +
    scale_fill_manual(values = category_cols, guide = "none") +
    labs(title = title,
         x = expression(log[2]~fold~change),
         y = expression(-log[10]~adjusted~italic(p))) +
    theme_cowplot(font_size = 14) +
    theme(plot.title       = element_text(face = "bold", size = 14),
          legend.position  = "right",
          plot.background  = element_rect(fill = "white", color = NA)) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
  
  # Direction arrows below x-axis -- EXACTLY the human-template pattern.
  # grobTree + annotation_custom(xmin=-Inf, xmax=Inf) + clip="off". With
  # y_cap hard-capped above, -y_cap * 0.22 stays in safe range.
  arrow_right <- grobTree(
    linesGrob(x = unit(c(0.52, 0.95), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
              gp = gpar(col = "#E76F51", lwd = 4)),                # salmon
    textGrob(dir_labels[1], x = unit(0.74, "npc"), y = unit(0, "npc"),
             gp = gpar(col = "#E76F51", fontsize = 12, fontface = "bold")))
  arrow_left <- grobTree(
    linesGrob(x = unit(c(0.48, 0.05), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
              gp = gpar(col = "#2A9D8F", lwd = 4)),                # teal
    textGrob(dir_labels[2], x = unit(0.26, "npc"), y = unit(0, "npc"),
             gp = gpar(col = "#2A9D8F", fontsize = 12, fontface = "bold")))
  
  p_final <- p +
    theme(plot.margin = margin(10, 10, 70, 10)) +
    annotation_custom(grob = arrow_right, xmin = -Inf, xmax = Inf,
                      ymin = -y_cap * 0.22, ymax = -y_cap * 0.09) +
    annotation_custom(grob = arrow_left, xmin = -Inf, xmax = Inf,
                      ymin = -y_cap * 0.22, ymax = -y_cap * 0.09) +
    coord_cartesian(ylim = c(0, y_cap), xlim = c(-x_lim, x_lim), clip = "off")
  
  ggsave(out_png, plot = p_final, width = 13, height = 10, dpi = 300, bg = "white")
  invisible(p_final)
}

contrast_pair <- c("miR29a_vs_EV", "miR29a_vs_Scr")

for (ct in contrast_pair) {
  ct_dir <- file.path(dir_volcano, ct)
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
  
  dir_labs <- if (ct == "miR29a_vs_EV") {
    c("Higher in miR-29a", "Higher in EV")
  } else {
    c("Higher in miR-29a", "Higher in Scramble")
  }
  
  for (cl in levels(obj$clusters)) {
    stem <- cluster_stem(cl)
    csv_path <- file.path(de_root, "DGE_MAST", "by_cluster", stem,
                          paste0(stem, "_", ct, ".csv"))
    if (!file.exists(csv_path)) {
      message("  missing: ", csv_path); next
    }
    
    de <- tryCatch(read.csv(csv_path, row.names = 1, check.names = FALSE),
                   error = function(e) NULL)
    if (is.null(de) || nrow(de) == 0) next
    
    cluster_label <- cluster_label_map[as.character(cl)]
    title <- sprintf("Cluster %s: %s  |  %s",
                     cl, cluster_label, gsub("_", " ", ct))
    
    make_volcano(
      de_df          = de,
      title          = title,
      out_png        = file.path(ct_dir, paste0(stem, ".png")),
      mir29a_targets = mir29a_targets_top200,
      dir_labels     = dir_labs
    )
  }
  message("  ", ct, ": volcanoes done")
}


# ============================================================================
# SECTION 4b: HIGHLIGHT VOLCANOES -- clusters 2, 4, 5 with email gene panels
# ============================================================================
# Per the wetlab collaborator, clusters 2 (Proliferative CD8), 4 (Non-cycling
# exhaustion-like CD8) and 5 (CD4 TCF1hi stem-like) have curated gene panels.
# These volcanoes label the panel genes regardless of significance (forced
# labels) with distinct down-vs-up colours so the metabolic / lipid / redox
# story reads off the plot directly.

dir_volcano_hl <- file.path(fig_base, "Volcano_Plots_Highlight")
dir.create(dir_volcano_hl, recursive = TRUE, showWarnings = FALSE)

panel_function_cols <- c(
  "Glycolysis / Hypoxia" = "#0072B2",   # deep blue   -- HIF axis, glucose import, lactate export
  "OXPHOS / FAO"         = "#D55E00",   # vermillion  -- mito energy, fatty-acid oxidation
  "Lipid metabolism"     = "#CC79A7",   # pink-purple -- lipid synthesis, membrane, peroxidation
  "Transcription factor" = "#009E73",   # green       -- TFs (Bhlhe40/41, Atf5, Nr4a3, Foxo3, Rxra, Tcf7)
  "Signaling / immune"   = "#E69F00",   # orange      -- adenosine, S1P, BTLA, kinases, immune
  "Other"                = "grey80"
)

# Per-gene function map -- covers all 32 unique genes across the c2/c4/c5
# panels. Categories chosen to make the metabolic-reprogramming story read
# off the plot (glycolysis down, OXPHOS/FAO up, lipid remodelling, TF shifts).
panel_function_map <- c(
  # Glycolysis / Hypoxia (HIF axis)
  Epas1 = "Glycolysis / Hypoxia", Slc2a6 = "Glycolysis / Hypoxia",
  Slc16a3 = "Glycolysis / Hypoxia", Pim3 = "Glycolysis / Hypoxia",
  
  # OXPHOS / mitochondrial / FAO
  Cox6a2 = "OXPHOS / FAO", `mt-Nd6` = "OXPHOS / FAO",
  Slc25a23 = "OXPHOS / FAO", Ldhb = "OXPHOS / FAO",
  Ppargc1b = "OXPHOS / FAO", Cpt1a = "OXPHOS / FAO",
  Mgll = "OXPHOS / FAO",
  
  # Lipid metabolism / membrane / redox
  Hacd1 = "Lipid metabolism", Agpat4 = "Lipid metabolism",
  Sgms1 = "Lipid metabolism", Atp8a2 = "Lipid metabolism",
  Atp10d = "Lipid metabolism", Abca3 = "Lipid metabolism",
  Gpx8 = "Lipid metabolism",
  
  # Transcription factors
  Bhlhe40 = "Transcription factor", Bhlhe41 = "Transcription factor",
  Atf5 = "Transcription factor",   Nr4a3 = "Transcription factor",
  Foxo3 = "Transcription factor",  Rxra = "Transcription factor",
  Tcf7 = "Transcription factor",
  
  # Signaling / immune modulation
  Adora2a = "Signaling / immune", Adora2b = "Signaling / immune",
  Irs2 = "Signaling / immune",    S1pr1 = "Signaling / immune",
  Btla = "Signaling / immune",    Pde4c = "Signaling / immune",
  Camkk1 = "Signaling / immune"
)

make_highlight_volcano <- function(de_df, title, out_png,
                                   panel_genes,
                                   dir_labels = c("Higher in miR-29a",
                                                  "Higher in control")) {
  if (is.null(de_df) || nrow(de_df) == 0) return(invisible(NULL))
  
  d <- de_df
  d$gene           <- rownames(d)
  d$neg_log10_padj <- -log10(d$p_val_adj + 1e-300)
  
  # Assign by FUNCTION, not direction (direction is already on the x-axis).
  # Genes outside the panel = "Other". Genes in the panel but missing from
  # the function map (shouldn't happen given the map covers all panel genes)
  # also fall to "Other".
  d$category <- "Other"
  panel_idx <- d$gene %in% panel_genes
  d$category[panel_idx] <- ifelse(
    d$gene[panel_idx] %in% names(panel_function_map),
    panel_function_map[d$gene[panel_idx]],
    "Other"
  )
  d$category <- factor(d$category, levels = names(panel_function_cols))
  
  # Label EVERY panel gene that's actually in the DE table (no sig threshold --
  # the email panel is the curated narrative, so all panel members get labelled
  # if they appear at all). Sig threshold still controls dashed reference line.
  panel_all <- panel_genes
  d$label <- ""
  to_label <- d$gene %in% panel_all
  d$label[to_label] <- d$gene[to_label]
  
  # Drop non-finite rows first (same fix as general volcano)
  d <- d %>% dplyr::filter(is.finite(avg_log2FC), is.finite(neg_log10_padj))
  if (nrow(d) == 0) return(invisible(NULL))
  
  # Axis caps with hard ceiling to keep arrow viewport math safe
  x_vals <- abs(d$avg_log2FC)
  x_lim  <- if (length(x_vals) > 0) quantile(x_vals, 0.995, na.rm = TRUE) * 1.15 else 2
  if (!is.finite(x_lim) || x_lim < 1) x_lim <- 2
  x_lim  <- min(x_lim, 5)
  d$avg_log2FC_plot <- pmax(pmin(d$avg_log2FC, x_lim * 0.98), -x_lim * 0.98)
  
  y_vals <- d$neg_log10_padj
  y_cap  <- if (length(y_vals) > 0) quantile(y_vals, 0.995, na.rm = TRUE) * 1.10 else 30
  if (!is.finite(y_cap) || y_cap < 10) y_cap <- 20
  y_cap  <- min(y_cap, 100)
  d$neg_log10_padj_plot <- pmin(d$neg_log10_padj, y_cap)
  
  p <- ggplot(d, aes(x = avg_log2FC_plot, y = neg_log10_padj_plot,
                     color = category)) +
    geom_point(data = d %>% dplyr::filter(category == "Other"),
               size = 1.1, alpha = 0.25) +
    geom_point(data = d %>% dplyr::filter(category != "Other"),
               size = 4.5, alpha = 0.95) +
    ggrepel::geom_label_repel(
      data = d %>% dplyr::filter(label != ""),
      aes(label = label, fill = category),
      color = "white",
      size = 6, fontface = "bold.italic", max.overlaps = 60,
      segment.size = 0.5, segment.color = "grey40",
      box.padding = 0.7, label.padding = 0.4,
      label.size = 0.5, show.legend = FALSE,
      alpha = 0.95, force = 3, force_pull = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = panel_function_cols,
                       name = "Gene category", drop = FALSE) +
    scale_fill_manual(values = panel_function_cols, guide = "none") +
    labs(title = title,
         x = expression(log[2]~fold~change),
         y = expression(-log[10]~adjusted~italic(p))) +
    theme_cowplot(font_size = 15) +
    theme(plot.title       = element_text(face = "bold", size = 15),
          legend.position  = "right",
          legend.text      = element_text(size = 12),
          plot.background  = element_rect(fill = "white", color = NA)) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
  
  # Same arrow style as general volcano -- human-template pattern, safe now
  # that y_cap is capped at 100.
  arrow_right <- grobTree(
    linesGrob(x = unit(c(0.52, 0.95), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
              gp = gpar(col = "#E76F51", lwd = 5)),
    textGrob(dir_labels[1], x = unit(0.74, "npc"), y = unit(0, "npc"),
             gp = gpar(col = "#E76F51", fontsize = 13, fontface = "bold")))
  arrow_left <- grobTree(
    linesGrob(x = unit(c(0.48, 0.05), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
              gp = gpar(col = "#2A9D8F", lwd = 5)),
    textGrob(dir_labels[2], x = unit(0.26, "npc"), y = unit(0, "npc"),
             gp = gpar(col = "#2A9D8F", fontsize = 13, fontface = "bold")))
  
  p_final <- p +
    theme(plot.margin = margin(10, 10, 75, 10)) +
    annotation_custom(grob = arrow_right, xmin = -Inf, xmax = Inf,
                      ymin = -y_cap * 0.22, ymax = -y_cap * 0.09) +
    annotation_custom(grob = arrow_left, xmin = -Inf, xmax = Inf,
                      ymin = -y_cap * 0.22, ymax = -y_cap * 0.09) +
    coord_cartesian(ylim = c(0, y_cap), xlim = c(-x_lim, x_lim), clip = "off")
  
  ggsave(out_png, plot = p_final, width = 14, height = 11, dpi = 300, bg = "white")
  invisible(p_final)
}

# Cluster -> panel mapping (panel_c2 has no $up_in_miR29a key for cluster 5)
highlight_panels <- list(
  "2" = list(down = panel_c2$down_in_miR29a, up = panel_c2$up_in_miR29a),
  "4" = list(down = panel_c4$down_in_miR29a, up = panel_c4$up_in_miR29a),
  "5" = list(down = character(0),            up = panel_c5$up_in_miR29a)
)

for (ct in contrast_pair) {
  ct_dir <- file.path(dir_volcano_hl, ct)
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
  
  dir_labs <- if (ct == "miR29a_vs_EV") {
    c("Higher in miR-29a", "Higher in EV")
  } else {
    c("Higher in miR-29a", "Higher in Scramble")
  }
  
  for (cl in names(highlight_panels)) {
    if (!cl %in% levels(obj$clusters)) next
    stem <- cluster_stem(cl)
    csv_path <- file.path(de_root, "DGE_MAST", "by_cluster", stem,
                          paste0(stem, "_", ct, ".csv"))
    if (!file.exists(csv_path)) next
    
    de <- tryCatch(read.csv(csv_path, row.names = 1, check.names = FALSE),
                   error = function(e) NULL)
    if (is.null(de) || nrow(de) == 0) next
    
    cluster_label <- cluster_label_map[as.character(cl)]
    title <- sprintf("Cluster %s: %s  |  %s  |  Gene panel",
                     cl, cluster_label, gsub("_", " ", ct))
    
    make_highlight_volcano(
      de_df       = de,
      title       = title,
      out_png     = file.path(ct_dir, paste0(stem, "_HIGHLIGHT.png")),
      panel_genes = c(highlight_panels[[cl]]$down, highlight_panels[[cl]]$up),
      dir_labels  = dir_labs
    )
  }
  message("  ", ct, ": highlight volcanoes done (clusters 2, 4, 5)")
}


# ============================================================================
# SECTION 5: FEATURE PLOTS -- one panel per gene set
# ============================================================================
message("\n=== Section 5: Feature plots ===")

pal_magma <- viridis(n = 10, option = "A")

plot_feature_panel <- function(genes, panel_name, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  genes_present <- intersect(genes, rownames(obj))
  missing <- setdiff(genes, genes_present)
  if (length(missing) > 0) {
    message("  [", panel_name, "] not found in data: ",
            paste(missing, collapse = ", "))
  }
  if (length(genes_present) == 0) return(invisible(NULL))
  
  # Per-gene UMAPs
  for (g in genes_present) {
    p <- FeaturePlot_scCustom(obj, reduction = umap_reduction, features = g,
                              colors_use = pal_magma, order = TRUE)
    ggsave(file.path(out_dir, paste0(g, ".png")),
           p, width = 8, height = 7, dpi = 300, bg = "white")
  }
  
  # Composite grid (no per-panel legend; saves space and removes redundancy)
  panels <- lapply(genes_present, function(g) {
    FeaturePlot_scCustom(obj, reduction = umap_reduction, features = g,
                         colors_use = pal_magma, order = TRUE) +
      theme(legend.position = "none",
            plot.title      = element_text(face = "italic"))
  })
  ncol <- min(4, length(panels))
  grid <- wrap_plots(panels, ncol = ncol)
  ggsave(file.path(out_dir, paste0("__composite_", panel_name, ".png")),
         grid, width = 5 * ncol,
         height = 4.5 * ceiling(length(panels) / ncol),
         dpi = 300, bg = "white", limitsize = FALSE)
  invisible(NULL)
}

# Cluster-specific panels (combined down+up for the feature plot view)
plot_feature_panel(
  genes      = c(panel_c2$down_in_miR29a, panel_c2$up_in_miR29a),
  panel_name = "Cluster_02_specific",
  out_dir    = file.path(dir_feature, "Cluster_02_specific")
)
plot_feature_panel(
  genes      = c(panel_c4$down_in_miR29a, panel_c4$up_in_miR29a),
  panel_name = "Cluster_04_specific",
  out_dir    = file.path(dir_feature, "Cluster_04_specific")
)
plot_feature_panel(
  genes      = c(panel_c5$up_in_miR29a),
  panel_name = "Cluster_05_specific",
  out_dir    = file.path(dir_feature, "Cluster_05_specific")
)
plot_feature_panel(
  genes      = panel_general,
  panel_name = "General_signature",
  out_dir    = file.path(dir_feature, "General_signature")
)


# ============================================================================
# SECTION 6: VIOLIN PLOTS -- per gene, split by condition, NUMBERED clusters
# ============================================================================
message("\n=== Section 6: Violin plots ===")

# x-axis = cluster NUMBER (not name) -- prevents truncation seen in earlier
# attempts; cluster_legend.csv maps numbers to names for figure captions.

plot_violin_panel <- function(genes, panel_name, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  genes_present <- intersect(genes, rownames(obj))
  if (length(genes_present) == 0) return(invisible(NULL))
  
  for (g in genes_present) {
    p <- VlnPlot2(obj, features = g,
                  group.by   = "clusters",        # numeric => short tick labels
                  split.by   = "condition",
                  cols       = "default",
                  stat.method = "wilcox.test") +
      theme(axis.text.x  = element_text(size = 10, angle = 45,
                                        hjust = 1, vjust = 1),  # rotated => no EV/Scr/miR29a collisions
            axis.title.y = element_text(face = "italic"),
            plot.title   = element_text(face = "bold.italic"),
            plot.margin  = margin(t = 8, r = 8, b = 55, l = 20, unit = "pt"))
    ggsave(file.path(out_dir, paste0(g, ".png")),
           p, width = 18, height = 8, dpi = 300, bg = "white")
  }
  invisible(NULL)
}

plot_violin_panel(c(panel_c2$down_in_miR29a, panel_c2$up_in_miR29a),
                  "Cluster_02_specific",
                  file.path(dir_violin, "Cluster_02_specific"))
plot_violin_panel(c(panel_c4$down_in_miR29a, panel_c4$up_in_miR29a),
                  "Cluster_04_specific",
                  file.path(dir_violin, "Cluster_04_specific"))
plot_violin_panel(c(panel_c5$up_in_miR29a),
                  "Cluster_05_specific",
                  file.path(dir_violin, "Cluster_05_specific"))
plot_violin_panel(panel_general,
                  "General_signature",
                  file.path(dir_violin, "General_signature"))


# ============================================================================
# Done
# ============================================================================
message("\nDone.")
message("Outputs: ", fig_base)
message("  Cluster_Abundance/      propeller stats + per-replicate bar plot")
message("  Split_UMAP/             condition-split UMAP with C-prefix labels")
message("  Cluster_Distribution/   SeuratExtend ClusterDistr plot")
message("  Volcano_Plots/<contrast>/        general volcanoes (14 clusters x 2 contrasts)")
message("  Volcano_Plots_Highlight/<contrast>/  email-panel volcanoes (clusters 2, 4, 5)")
message("  Feature_Plots/<panel>/      UMAP feature plots per gene panel")
message("  Violin_Plots/<panel>/       violin plots per gene, numbered clusters")
message("  cluster_legend.csv          number -> name mapping for captions")
message("\nNote: this script treats miR29a_vs_EV as the PRIMARY contrast.")
message("Older scripts (05, 08) use miR29a_vs_Scr as primary -- both outputs")
message("still on disk and valid; only the manuscript narrative shifts.")
###############################################################################