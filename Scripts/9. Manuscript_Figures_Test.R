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

# ---- Shift cluster numbering by +1 for manuscript readability ----
# Internally and in plots, clusters now display as 1..N instead of 0..(N-1).
# DE CSVs from script 05 are still on disk with 0-based stems; the
# cluster_stem_lookup() helper below remaps a new (1-based) cluster ID to
# the corresponding 0-based stem so file lookups still work without
# touching the upstream pipeline.
levels(obj$clusters) <- as.character(as.integer(levels(obj$clusters)) + 1L)

umap_reduction <- if ("umap.harmony" %in% Reductions(obj)) "umap.harmony" else
  Reductions(obj)[grepl("umap", Reductions(obj), ignore.case = TRUE)][1]

# Cluster name mapping + safe filename stem (rebuilt AFTER the +1 shift so
# keys are the new 1-based IDs).
cluster_label_map <- tapply(as.character(obj$tentative_state),
                            obj$clusters, function(x) x[1])
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# OUTPUT stem -- new 1-based numbering for figures and folders
cluster_stem <- function(cl_id) {
  lbl <- cluster_label_map[as.character(cl_id)]
  if (is.na(lbl) || nchar(lbl) == 0) return(paste0("cluster_", cl_id))
  safe_name(paste0(sprintf("%02d", as.integer(as.character(cl_id))), "_", lbl))
}

# LOOKUP stem -- 0-based; matches script 05 DE CSV file naming on disk
cluster_stem_lookup <- function(cl_id) {
  lbl <- cluster_label_map[as.character(cl_id)]
  if (is.na(lbl) || nchar(lbl) == 0) return(paste0("cluster_", as.integer(cl_id) - 1L))
  safe_name(paste0(sprintf("%02d", as.integer(as.character(cl_id)) - 1L), "_", lbl))
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
# Wetlab collaborator's curated cluster-specific signatures.
# Clusters 1 + 8 added per collaborator follow-up: cluster 1 has the largest
# miR-29a abundance increase; cluster 8 chosen to test if Epas1 / hypoxia
# repression generalises beyond clusters 3/5/6.
panel_c1 <- list(
  down_in_miR29a = c("Epas1","Slc2a6","Atf3","Xdh","Batf3","Ccr5","Gzma",
                     "Ly6a","Serpina3g"),
  up_in_miR29a   = c("Acss1","Nqo1","Nme4","Clybl","Cox7a1","Tcf7","Rtkn2",
                     "Cd200","Cd200r1")
)
panel_c3 <- list(
  down_in_miR29a = c("Epas1","Gpx8","Slc2a6","Cox6a2","Pim3","Hacd1","Atf5",
                     "Adora2a","Adora2b","Irs2","Agpat4","Bhlhe40","Bhlhe41","Nr4a3"),
  up_in_miR29a   = c("Slc25a23","Ldhb","Rxra","Camkk1","Sgms1","Abca3",
                     "Atp8a2","Atp10d")
)
panel_c5 <- list(
  down_in_miR29a = c("Epas1","Gpx8","Pim3","Slc2a6","Slc16a3","Agpat4","Hacd1",
                     "Atf5","Adora2a","Bhlhe40","Nr4a3"),
  up_in_miR29a   = c("Foxo3","Cpt1a","Ldhb","Mgll","Rxra","S1pr1")
)
# Cluster 5 expanded with downregulated panel from collaborator follow-up.
panel_c6 <- list(
  down_in_miR29a = c("Epas1","Slc2a6","Atf3","Ccr5","Gzma","Gzmb","Ly6a",
                     "Ly6c1","Ly6c2","Serpina3g","Batf","Cycs"),
  up_in_miR29a   = c("Tcf7","S1pr1","Ppargc1b","mt-Nd6","Pde4c","Btla")
)
panel_c8 <- list(
  down_in_miR29a = c("Epas1","Slc2a6","Atf3","Serpina3g","Ccr5","Ly6a","Ly6c1",
                     "Ly6c2","P2ry14","Slc4a7","Tgm2","Mir155hg","Eomes"),
  up_in_miR29a   = c("S1pr1","Tcf7","Klf2","Clybl","Gpd2","Slc25a24","Faah",
                     "Hvcn1","Mt1","Smad7","Cd200r1")
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

# Hypoxia / central-metabolism priority list. When labelling top-N TargetScan
# miR-29a targets in the highlight volcanoes, any non-panel target that's in
# this list AND sig in the contrast is auto-labelled regardless of its
# p-value rank, so the metabolic narrative (the "blue box" + OXPHOS / FAO
# story) survives the 10-label cap even when other targets dominate by p.
hypoxia_metabolic_priority <- c(
  # HIF axis
  "Hif1a","Epas1","Hif3a","Arnt","Vhl","Hif1an",
  # HIF-target / hypoxia response
  "Vegfa","Vegfb","Bnip3","Bnip3l","Ca9","Plod2","Loxl2","P4ha1","P4ha2",
  # Glycolysis enzymes / glucose & lactate transporters
  "Hk1","Hk2","Hk3","Pfkfb3","Pfkfb4","Pkm","Ldha","Ldhb","Pdk1","Pdk2",
  "Slc2a1","Slc2a3","Slc2a6","Slc16a1","Slc16a3","Pim3",
  # OXPHOS / mito biogenesis / dynamics
  "Sdha","Sdhb","Sdhc","Sdhd","Cycs","Cox6a2","Cox7a1","Ndufa1","Atp5b",
  "Ppargc1a","Ppargc1b","Tfam","Nrf1","Mfn1","Mfn2","Opa1","Drp1",
  # FAO + lipid handling
  "Cpt1a","Cpt1b","Cpt2","Acadm","Acadl","Acox1","Mgll","Hadha","Hadhb"
)

# Per-cluster wetlab panel (flat vector for volcano highlighting)
cluster_specific_panel <- function(cl_id) {
  switch(as.character(cl_id),
         "1" = c(panel_c1$down_in_miR29a, panel_c1$up_in_miR29a),
         "3" = c(panel_c3$down_in_miR29a, panel_c3$up_in_miR29a),
         "5" = c(panel_c5$down_in_miR29a, panel_c5$up_in_miR29a),
         "6" = c(panel_c6$down_in_miR29a, panel_c6$up_in_miR29a),
         "8" = c(panel_c8$down_in_miR29a, panel_c8$up_in_miR29a),
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

p_abund_all <- ggplot(agg, aes(x = factor(cluster_id), y = mean_prop, fill = condition)) +
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
  labs(title = "Cluster abundance by condition (all three groups)",
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

ggsave(file.path(dir_abundance, "cluster_proportions_all_conditions.png"),
       p_abund_all, width = 14, height = 7, dpi = 300, bg = "white")

# ---- Primary version per collaborator request: EV vs miR29a only, boxplot
# without individual replicate dots. With n=2 the boxplot is degenerate
# (median = midpoint, no IQR), but matches the requested visual style.
prop_df_2cond <- prop_df %>% dplyr::filter(condition %in% c("EV", "miR29a")) %>%
  dplyr::mutate(condition = factor(condition, levels = c("EV", "miR29a")))

sig_annot_2cond <- sig_annot %>%
  dplyr::transmute(cluster_id,
                   y_pos = y_primary,
                   sig_lab = sig_lab_primary)

cond_cols_2 <- cond_cols[c("EV", "miR29a")]

p_abund_2cond <- ggplot(prop_df_2cond,
                        aes(x = factor(cluster_id), y = prop, fill = condition)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7, color = "grey20", linewidth = 0.4,
               outlier.shape = NA) +
  geom_text(data = sig_annot_2cond,
            aes(x = factor(cluster_id), y = y_pos, label = sig_lab),
            inherit.aes = FALSE, size = 5, fontface = "bold", color = "#5E60CE") +
  scale_fill_manual(values = cond_cols_2) +
  labs(title = "Cluster abundance: EV vs miR-29a",
       subtitle = paste0("Boxplots over biological replicates (n=2 each). ",
                         "Purple stars: propeller FDR (miR29a vs EV)."),
       x = "Cluster",
       y = "Proportion of cells",
       fill = "Condition") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40", size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "right")

ggsave(file.path(dir_abundance, "cluster_proportions_EV_vs_miR29a.png"),
       p_abund_2cond, width = 14, height = 7, dpi = 300, bg = "white")


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
# clusters 3, 5, 6 with their own colour scheme.
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
  # killer for cluster 11 (huge -log10 p-values cascading into bad viewports).
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
  # Salmon up / deep ocean blue down; thicker (lwd=6) for high-DPI clarity.
  arrow_right <- grobTree(
    linesGrob(x = unit(c(0.52, 0.95), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
              gp = gpar(col = "#E76F51", lwd = 6)),                # salmon
    textGrob(dir_labels[1], x = unit(0.74, "npc"), y = unit(0, "npc"),
             gp = gpar(col = "#E76F51", fontsize = 12, fontface = "bold")))
  arrow_left <- grobTree(
    linesGrob(x = unit(c(0.48, 0.05), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
              gp = gpar(col = "#1B4965", lwd = 6)),                # deep ocean blue
    textGrob(dir_labels[2], x = unit(0.26, "npc"), y = unit(0, "npc"),
             gp = gpar(col = "#1B4965", fontsize = 12, fontface = "bold")))
  
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
    stem_in  <- cluster_stem_lookup(cl)   # 0-based: matches script 05 file names
    stem_out <- cluster_stem(cl)          # 1-based: new manuscript naming
    csv_path <- file.path(de_root, "DGE_MAST", "by_cluster", stem_in,
                          paste0(stem_in, "_", ct, ".csv"))
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
      out_png        = file.path(ct_dir, paste0(stem_out, ".png")),
      mir29a_targets = mir29a_targets_top200,
      dir_labels     = dir_labs
    )
  }
  message("  ", ct, ": volcanoes done")
}


# ============================================================================
# SECTION 4b: HIGHLIGHT VOLCANOES -- clusters 1, 3, 5, 6, 8 with email gene panels
# ============================================================================
# Per the wetlab collaborator, the original highlights were clusters 3
# (Proliferative CD8), 5 (Non-cycling exhaustion-like CD8) and 6 (CD4 TCF1hi
# stem-like). Cluster 1 (largest miR-29a abundance increase) and cluster 8
# were added in the follow-up so we can test if hypoxia/Epas1 repression
# generalises. Volcanoes label panel genes regardless of significance
# (forced labels) coloured by function so the metabolic / cytotoxic /
# lipid / redox story reads off the plot directly.

dir_volcano_hl <- file.path(fig_base, "Volcano_Plots_Highlight")
dir.create(dir_volcano_hl, recursive = TRUE, showWarnings = FALSE)

panel_function_cols <- c(
  "Glycolysis / Hypoxia"     = "#0072B2",   # deep blue   -- HIF axis, glucose import, lactate export
  "OXPHOS / FAO"             = "#D55E00",   # vermillion  -- mito energy, fatty-acid oxidation, redox-active mito enzymes
  "Lipid / Redox"            = "#CC79A7",   # pink-purple -- lipid synth, membrane, peroxidation, antioxidant defence
  "miR-29a x metabolic"      = "#F50057",   # vivid magenta -- INTERSECTION: gene is a TargetScan miR-29a target AND in a metabolic category. Overrides the metabolic colour so direct miR-29a-driven metabolic hits are trackable separately.
  "Transcription factor"     = "#009E73",   # green       -- TFs (effector, memory, stress, activation)
  "Cytotoxicity / Effector"  = "#6A3D9A",   # purple      -- Gzma/Gzmb, Ly6 activation cluster, chemokine receptors, IFN-induced
  "Signaling / immune"       = "#E69F00",   # orange      -- adenosine, S1P, BTLA, kinases, immune modulation
  "miR-29a target"           = "#8B4513",   # saddle brown-- non-metabolic TargetScan-predicted miR-29a targets (overlay)
  "Other"                    = "grey80"
)

# Per-gene function map. Now BROADER than the wetlab panel set: any canonical
# member of these categories gets coloured by function regardless of whether
# it's in the cluster's email panel. This is what makes e.g. Gzma/Gzmb appear
# purple in cluster 3 (Proliferative CD8) even though they're not in that
# cluster's panel -- previously the Cytotoxicity legend entry was blank.
panel_function_map <- c(
  # Glycolysis / Hypoxia (HIF axis, glucose import, lactate export)
  Epas1 = "Glycolysis / Hypoxia", Hif1a = "Glycolysis / Hypoxia",
  Hif3a = "Glycolysis / Hypoxia", Arnt = "Glycolysis / Hypoxia",
  Vhl = "Glycolysis / Hypoxia",   Hif1an = "Glycolysis / Hypoxia",
  Vegfa = "Glycolysis / Hypoxia", Vegfb = "Glycolysis / Hypoxia",
  Bnip3 = "Glycolysis / Hypoxia", Bnip3l = "Glycolysis / Hypoxia",
  Ca9 = "Glycolysis / Hypoxia",   Plod2 = "Glycolysis / Hypoxia",
  Loxl2 = "Glycolysis / Hypoxia", P4ha1 = "Glycolysis / Hypoxia",
  P4ha2 = "Glycolysis / Hypoxia",
  Hk1 = "Glycolysis / Hypoxia",   Hk2 = "Glycolysis / Hypoxia",
  Hk3 = "Glycolysis / Hypoxia",   Pfkfb3 = "Glycolysis / Hypoxia",
  Pfkfb4 = "Glycolysis / Hypoxia", Pkm = "Glycolysis / Hypoxia",
  Ldha = "Glycolysis / Hypoxia", Pdk1 = "Glycolysis / Hypoxia",
  Pdk2 = "Glycolysis / Hypoxia",
  Slc2a1 = "Glycolysis / Hypoxia", Slc2a3 = "Glycolysis / Hypoxia",
  Slc2a6 = "Glycolysis / Hypoxia", Slc16a1 = "Glycolysis / Hypoxia",
  Slc16a3 = "Glycolysis / Hypoxia", Pim3 = "Glycolysis / Hypoxia",
  
  # OXPHOS / mitochondrial / FAO
  Cox6a2 = "OXPHOS / FAO", Cox7a1 = "OXPHOS / FAO",
  `mt-Nd6` = "OXPHOS / FAO", Ndufa1 = "OXPHOS / FAO",
  Ndufa6 = "OXPHOS / FAO", Atp5b = "OXPHOS / FAO",
  Slc25a23 = "OXPHOS / FAO", Slc25a24 = "OXPHOS / FAO",
  Ldhb = "OXPHOS / FAO",    Acss1 = "OXPHOS / FAO",
  Nme4 = "OXPHOS / FAO",    Clybl = "OXPHOS / FAO",
  Cycs = "OXPHOS / FAO",    Gpd2 = "OXPHOS / FAO",
  Xdh = "OXPHOS / FAO",
  Sdha = "OXPHOS / FAO",    Sdhb = "OXPHOS / FAO",
  Sdhc = "OXPHOS / FAO",    Sdhd = "OXPHOS / FAO",
  Ppargc1a = "OXPHOS / FAO", Ppargc1b = "OXPHOS / FAO",
  Tfam = "OXPHOS / FAO",    Nrf1 = "OXPHOS / FAO",
  Mfn1 = "OXPHOS / FAO",    Mfn2 = "OXPHOS / FAO",
  Opa1 = "OXPHOS / FAO",    Drp1 = "OXPHOS / FAO",
  Cpt1a = "OXPHOS / FAO",   Cpt1b = "OXPHOS / FAO",
  Cpt2 = "OXPHOS / FAO",    Acadm = "OXPHOS / FAO",
  Acadl = "OXPHOS / FAO",   Acox1 = "OXPHOS / FAO",
  Hadha = "OXPHOS / FAO",   Hadhb = "OXPHOS / FAO",
  Mgll = "OXPHOS / FAO",
  
  # Lipid metabolism / membrane / redox defence
  Hacd1 = "Lipid / Redox",  Agpat4 = "Lipid / Redox",
  Sgms1 = "Lipid / Redox",  Atp8a2 = "Lipid / Redox",
  Atp10d = "Lipid / Redox", Abca3 = "Lipid / Redox",
  Gpx8 = "Lipid / Redox",   Faah = "Lipid / Redox",
  Nqo1 = "Lipid / Redox",   Mt1 = "Lipid / Redox",
  Gpx1 = "Lipid / Redox",   Gpx4 = "Lipid / Redox",
  Sod1 = "Lipid / Redox",   Sod2 = "Lipid / Redox",
  Cat = "Lipid / Redox",    Txn1 = "Lipid / Redox",
  Prdx1 = "Lipid / Redox",
  
  # Transcription factors (effector, memory, stress, activation, exhaustion)
  Bhlhe40 = "Transcription factor", Bhlhe41 = "Transcription factor",
  Atf5 = "Transcription factor",    Nr4a3 = "Transcription factor",
  Foxo3 = "Transcription factor",   Rxra = "Transcription factor",
  Tcf7 = "Transcription factor",    Atf3 = "Transcription factor",
  Batf = "Transcription factor",    Batf3 = "Transcription factor",
  Eomes = "Transcription factor",   Klf2 = "Transcription factor",
  Tbx21 = "Transcription factor",   Bach2 = "Transcription factor",
  Lef1 = "Transcription factor",    Myb = "Transcription factor",
  Tox = "Transcription factor",
  
  # Cytotoxicity / Effector / activation-stress program
  Gzma = "Cytotoxicity / Effector", Gzmb = "Cytotoxicity / Effector",
  Gzmk = "Cytotoxicity / Effector", Gzmh = "Cytotoxicity / Effector",
  Prf1 = "Cytotoxicity / Effector", Nkg7 = "Cytotoxicity / Effector",
  Ifng = "Cytotoxicity / Effector",
  Ccr5 = "Cytotoxicity / Effector", Ly6a = "Cytotoxicity / Effector",
  Ly6c1 = "Cytotoxicity / Effector", Ly6c2 = "Cytotoxicity / Effector",
  Serpina3g = "Cytotoxicity / Effector", Mir155hg = "Cytotoxicity / Effector",
  
  # Signaling / immune modulation / checkpoints
  Adora2a = "Signaling / immune", Adora2b = "Signaling / immune",
  Irs2 = "Signaling / immune",    S1pr1 = "Signaling / immune",
  Btla = "Signaling / immune",    Pde4c = "Signaling / immune",
  Camkk1 = "Signaling / immune",  Cd200 = "Signaling / immune",
  Cd200r1 = "Signaling / immune", Rtkn2 = "Signaling / immune",
  P2ry14 = "Signaling / immune",  Slc4a7 = "Signaling / immune",
  Tgm2 = "Signaling / immune",    Hvcn1 = "Signaling / immune",
  Smad7 = "Signaling / immune",   Pdcd1 = "Signaling / immune",
  Havcr2 = "Signaling / immune",  Lag3 = "Signaling / immune",
  Tigit = "Signaling / immune",   Entpd1 = "Signaling / immune",
  Ctla4 = "Signaling / immune",   Icos = "Signaling / immune"
)

make_highlight_volcano <- function(de_df, title, out_png,
                                   panel_genes,
                                   mir29a_targets  = character(0),
                                   n_target_labels = 10,
                                   priority_genes  = hypoxia_metabolic_priority,
                                   dir_labels = c("Higher in miR-29a",
                                                  "Higher in control")) {
  if (is.null(de_df) || nrow(de_df) == 0) return(invisible(NULL))
  
  d <- de_df
  d$gene           <- rownames(d)
  d$neg_log10_padj <- -log10(d$p_val_adj + 1e-300)
  
  # Category priority (highest wins last assignment):
  #   1. miR-29a target overlay applies first (everything in mir29a_targets)
  #   2. Function map overrides -- ANY gene mapped to a function category
  #      gets its category colour, regardless of panel membership.
  #   3. INTERSECTION: gene is BOTH a miR-29a target AND in a metabolic
  #      function category -> the "miR-29a x metabolic" combined category.
  #      Overrides everything so the highest-value hits (direct miR-29a
  #      targets that drive metabolic reprogramming) are trackable.
  d$category <- "Other"
  d$category[d$gene %in% mir29a_targets] <- "miR-29a target"
  in_map <- d$gene %in% names(panel_function_map)
  d$category[in_map] <- panel_function_map[d$gene[in_map]]
  
  metabolic_cats <- c("Glycolysis / Hypoxia", "OXPHOS / FAO", "Lipid / Redox")
  gene_func <- panel_function_map[d$gene]
  is_target_metabolic <- d$gene %in% mir29a_targets &
    !is.na(gene_func) &
    gene_func %in% metabolic_cats
  d$category[is_target_metabolic] <- "miR-29a x metabolic"
  
  d$category <- factor(d$category, levels = names(panel_function_cols))
  
  # Labels: every panel gene present in DE table (no sig threshold), plus
  # every gene in the miR-29a x metabolic intersection that's sig (no cap --
  # these are the manuscript's headline genes), plus up to n_target_labels
  # non-panel miR-29a TargetScan hits. Hypoxia / metabolic priority genes
  # are auto-labelled first so the "blue box" metabolic story survives the
  # cap even if other targets dominate by p-value alone.
  d$label <- ""
  d$label[d$gene %in% panel_genes] <- d$gene[d$gene %in% panel_genes]
  
  # ALWAYS label sig genes in the miR-29a x metabolic intersection
  intersection_genes <- d %>%
    dplyr::filter(category == "miR-29a x metabolic",
                  !is.na(p_val_adj), p_val_adj < 0.05) %>%
    dplyr::pull(gene)
  d$label[d$gene %in% intersection_genes] <-
    d$gene[d$gene %in% intersection_genes]
  
  if (length(mir29a_targets) > 0 && n_target_labels > 0) {
    candidates <- d %>%
      dplyr::filter(gene %in% mir29a_targets,
                    !gene %in% panel_genes,
                    category != "miR-29a x metabolic",   # already auto-labelled above
                    !is.na(p_val_adj), p_val_adj < 0.05) %>%
      dplyr::arrange(dplyr::desc(abs(avg_log2FC)))   # effect size first, p only gates
    
    # Pass 1: hypoxia/metabolic priority hits among the sig candidates
    priority_hits <- candidates %>%
      dplyr::filter(gene %in% priority_genes)
    
    # Pass 2: fill remaining slots with top |log2FC| from non-priority candidates
    remaining_slots <- max(0, n_target_labels - nrow(priority_hits))
    other_hits <- candidates %>%
      dplyr::filter(!gene %in% priority_hits$gene) %>%
      dplyr::slice_head(n = remaining_slots)
    
    # Cap total at n_target_labels (if priority alone exceeds cap, keep top
    # n_target_labels by |log2FC| among the priority hits).
    target_to_label <- dplyr::bind_rows(priority_hits, other_hits) %>%
      dplyr::slice_head(n = n_target_labels)
    
    d$label[d$gene %in% target_to_label$gene] <-
      d$gene[d$gene %in% target_to_label$gene]
  }
  
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
  
  # Ghost row per category at NA coords -- makes the legend always render a
  # swatch for every level in panel_function_cols, even when no real data
  # falls into that category (fixes the "Cytotoxicity / Effector" entry
  # showing text-only with no colour dot).
  ghost_df <- data.frame(
    avg_log2FC_plot      = NA_real_,
    neg_log10_padj_plot  = NA_real_,
    category             = factor(names(panel_function_cols),
                                  levels = names(panel_function_cols))
  )
  
  p <- ggplot(d, aes(x = avg_log2FC_plot, y = neg_log10_padj_plot,
                     color = category)) +
    geom_point(data = ghost_df,
               aes(color = category),
               size = 4.5, alpha = 1, na.rm = TRUE) +
    geom_point(data = d %>% dplyr::filter(category == "Other"),
               size = 1.1, alpha = 0.25) +
    geom_point(data = d %>% dplyr::filter(category == "miR-29a target"),
               size = 2.8, alpha = 0.7) +
    geom_point(data = d %>% dplyr::filter(!category %in% c("Other", "miR-29a target")),
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
  
  # Arrow style: salmon up / deep ocean blue down. Thickness bumped so the
  # arrows read clearly under high-DPI rendering.
  arrow_right <- grobTree(
    linesGrob(x = unit(c(0.52, 0.95), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(0.45, "cm"), type = "closed"),
              gp = gpar(col = "#E76F51", lwd = 7)),
    textGrob(dir_labels[1], x = unit(0.74, "npc"), y = unit(0, "npc"),
             gp = gpar(col = "#E76F51", fontsize = 13, fontface = "bold")))
  arrow_left <- grobTree(
    linesGrob(x = unit(c(0.48, 0.05), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(0.45, "cm"), type = "closed"),
              gp = gpar(col = "#1B4965", lwd = 7)),
    textGrob(dir_labels[2], x = unit(0.26, "npc"), y = unit(0, "npc"),
             gp = gpar(col = "#1B4965", fontsize = 13, fontface = "bold")))
  
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

# Cluster -> panel mapping. Clusters 1 + 8 added per collaborator follow-up;
# cluster 6 now also has a down panel.
highlight_panels <- list(
  "1" = list(down = panel_c1$down_in_miR29a, up = panel_c1$up_in_miR29a),
  "3" = list(down = panel_c3$down_in_miR29a, up = panel_c3$up_in_miR29a),
  "5" = list(down = panel_c5$down_in_miR29a, up = panel_c5$up_in_miR29a),
  "6" = list(down = panel_c6$down_in_miR29a, up = panel_c6$up_in_miR29a),
  "8" = list(down = panel_c8$down_in_miR29a, up = panel_c8$up_in_miR29a)
)

# Highlight loop runs on three contrasts. EV_vs_Scr is a negative control:
# miR-29a targets should be flat between two non-miR29a conditions.
highlight_contrasts <- c("miR29a_vs_EV", "miR29a_vs_Scr", "EV_vs_Scr")

for (ct in highlight_contrasts) {
  ct_dir <- file.path(dir_volcano_hl, ct)
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
  
  dir_labs <- switch(ct,
                     "miR29a_vs_EV"  = c("Higher in miR-29a", "Higher in EV"),
                     "miR29a_vs_Scr" = c("Higher in miR-29a", "Higher in Scramble"),
                     "EV_vs_Scr"     = c("Higher in EV",      "Higher in Scramble")
  )
  
  for (cl in names(highlight_panels)) {
    if (!cl %in% levels(obj$clusters)) next
    stem_in  <- cluster_stem_lookup(cl)   # 0-based for DE CSV lookup
    stem_out <- cluster_stem(cl)          # 1-based for output PNG
    csv_path <- file.path(de_root, "DGE_MAST", "by_cluster", stem_in,
                          paste0(stem_in, "_", ct, ".csv"))
    if (!file.exists(csv_path)) next
    
    de <- tryCatch(read.csv(csv_path, row.names = 1, check.names = FALSE),
                   error = function(e) NULL)
    if (is.null(de) || nrow(de) == 0) next
    
    cluster_label <- cluster_label_map[as.character(cl)]
    title <- sprintf("Cluster %s: %s  |  %s  |  Gene panel",
                     cl, cluster_label, gsub("_", " ", ct))
    
    make_highlight_volcano(
      de_df          = de,
      title          = title,
      out_png        = file.path(ct_dir, paste0(stem_out, "_HIGHLIGHT.png")),
      panel_genes    = c(highlight_panels[[cl]]$down, highlight_panels[[cl]]$up),
      mir29a_targets = mir29a_targets_top200,
      dir_labels     = dir_labs
    )
  }
  message("  ", ct, ": highlight volcanoes done (clusters 1, 3, 5, 6, 8)")
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
  genes      = c(panel_c1$down_in_miR29a, panel_c1$up_in_miR29a),
  panel_name = "Cluster_01_specific",
  out_dir    = file.path(dir_feature, "Cluster_01_specific")
)
plot_feature_panel(
  genes      = c(panel_c3$down_in_miR29a, panel_c3$up_in_miR29a),
  panel_name = "Cluster_03_specific",
  out_dir    = file.path(dir_feature, "Cluster_03_specific")
)
plot_feature_panel(
  genes      = c(panel_c5$down_in_miR29a, panel_c5$up_in_miR29a),
  panel_name = "Cluster_05_specific",
  out_dir    = file.path(dir_feature, "Cluster_05_specific")
)
plot_feature_panel(
  genes      = c(panel_c6$down_in_miR29a, panel_c6$up_in_miR29a),
  panel_name = "Cluster_06_specific",
  out_dir    = file.path(dir_feature, "Cluster_06_specific")
)
plot_feature_panel(
  genes      = c(panel_c8$down_in_miR29a, panel_c8$up_in_miR29a),
  panel_name = "Cluster_08_specific",
  out_dir    = file.path(dir_feature, "Cluster_08_specific")
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

plot_violin_panel <- function(genes, panel_name, out_dir,
                              conditions = NULL) {
  # If conditions is provided, subset obj to those conditions only (drops
  # unused factor levels). Used for the EV-vs-miR-29a-only violins.
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  genes_present <- intersect(genes, rownames(obj))
  if (length(genes_present) == 0) return(invisible(NULL))
  
  obj_use <- obj
  if (!is.null(conditions)) {
    keep <- obj_use$condition %in% conditions
    obj_use <- obj_use[, keep]
    obj_use$condition <- factor(as.character(obj_use$condition), levels = conditions)
  }
  
  for (g in genes_present) {
    p <- VlnPlot2(obj_use, features = g,
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

# Cluster-specific panels (all three conditions)
plot_violin_panel(c(panel_c1$down_in_miR29a, panel_c1$up_in_miR29a),
                  "Cluster_01_specific",
                  file.path(dir_violin, "Cluster_01_specific"))
plot_violin_panel(c(panel_c3$down_in_miR29a, panel_c3$up_in_miR29a),
                  "Cluster_03_specific",
                  file.path(dir_violin, "Cluster_03_specific"))
plot_violin_panel(c(panel_c5$down_in_miR29a, panel_c5$up_in_miR29a),
                  "Cluster_05_specific",
                  file.path(dir_violin, "Cluster_05_specific"))
plot_violin_panel(c(panel_c6$down_in_miR29a, panel_c6$up_in_miR29a),
                  "Cluster_06_specific",
                  file.path(dir_violin, "Cluster_06_specific"))
plot_violin_panel(c(panel_c8$down_in_miR29a, panel_c8$up_in_miR29a),
                  "Cluster_08_specific",
                  file.path(dir_violin, "Cluster_08_specific"))
plot_violin_panel(panel_general,
                  "General_signature",
                  file.path(dir_violin, "General_signature"))

# Top miR-29a TargetScan targets, EV vs miR-29a only -- per collaborator
# follow-up. n=top 40 to keep file count manageable; can extend if needed.
message("  miR-29a target violins (EV vs miR-29a only)...")
plot_violin_panel(
  genes      = head(mir29a_targets_top200, 40),
  panel_name = "miR29a_targets_EV_vs_miR29a",
  out_dir    = file.path(dir_violin, "miR29a_targets_EV_vs_miR29a"),
  conditions = c("EV", "miR29a")
)


# ============================================================================
# SECTION 7: HYPOXIA PATHWAY ENRICHMENT across clusters (miR29a vs EV)
# ============================================================================
# Pulls EnrichR pathway results written by script 06 and extracts hypoxia /
# HIF gene sets across clusters for the primary contrast (miR29a_vs_EV).
# If hypoxia is a true pathway-level signal, it should be enriched on the
# DOWN-in-miR-29a side across many clusters (consistent with Epas1 / Slc2a6
# repression).
message("\n=== Section 7: Hypoxia pathway enrichment ===")

dir_hypoxia <- file.path(fig_base, "Hypoxia_Enrichment")
dir.create(dir_hypoxia, recursive = TRUE, showWarnings = FALSE)

# Path to EnrichR outputs written by script 06. Adjust if directory layout
# differs in your project (the structure here matches the documented script
# 06 convention: Pathway_Analysis_EnrichR/<level>/<stem>/<contrast>/<db>.csv).
enrichr_root <- file.path(project_dir, "Pathway_Analysis_EnrichR")

# Pathway databases that contain hypoxia/HIF gene sets we care about
hypoxia_keywords <- c("hypoxi", "hif", "oxygen", "HIF-1")    # case-insensitive grep

collect_hypoxia <- function(level = "by_cluster",
                            contrast = "miR29a_vs_EV") {
  level_dir <- file.path(enrichr_root, level)
  if (!dir.exists(level_dir)) {
    message("  No EnrichR directory found at ", level_dir, " -- skipping.")
    return(NULL)
  }
  
  stems <- list.dirs(level_dir, recursive = FALSE, full.names = FALSE)
  if (length(stems) == 0) return(NULL)
  
  rows <- list()
  for (stem in stems) {
    contrast_dir <- file.path(level_dir, stem, contrast)
    if (!dir.exists(contrast_dir)) next
    csv_files <- list.files(contrast_dir, pattern = "\\.csv$", full.names = TRUE)
    if (length(csv_files) == 0) next
    
    for (csvf in csv_files) {
      db_name <- tools::file_path_sans_ext(basename(csvf))
      en <- tryCatch(read.csv(csvf, check.names = FALSE),
                     error = function(e) NULL)
      if (is.null(en) || nrow(en) == 0) next
      # Identify the term column heuristically
      term_col <- intersect(c("Term", "term", "Description", "pathway"),
                            colnames(en))[1]
      if (is.na(term_col)) next
      pval_col <- intersect(c("Adjusted.P.value", "adj_p_value", "padj",
                              "p_adj", "P.value"), colnames(en))[1]
      if (is.na(pval_col)) next
      odds_col <- intersect(c("Odds.Ratio", "odds_ratio", "OR"),
                            colnames(en))[1]
      dir_col  <- intersect(c("direction", "Direction", "side"),
                            colnames(en))[1]
      
      hits <- en[grepl(paste(hypoxia_keywords, collapse = "|"),
                       en[[term_col]], ignore.case = TRUE), , drop = FALSE]
      if (nrow(hits) == 0) next
      
      rows[[length(rows) + 1]] <- data.frame(
        stem      = stem,
        database  = db_name,
        term      = hits[[term_col]],
        p_adj     = as.numeric(hits[[pval_col]]),
        odds      = if (!is.na(odds_col)) as.numeric(hits[[odds_col]]) else NA_real_,
        direction = if (!is.na(dir_col))  as.character(hits[[dir_col]])  else NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(rows) == 0) return(NULL)
  dplyr::bind_rows(rows)
}

hyp_cluster <- collect_hypoxia("by_cluster",  "miR29a_vs_EV")
hyp_lineage <- collect_hypoxia("by_lineage",  "miR29a_vs_EV")

if (!is.null(hyp_cluster) && nrow(hyp_cluster) > 0) {
  write.csv(hyp_cluster, file.path(dir_hypoxia, "hypoxia_terms_by_cluster.csv"),
            row.names = FALSE)
  
  # Dot plot: cluster x term, fill = -log10 padj, size = odds ratio
  hp <- hyp_cluster %>%
    dplyr::mutate(neg_log10_padj = -log10(pmax(p_adj, 1e-50)),
                  cluster_id = as.integer(sub("^(\\d+)_.*$", "\\1", stem)) + 1L)
  
  hp_top_terms <- hp %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(min_padj = min(p_adj, na.rm = TRUE),
                     n_sig    = sum(p_adj < 0.05, na.rm = TRUE),
                     .groups  = "drop") %>%
    dplyr::arrange(min_padj) %>%
    head(15) %>%
    dplyr::pull(term)
  
  hp_plot <- hp %>% dplyr::filter(term %in% hp_top_terms)
  
  if (nrow(hp_plot) > 0) {
    p_hyp <- ggplot(hp_plot,
                    aes(x = factor(as.integer(cluster_id)),
                        y = term,
                        size = pmin(odds, 20),
                        fill = neg_log10_padj)) +
      geom_point(shape = 21, color = "grey20", stroke = 0.3) +
      scale_fill_gradientn(colors = c("#FFF7BC", "#FEC44F", "#D95F0E", "#7F2704"),
                           name = expression(-log[10]~p[BH])) +
      scale_size_continuous(range = c(2, 10), name = "Odds ratio") +
      labs(title = "Hypoxia / HIF pathway enrichment, miR29a vs EV, per cluster",
           subtitle = paste0("Terms grep-matching hypoxi/HIF/oxygen across all ",
                             "Pathway databases (top 15 by min p_adj). ",
                             "Down-in-miR-29a side is the biologically expected direction."),
           x = "Cluster",
           y = NULL) +
      theme_minimal(base_size = 11) +
      theme(plot.title    = element_text(face = "bold"),
            plot.subtitle = element_text(color = "grey40", size = 9),
            axis.text.y   = element_text(size = 9),
            panel.grid.minor = element_blank(),
            legend.position = "right")
    
    ggsave(file.path(dir_hypoxia, "hypoxia_enrichment_by_cluster.png"),
           p_hyp, width = 14, height = 7, dpi = 300, bg = "white")
  }
} else {
  message("  No hypoxia terms found in by_cluster EnrichR outputs.")
}

if (!is.null(hyp_lineage) && nrow(hyp_lineage) > 0) {
  write.csv(hyp_lineage, file.path(dir_hypoxia, "hypoxia_terms_by_lineage.csv"),
            row.names = FALSE)
}


# ============================================================================
# Done
# ============================================================================
message("\nDone.")
message("Outputs: ", fig_base)
message("  Cluster_Abundance/      propeller stats; bars+dots (all 3 conditions) AND boxplot (EV vs miR29a)")
message("  Split_UMAP/             condition-split UMAP with C-prefix labels")
message("  Cluster_Distribution/   SeuratExtend ClusterDistr plot")
message("  Volcano_Plots/<contrast>/        general volcanoes (14 clusters x 2 contrasts)")
message("  Volcano_Plots_Highlight/<contrast>/  function-coloured panels for clusters 1/3/5/6/8,")
message("                                     3 contrasts (incl. EV_vs_Scr control),")
message("                                     miR-29a TargetScan targets overlaid")
message("  Feature_Plots/<panel>/      UMAP feature plots per gene panel (incl. clusters 1/6/8)")
message("  Violin_Plots/<panel>/       violin plots per gene; new miR29a_targets_EV_vs_miR29a/ subfolder")
message("  Hypoxia_Enrichment/         hypoxia/HIF terms grepped from EnrichR pathway outputs")
message("  cluster_legend.csv          number -> name mapping for captions")
message("\nNote: this script treats miR29a_vs_EV as the PRIMARY contrast.")
message("Older scripts (05, 08) use miR29a_vs_Scr as primary -- both outputs")
message("still on disk and valid; only the manuscript narrative shifts.")
###############################################################################