###############################################################################
# 07 -- Module Scores & Feature Visualizations                                #
# Mouse pre-infusion CAR-T scRNA-seq, miR-29a project                         #
###############################################################################
# Implements wetlab Sheet 2 (NK_AI_06_04_26_scRNAseq_Updates.xlsx):
#   - Reads gene lists from Resources/Module_Gene_Lists.xlsx
#     (16 modules: metabolism, FOXO axis, T-cell state)
#   - Reads the conserved miR-29a-3p TargetScan list
#     (Resources/miR29a_targetscan_conserved.csv) and builds a target module
#     from the top-N strongest predicted targets.
#   - Computes AddModuleScore for every module; writes scores back as object
#     metadata so the referee (script 08) and any downstream plotting can use
#     them directly.
#
# Outputs per module:
#   - FeaturePlot on UMAP
#   - VlnPlot by tentative_state (which clusters express the module)
#   - VlnPlot by condition within each cluster (the miR-29a-effect view)
#
# A composite grid (`UMAP_priority_modules.png`) shows the 9 modules the
# wetlab explicitly requested for UMAP scoring (priority_UMAP = TRUE in the
# Excel). All other modules are still computed and plotted individually.
#
# Headline test:
#   miR-29a target module score is compared across conditions with a Wilcoxon
#   test (overall and per cluster). Lower scores in miR29a vs control =
#   direct evidence the targets are being repressed.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(readxl); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(viridis)
  library(SeuratExtend); library(scCustomize)
})

# ============================ Paths ==========================================
project_dir   <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a"
saved_dir     <- file.path(project_dir, "saved_R_data")
resources_dir <- file.path(project_dir, "Resources")
mod_xlsx      <- file.path(resources_dir, "Module_Gene_Lists.xlsx")
mir29a_csv    <- file.path(resources_dir, "miR29a_targetscan_conserved.csv")

out_base      <- file.path(project_dir, "Module_Scores")
plot_dir      <- file.path(out_base, "Plots")
umap_dir      <- file.path(plot_dir, "UMAP")
vln_clu_dir   <- file.path(plot_dir, "VlnByCluster")
vln_cnd_dir   <- file.path(plot_dir, "VlnByCondCluster")
data_dir      <- file.path(out_base, "Tables")
for (d in c(plot_dir, umap_dir, vln_clu_dir, vln_cnd_dir, data_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ============================ Load ===========================================
obj <- qs_read(file.path(saved_dir, "Mouse_CARTmiR29a_PreAnnotation.qs2"))
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")
if (!"data" %in% Layers(obj[["RNA"]]))
  obj <- NormalizeData(obj, verbose = FALSE)

# Keep ALL clusters for module-score visualization, including the wetlab-flagged
# "Cycling non-T" (8) and "NK-like / innate-like" (13). The exclude_cluster flag
# only governs DE / statistical inference (script 05); for plotting metabolic,
# FOXO, and target-score modules, every cluster is part of the landscape.
obj$clusters          <- droplevels(factor(obj$clusters))
obj$tentative_lineage <- droplevels(factor(obj$tentative_lineage))
obj$tentative_state   <- droplevels(factor(obj$tentative_state))
obj$condition         <- factor(obj$condition, levels = c("EV", "Scr", "miR29a"))
umap_reduction        <- "umap.harmony"
if (!umap_reduction %in% Reductions(obj))
  umap_reduction <- Reductions(obj)[grepl("umap", Reductions(obj), ignore.case = TRUE)][1]

# ============================ Read module lists ==============================
mods_long <- read_excel(mod_xlsx, sheet = "Module_Gene_Lists")
stopifnot(all(c("module", "gene", "priority_UMAP") %in% colnames(mods_long)))

# split into named list-of-vectors
module_lists <- split(mods_long$gene, mods_long$module)

# read miR-29a-3p target list and build target module (top-N strongest)
mir29 <- read.csv(mir29a_csv, stringsAsFactors = FALSE)
stopifnot(all(c("gene", "cum_weighted_context_score") %in% colnames(mir29)))
mir29 <- mir29 %>% arrange(cum_weighted_context_score)   # most negative first
top_n_targets    <- 100
module_lists$miR29a_targets_top100 <- head(mir29$gene, top_n_targets)
# also a broader version (everything with score < -0.2 = "high confidence")
module_lists$miR29a_targets_score_lt_neg0.2 <-
  mir29$gene[mir29$cum_weighted_context_score < -0.2]

# ============================ Validate gene presence =========================
present <- rownames(obj[["RNA"]])
audit <- lapply(names(module_lists), function(m) {
  g <- module_lists[[m]]; ok <- intersect(g, present)
  data.frame(module = m, n_input = length(g), n_in_data = length(ok),
             missing = paste(setdiff(g, present), collapse = ", "),
             stringsAsFactors = FALSE)
}) |> dplyr::bind_rows()
write.csv(audit, file.path(data_dir, "module_gene_presence_audit.csv"),
          row.names = FALSE)
message("Module gene presence audit written.")
# filter module lists to genes present in data
module_lists <- lapply(module_lists, function(g) intersect(g, present))
module_lists <- module_lists[lengths(module_lists) >= 3]   # need >=3 for a sensible score

# write a long-format CSV of the actual gene sets used (after presence filter)
# -- single source of truth for the methods section and supplement.
modules_used <- do.call(rbind, lapply(names(module_lists), function(m) {
  data.frame(module = m, gene = module_lists[[m]], stringsAsFactors = FALSE)
}))
write.csv(modules_used, file.path(data_dir, "module_genes_used.csv"),
          row.names = FALSE)
message("Module gene lists (used) written: module_genes_used.csv")

# ============================ AddModuleScore =================================
# Seurat appends "1" to the supplied name; we rename to drop the suffix
message("Computing module scores: ", length(module_lists), " modules")
for (m in names(module_lists)) {
  obj <- AddModuleScore(obj, features = list(module_lists[[m]]),
                        name = paste0(m, "__"), assay = "RNA",
                        ctrl = max(10, min(100, length(module_lists[[m]]) * 5)))
  # column is now "<m>__1"; rename to "score_<m>"
  raw <- paste0(m, "__1"); new <- paste0("score_", m)
  obj@meta.data[[new]] <- obj@meta.data[[raw]]
  obj@meta.data[[raw]] <- NULL
}
score_cols <- paste0("score_", names(module_lists))

# Display-friendly mirror columns (no "score_" prefix) so VlnPlot2 and
# FeaturePlot_scCustom use the module name as the natural plot title.
# Avoids the double-title (Exhaustion + score_Exhaustion) issue from stacking
# ggtitle() on top of VlnPlot2's auto-title.
for (m in names(module_lists)) {
  obj@meta.data[[m]] <- obj@meta.data[[paste0("score_", m)]]
}

# ============================ Plot helpers ===================================
priority_modules <- mods_long %>%
  dplyr::filter(priority_UMAP == TRUE) %>%
  dplyr::pull(module) %>% unique()
priority_modules <- intersect(priority_modules, names(module_lists))

# palette matches script 04 (viridis magma) for visual consistency across the project
pal <- viridis(n = 10, option = "A")

# per-module FeaturePlot of the score on UMAP -- single title from feature name
plot_score_umap <- function(feat) {
  FeaturePlot_scCustom(obj, reduction = umap_reduction, features = feat,
                       colors_use = pal, order = TRUE)
}

# per-module VlnPlot by cluster state -- matches script 04 VlnPlot conventions
plot_score_vln_cluster <- function(feat) {
  VlnPlot2(obj, features = feat, group.by = "tentative_state",
           cols = "default", show.mean = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
}

# per-module VlnPlot by condition WITHIN each cluster (the miR-29a effect view)
# Wilcoxon stat annotations match script 04's by-condition pattern.
plot_score_vln_cond_x_cluster <- function(feat) {
  VlnPlot2(obj, features = feat, group.by = "tentative_state",
           split.by = "condition", cols = "default",
           stat.method = "wilcox.test") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
}

safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# ============================ Per-module plots ===============================
message("Saving per-module plots...")
for (m in names(module_lists)) {
  ggsave(file.path(umap_dir, paste0(safe(m), ".png")),
         plot_score_umap(m),
         width = 8, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(vln_clu_dir, paste0(safe(m), ".png")),
         plot_score_vln_cluster(m),
         width = 14, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(vln_cnd_dir, paste0(safe(m), ".png")),
         plot_score_vln_cond_x_cluster(m),
         width = 14, height = 8, dpi = 300, bg = "white")
}

# ============================ Composite priority-UMAP grid ===================
if (length(priority_modules) > 0) {
  message("Composite UMAP grid for ", length(priority_modules), " priority modules")
  panels <- lapply(priority_modules, function(m) plot_score_umap(m))
  grid <- wrap_plots(panels, ncol = 3)
  ggsave(file.path(plot_dir, "UMAP_priority_modules.png"),
         grid, width = 21, height = 7 * ceiling(length(priority_modules) / 3),
         dpi = 300, bg = "white", limitsize = FALSE)
}

# ============================ miR-29a target headline tests ==================
# Wilcoxon: is the target module score lower in miR29a vs each control?
mir_tests <- list()
for (mir_col in c("score_miR29a_targets_top100",
                  "score_miR29a_targets_score_lt_neg0.2")) {
  if (!mir_col %in% colnames(obj@meta.data)) next
  md <- obj@meta.data
  # overall (across all cells)
  for (ctrl in c("EV", "Scr")) {
    x <- md[[mir_col]][md$condition == "miR29a"]
    y <- md[[mir_col]][md$condition == ctrl]
    if (length(x) > 10 && length(y) > 10) {
      w <- wilcox.test(x, y, alternative = "less")  # H1: miR29a < ctrl
      mir_tests[[length(mir_tests)+1]] <- data.frame(
        module = mir_col, level = "OVERALL", group = paste0("miR29a_vs_", ctrl),
        median_miR29a = median(x), median_ctrl = median(y),
        delta = median(x) - median(y), W = unname(w$statistic), p_value = w$p.value
      )
    }
  }
  # per cluster (tentative_state) -- where is repression strongest?
  for (cl in levels(md$tentative_state)) {
    sub <- md[md$tentative_state == cl, ]
    for (ctrl in c("EV", "Scr")) {
      x <- sub[[mir_col]][sub$condition == "miR29a"]
      y <- sub[[mir_col]][sub$condition == ctrl]
      if (length(x) > 10 && length(y) > 10) {
        w <- wilcox.test(x, y, alternative = "less")
        mir_tests[[length(mir_tests)+1]] <- data.frame(
          module = mir_col, level = as.character(cl),
          group = paste0("miR29a_vs_", ctrl),
          median_miR29a = median(x), median_ctrl = median(y),
          delta = median(x) - median(y), W = unname(w$statistic), p_value = w$p.value
        )
      }
    }
  }
}
mir_test_df <- dplyr::bind_rows(mir_tests)
if (nrow(mir_test_df) > 0) {
  mir_test_df <- mir_test_df %>%
    dplyr::group_by(module) %>%
    dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup()
  write.csv(mir_test_df,
            file.path(data_dir, "miR29a_target_module_wilcoxon.csv"),
            row.names = FALSE)
  message("Wilcoxon results written: miR29a_target_module_wilcoxon.csv")
}

# ============================ Re-save object =================================
# Save module scores to a NEW checkpoint -- do NOT overwrite PreAnnotation.qs2
# (that file is script 04's output; module scoring is a downstream stage).
qs_save(obj, file = file.path(saved_dir, "Mouse_CARTmiR29a_WithModuleScores.qs2"))
message("\nDone. Module scores saved into object metadata (columns prefixed `score_`).")
message("Object written to: Mouse_CARTmiR29a_WithModuleScores.qs2")
message("See ", out_base)
###############################################################################