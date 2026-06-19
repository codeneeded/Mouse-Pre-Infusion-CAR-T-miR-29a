# =============================================================================
# 14_metabolic_followup_CD8_mouse_CARTmiR29a.R
# -----------------------------------------------------------------------------
# Follow-up analyses on the UCell-scored CD8 subclustered object (script 13),
# built to develop the miR29a-vs-EV metabolic story:
#
#   1. DIFFERENTIAL ABUNDANCE  - does miR29a shift the stem<->effector balance?
#                                (per-replicate proportions; descriptive, n=2/cond)
#   2. DIFFERENTIAL-DENSITY UMAP- where on the harmony embedding does miR29a
#                                accumulate / deplete cells? (miR29a - EV density)
#   3. BALANCE AXES            - the collaborator's yin/yang pairs as derived
#                                per-cell scores: resilience-stress (M2-M1) and
#                                mito maintenance-mitophagy (M3A-M3B)
#   4. MECHANISM LINK          - canonical miR-29a target repression in the
#                                subclusters + target score vs resilience
#   5. LEADING-EDGE GENES      - which member genes drive the top-moving pathways
#   6. SUMMARY DUMBBELLS       - EV->miR29a shift across all modules per subcluster
#
# Reads the scored object from script 13. UCell scale throughout. EV is control.
# scCustomize for FeaturePlots; SeuratExtend otherwise. Harmony UMAP only.
# Outputs -> Subclustering/CD8/metabolic_scoring/followup/<section>/
# =============================================================================

if (!requireNamespace("ggridges", quietly = TRUE)) install.packages("ggridges")

library(Seurat)
library(SeuratExtend)
library(scCustomize)
library(UCell)
library(qs2)
library(dplyr); library(tidyr); library(tibble)
library(ggplot2); library(patchwork); library(ggridges)
library(pheatmap)
library(MASS)            # kde2d for differential density

set.seed(42)

# ---- config -----------------------------------------------------------------
project_dir <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/10x_scRNAseq"
saved_dir   <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/saved_R_data"
scored_path <- file.path(saved_dir, "Mouse_CARTmiR29a_CD8_metabolic_scored.qs2")

out_base <- file.path(project_dir, "Subclustering", "CD8", "metabolic_scoring", "followup")
dirs <- list(abund = "abundance", dens = "density", bal = "balance",
             mech = "mechanism", lead = "leading_edge", summ = "summary",
             tab = "tables")
dirs <- lapply(dirs, function(x) file.path(out_base, x))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

cond_col      <- "condition"
contrast_test <- "miR29a"
contrast_ctrl <- "EV"
cond_levels   <- c("EV", "Scr", "miR29a")
cond_cols     <- c(EV = "grey60", Scr = "#4C72B0", miR29a = "#C44E52")
cond_cols2    <- cond_cols[c(contrast_ctrl, contrast_test)]
coarse_levels <- c("Effector 1", "Effector 2", "Stem-like")
coarse_cols   <- c("Effector 1" = "#4C72B0", "Effector 2" = "#DD8452", "Stem-like" = "#55A868")
red_use       <- "umap.harmony"
rep_col       <- "replicate"

# MASS::select() masks dplyr::select -- be explicit downstream.
sel <- dplyr::select

# ---- load -------------------------------------------------------------------
message("=== load scored CD8 object ===")
seu <- qs_read(scored_path)
DefaultAssay(seu) <- "RNA"
seu@meta.data[[cond_col]] <- factor(as.character(seu@meta.data[[cond_col]]), levels = cond_levels)
if (!red_use %in% Reductions(seu))
  stop("Reduction '", red_use, "' not found (have: ", paste(Reductions(seu), collapse = ", "), ").")
if (!rep_col %in% colnames(seu@meta.data)) rep_col <- ".sample"

# resolve UCell score columns and a basename->column lookup
score_cols <- grep("_UCell$", colnames(seu@meta.data), value = TRUE)
sbase <- function(x) sub("_UCell$", "", x)
getcol <- function(base) {
  hit <- score_cols[sbase(score_cols) == base]
  if (!length(hit)) stop("score column not found for: ", base); hit[1]
}
pretty_label <- function(x) {
  x <- sbase(x); x <- gsub("HALLMARK_", "H: ", x); x <- gsub("^GO_", "GO: ", x); gsub("_", " ", x)
}

m1 <- getcol("Module1_MetabolicStress_Hypoxia")
m2 <- getcol("Module2_MetabolicResilience_OxFitness")
m3a <- getcol("Module3A_MitoOrganization")
m3b <- getcol("Module3B_Mitophagy_Damage")

# two-condition subset (EV vs miR29a) reused throughout
seu2 <- subset(seu, subset = condition %in% c(contrast_ctrl, contrast_test))
seu2$condition <- factor(as.character(seu2$condition), levels = c(contrast_ctrl, contrast_test))

# ============================================================================
# 1. DIFFERENTIAL ABUNDANCE (per-replicate proportions; descriptive)
# ============================================================================
message("\n=== 1. differential abundance ===")
prop_tbl <- seu@meta.data %>%
  sel(coarse_annot, !!cond_col, !!rep_col) %>%
  dplyr::count(coarse_annot, !!rlang::sym(cond_col), !!rlang::sym(rep_col), name = "n") %>%
  dplyr::group_by(!!rlang::sym(cond_col), !!rlang::sym(rep_col)) %>%
  dplyr::mutate(prop = n / sum(n)) %>% dplyr::ungroup()
write.csv(prop_tbl, file.path(dirs$tab, "coarse_proportions_per_replicate.csv"), row.names = FALSE)

p_abund <- ggplot(prop_tbl, aes(x = coarse_annot, y = prop, fill = .data[[cond_col]])) +
  stat_summary(fun = mean, geom = "col", position = position_dodge(0.8),
               width = 0.7, alpha = 0.85) +
  geom_point(aes(group = .data[[cond_col]]), position = position_dodge(0.8),
             size = 1.8, colour = "black", show.legend = FALSE) +
  scale_fill_manual(values = cond_cols) +
  labs(title = "CD8 subcluster composition (bars = mean, dots = replicates)",
       x = NULL, y = "fraction of cells", fill = "condition") +
  theme_classic(base_size = 12)
ggsave(file.path(dirs$abund, "coarse_proportions.png"), p_abund,
       width = 8, height = 5, dpi = 300, bg = "white")

# miR29a - EV mean proportion shift per coarse group (descriptive; n=2/cond)
prop_shift <- prop_tbl %>%
  dplyr::group_by(coarse_annot, !!rlang::sym(cond_col)) %>%
  dplyr::summarise(mean_prop = mean(prop), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = !!cond_col, values_from = mean_prop) %>%
  dplyr::mutate(delta_miR_minus_EV = .data[[contrast_test]] - .data[[contrast_ctrl]])
write.csv(prop_shift, file.path(dirs$tab, "coarse_proportion_shift.csv"), row.names = FALSE)
print(prop_shift)

# fine-cluster proportions (supplementary)
fine_tbl <- seu@meta.data %>%
  dplyr::count(annot_cluster, !!rlang::sym(cond_col), !!rlang::sym(rep_col), name = "n") %>%
  dplyr::group_by(!!rlang::sym(cond_col), !!rlang::sym(rep_col)) %>%
  dplyr::mutate(prop = n / sum(n)) %>% dplyr::ungroup()
write.csv(fine_tbl, file.path(dirs$tab, "fine_proportions_per_replicate.csv"), row.names = FALSE)

# ============================================================================
# 2. DIFFERENTIAL-DENSITY UMAP (miR29a - EV cell density on harmony embedding)
# ============================================================================
message("\n=== 2. differential-density UMAP ===")
emb <- as.data.frame(Embeddings(seu2, red_use))[, 1:2]
colnames(emb) <- c("UMAP_1", "UMAP_2")
emb$condition <- seu2$condition
# Estimate each condition's 2D density on a common grid, then read the
# (miR29a - EV) difference back onto EVERY cell's own coordinates. Plotting the
# cells themselves (not a filled grid) means this panel is the SAME scatter as
# the coarse_annot DimPlot -- it overlays 1:1 instead of painting empty corners.
all_emb <- as.data.frame(Embeddings(seu, red_use))[, 1:2]
colnames(all_emb) <- c("UMAP_1", "UMAP_2")
rng_x <- range(all_emb$UMAP_1); rng_y <- range(all_emb$UMAP_2); ngrid <- 200
kde_z <- function(cc) {
  MASS::kde2d(emb$UMAP_1[emb$condition == cc], emb$UMAP_2[emb$condition == cc],
              n = ngrid, lims = c(rng_x, rng_y))$z
}
gx <- MASS::kde2d(emb$UMAP_1, emb$UMAP_2, n = ngrid, lims = c(rng_x, rng_y))
diffmat <- kde_z(contrast_test) / sum(kde_z(contrast_test)) -
  kde_z(contrast_ctrl) / sum(kde_z(contrast_ctrl))   # normalised surfaces

# look up the difference at each cell (nearest grid node), for ALL cells
ix <- pmin(pmax(findInterval(all_emb$UMAP_1, gx$x), 1), ngrid)
iy <- pmin(pmax(findInterval(all_emb$UMAP_2, gx$y), 1), ngrid)
all_emb$dens_diff <- diffmat[cbind(ix, iy)]
all_emb <- all_emb[order(abs(all_emb$dens_diff)), ]   # extreme cells drawn on top
lim <- max(abs(all_emb$dens_diff))

p_dens <- ggplot(all_emb, aes(UMAP_1, UMAP_2, colour = dens_diff)) +
  geom_point(size = 0.4, stroke = 0) +
  scale_colour_gradient2(low = "#3B4CC0", mid = "grey90", high = "#B40426",
                         limits = c(-lim, lim), name = "miR29a - EV\ndensity") +
  ggtitle("miR29a - EV cell density") +
  theme_umap_arrows() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(colour = guide_colourbar(barheight = 8))
ggsave(file.path(dirs$dens, "differential_density_umap.png"), p_dens,
       width = 7, height = 6, dpi = 300, bg = "white")

# ============================================================================
# 3. BALANCE AXES (derived per-cell scores)
# ============================================================================
message("\n=== 3. balance axes ===")
seu$resilience_minus_stress <- seu@meta.data[[m2]]  - seu@meta.data[[m1]]
seu$mito_maint_minus_phagy  <- seu@meta.data[[m3a]] - seu@meta.data[[m3b]]
seu2$resilience_minus_stress <- seu2@meta.data[[m2]]  - seu2@meta.data[[m1]]
seu2$mito_maint_minus_phagy  <- seu2@meta.data[[m3a]] - seu2@meta.data[[m3b]]

balance_axes <- c(resilience_minus_stress = "Resilience - Stress  (M2 - M1)",
                  mito_maint_minus_phagy  = "Mito maintenance - mitophagy  (M3A - M3B)")
for (ax in names(balance_axes)) {
  # split violins by subcluster x condition with Wilcoxon
  pv <- VlnPlot2(seu2, features = ax, group.by = "coarse_annot", split.by = cond_col,
                 cols = cond_cols2, stat.method = "wilcox.test", hide.points = TRUE) +
    ggtitle(balance_axes[ax]) + ylab("UCell score difference") +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40")
  ggsave(file.path(dirs$bal, paste0("vln_", ax, ".png")), pv,
         width = 7.5, height = 5, dpi = 300, bg = "white")
  # ridge plot: distribution shift by condition, faceted by subcluster
  rd <- seu2@meta.data %>% sel(coarse_annot, !!cond_col, !!ax)
  pr <- ggplot(rd, aes(x = .data[[ax]], y = coarse_annot, fill = .data[[cond_col]])) +
    geom_density_ridges(alpha = 0.6, scale = 1.1, colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    scale_fill_manual(values = cond_cols2) +
    labs(title = balance_axes[ax], x = "score difference", y = NULL, fill = "condition") +
    theme_ridges(font_size = 11)
  ggsave(file.path(dirs$bal, paste0("ridge_", ax, ".png")), pr,
         width = 8, height = 5, dpi = 300, bg = "white")
}

# 2D stress(M1) vs resilience(M2) landscape, density contours by condition
land <- seu2@meta.data %>% sel(coarse_annot, !!cond_col, !!m1, !!m2)
colnames(land)[3:4] <- c("Stress_M1", "Resilience_M2")
p_land <- ggplot(land, aes(Stress_M1, Resilience_M2, colour = .data[[cond_col]])) +
  geom_density_2d(linewidth = 0.4) +
  scale_colour_manual(values = cond_cols2) +
  facet_wrap(~ coarse_annot, scales = "free") +
  labs(title = "Stress vs resilience landscape (density contours, miR29a vs EV)",
       x = "Module 1: metabolic stress / hypoxia", y = "Module 2: resilience / ox-fitness",
       colour = "condition") +
  theme_bw(base_size = 11)
ggsave(file.path(dirs$bal, "landscape_M1_vs_M2.png"), p_land,
       width = 12, height = 4.5, dpi = 300, bg = "white")

# ============================================================================
# 4. MECHANISM LINK: canonical miR-29a target repression -> resilience
# ============================================================================
message("\n=== 4. mechanism link ===")
targets <- c("Eomes", "Dnmt3a", "Dnmt3b", "Tet2", "Tet3", "Mycn")
targets <- intersect(targets, rownames(seu2))

# confirm target repression in the subclusters (miR29a vs EV)
pt <- VlnPlot2(seu2, features = targets, group.by = "coarse_annot", split.by = cond_col,
               cols = cond_cols2, stat.method = "wilcox.test", hide.points = TRUE, ncol = 3)
ggsave(file.path(dirs$mech, "canonical_target_violins.png"), pt,
       width = 13, height = 8, dpi = 300, bg = "white")

# target module score vs resilience axis, coloured by condition
seu2 <- AddModuleScore_UCell(seu2, features = list(miR29a_targets = targets),
                             name = "_UCell", ncores = 1, maxRank = 2000)
md <- seu2@meta.data
md$target_score <- md[["miR29a_targets_UCell"]]
md$resilience   <- md[[m2]]
corr_tbl <- md %>% dplyr::group_by(coarse_annot, !!rlang::sym(cond_col)) %>%
  dplyr::summarise(r_target_vs_resilience = cor(target_score, resilience, method = "spearman"),
                   mean_target = mean(target_score), mean_resilience = mean(resilience),
                   .groups = "drop")
write.csv(corr_tbl, file.path(dirs$tab, "target_vs_resilience_correlation.csv"), row.names = FALSE)
print(corr_tbl)

p_scatter <- ggplot(md, aes(target_score, resilience, colour = .data[[cond_col]])) +
  geom_point(size = 0.2, alpha = 0.25) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = cond_cols2) +
  facet_wrap(~ coarse_annot, scales = "free") +
  labs(title = "miR-29a target score vs metabolic resilience (per cell)",
       x = "miR-29a target module (UCell)", y = "Module 2 resilience (UCell)",
       colour = "condition") +
  theme_bw(base_size = 11)
ggsave(file.path(dirs$mech, "target_vs_resilience_scatter.png"), p_scatter,
       width = 12, height = 4.5, dpi = 300, bg = "white")

# ============================================================================
# 5. LEADING-EDGE GENES for the top-moving pathways (in Effector 1, powered)
# ============================================================================
message("\n=== 5. leading-edge genes ===")
# member-gene sets for the pathways that moved most (edit as needed)
lead_sets <- list(
  TNFA_NFKB   = c("Nfkbia","Nfkbie","Tnfaip3","Rel","Relb","Cd83","Bcl2a1d","Bcl2a1b",
                  "Nr4a1","Nr4a2","Nr4a3","Egr1","Junb","Ier3","Gadd45b","Tnf","Ccl4","Ccl3"),
  Apoptosis   = c("Bax","Bak1","Bid","Bcl2l11","Pmaip1","Casp3","Casp8","Cflar","Gsn",
                  "Bcl2","Mcl1","Birc3","Cd44","Sod1","Gadd45a"),
  Module1_Stress = c("Epas1","Hif1a","Atf3","Atf4","Atf5","Ddit3","Slc2a6","Slc16a3",
                     "Ldha","Pdk1","Pdk3","Gpx8","Hmox1","Bnip3","Bnip3l","Irs2","Pim3",
                     "Adora2a","Bhlhe40","Nr4a3"))

eff1 <- subset(seu2, subset = coarse_annot == "Effector 1")
Idents(eff1) <- cond_col
avg <- AverageExpression(eff1, assays = "RNA", layer = "data", group.by = cond_col)$RNA
colnames(avg) <- sub("^g", "", colnames(avg))
for (nm in names(lead_sets)) {
  g <- intersect(lead_sets[[nm]], rownames(avg))
  if (length(g) < 3) next
  mat <- log1p(as.matrix(avg[g, c(contrast_ctrl, contrast_test), drop = FALSE]))
  lfc <- mat[, contrast_test] - mat[, contrast_ctrl]
  mat <- mat[order(lfc), , drop = FALSE]            # order by EV->miR direction
  pheatmap(mat, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
           main = paste0("Effector 1 leading edge: ", nm, " (log mean expr)"),
           color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(101),
           filename = file.path(dirs$lead, paste0("leadingedge_", nm, ".png")),
           width = 4.5, height = max(4, 0.28 * nrow(mat)))
}

# ============================================================================
# 6. SUMMARY DUMBBELLS: EV -> miR29a shift across all modules, per subcluster
# ============================================================================
message("\n=== 6. summary dumbbells ===")
mean_long <- seu2@meta.data %>%
  sel(coarse_annot, !!cond_col, dplyr::all_of(score_cols)) %>%
  tidyr::pivot_longer(dplyr::all_of(score_cols), names_to = "set", values_to = "score") %>%
  dplyr::group_by(coarse_annot, !!rlang::sym(cond_col), set) %>%
  dplyr::summarise(mean_score = mean(score), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = !!cond_col, values_from = mean_score) %>%
  dplyr::mutate(set = pretty_label(set),
                delta = .data[[contrast_test]] - .data[[contrast_ctrl]])
write.csv(mean_long, file.path(dirs$tab, "mean_scores_wide_with_delta.csv"), row.names = FALSE)

for (grp in coarse_levels) {
  d <- mean_long %>% dplyr::filter(coarse_annot == grp) %>% dplyr::arrange(delta)
  d$set <- factor(d$set, levels = d$set)
  pd <- ggplot(d) +
    geom_segment(aes(x = .data[[contrast_ctrl]], xend = .data[[contrast_test]],
                     y = set, yend = set), colour = "grey70", linewidth = 1) +
    geom_point(aes(x = .data[[contrast_ctrl]], y = set), colour = cond_cols2[[1]], size = 2.6) +
    geom_point(aes(x = .data[[contrast_test]], y = set), colour = cond_cols2[[2]], size = 2.6) +
    labs(title = paste0(grp, ": EV (grey) -> miR29a (red) per module"),
         x = "mean UCell score", y = NULL) +
    theme_classic(base_size = 11)
  ggsave(file.path(dirs$summ, paste0("dumbbell_", gsub("[^A-Za-z0-9]", "_", grp), ".png")),
         pd, width = 8, height = 7, dpi = 300, bg = "white")
}

qs_save(seu, file.path(saved_dir, "Mouse_CARTmiR29a_CD8_metabolic_scored.qs2"))  # adds balance axes
message("\nDone. Follow-up outputs in: ", out_base)