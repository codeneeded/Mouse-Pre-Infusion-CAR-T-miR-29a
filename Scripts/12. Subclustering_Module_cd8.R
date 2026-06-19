# =============================================================================
# 13_metabolic_scoring_CD8_mouse_CARTmiR29a.R
# -----------------------------------------------------------------------------
# Metabolic / stress module + pathway scoring on the CD8 SUBCLUSTERED object,
# contrasting miR29a vs EV within each (collaborator-defined) coarse subcluster.
#
# Collaborator's CD8 grouping (follow the dendrogram):
#     annot_cluster 0-7  -> "Effector 1"
#     annot_cluster 8    -> "Effector 2"
#     annot_cluster 9    -> "Stem-like"
#
# Scoring: UCell (rank-based, sequencing-depth robust -- the right tool for
# comparing module activity ACROSS conditions; AddModuleScore is centering-based
# and less comparable between groups). One unified UCell scale for the custom
# modules AND the requested Hallmark/GO gene sets.
#
# Plots: scCustomize for FeaturePlots; SeuratExtend (DimPlot2/VlnPlot2/Heatmap)
# for everything else. PNG only, 300 dpi, white bg.
#
# Outputs -> Subclustering/CD8/metabolic_scoring/
#     feature_plots/   per-set FeaturePlot_scCustom, split by condition
#     violins/         per-set VlnPlot2, subcluster x condition (miR29a vs EV)
#     heatmaps/        mean-score heatmap + the miR29a-EV DELTA heatmap (headline)
#     tables/          gene-set membership, mean scores, pseudobulk, Wilcoxon
# =============================================================================

# ---- packages ---------------------------------------------------------------
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("UCell",   quietly = TRUE)) remotes::install_github("carmonalab/UCell")
if (!requireNamespace("msigdbr", quietly = TRUE)) install.packages("msigdbr")

library(Seurat)
library(SeuratExtend)
library(scCustomize)
library(UCell)
library(msigdbr)
library(qs2)
library(dplyr); library(tidyr); library(tibble)
library(ggplot2); library(patchwork)
library(pheatmap)

set.seed(42)

# ---- config -----------------------------------------------------------------
project_dir <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/10x_scRNAseq"
saved_dir   <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/saved_R_data"

cd8_obj_path <- file.path(saved_dir, "Mouse_CARTmiR29a_CD8_annotated.qs2")

out_base <- file.path(project_dir, "Subclustering", "CD8", "metabolic_scoring")
dir_feat <- file.path(out_base, "feature_plots")
dir_vln  <- file.path(out_base, "violins")
dir_heat <- file.path(out_base, "heatmaps")
dir_tab  <- file.path(out_base, "tables")
for (d in c(dir_feat, dir_vln, dir_heat, dir_tab))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Conditions / contrast. EV is the control the collaborator asked to compare to.
cond_col      <- "condition"
contrast_test <- "miR29a"        # numerator
contrast_ctrl <- "EV"            # denominator (control)
cond_levels   <- c("EV", "Scr", "miR29a")
cond_cols     <- c(EV = "grey60", Scr = "#4C72B0", miR29a = "#C44E52")

# Collaborator's coarse grouping. C7 was previously flagged as a possible
# doublet/ambient cluster but the collaborator's map folds 0-7 into Effector 1;
# we follow that. Flip exclude_c7 to TRUE to drop C7 instead.
exclude_c7 <- FALSE
coarse_levels <- c("Effector 1", "Effector 2", "Stem-like")
coarse_cols   <- c("Effector 1" = "#4C72B0", "Effector 2" = "#DD8452",
                   "Stem-like"  = "#55A868")

# Integrated embedding from script 11 (RunUMAP on harmony). The subset object
# also carries a STALE parent "umap"; always plot on the harmony one.
red_use <- "umap.harmony"

# ---- 1. LOAD + BUILD COARSE ANNOTATION ======================================
message("=== load CD8 subclustered/annotated object ===")
seu <- qs_read(cd8_obj_path)
DefaultAssay(seu) <- "RNA"

stopifnot("annot_cluster" %in% colnames(seu@meta.data))
if (!red_use %in% Reductions(seu))
  stop("Reduction '", red_use, "' not found in object (have: ",
       paste(Reductions(seu), collapse = ", "), "). Check script 11 output.")
ac <- as.character(seu$annot_cluster)
coarse <- dplyr::case_when(
  ac %in% as.character(0:7) ~ "Effector 1",
  ac == "8"                 ~ "Effector 2",
  ac == "9"                 ~ "Stem-like",
  TRUE                      ~ NA_character_)
seu$coarse_annot <- factor(coarse, levels = coarse_levels)

if (exclude_c7) {
  keep <- !(ac == "7")
  message("  excluding C7 (", sum(!keep), " cells) per exclude_c7 = TRUE")
  seu <- subset(seu, cells = colnames(seu)[keep])
}

# tidy condition factor
seu@meta.data[[cond_col]] <- factor(as.character(seu@meta.data[[cond_col]]),
                                    levels = cond_levels)
Idents(seu) <- "coarse_annot"
message("  cells per coarse group x condition:")
print(table(seu$coarse_annot, seu@meta.data[[cond_col]]))

# replicate column for pseudobulk (fall back to .sample if absent)
rep_col <- if ("replicate" %in% colnames(seu@meta.data)) "replicate" else ".sample"

# ---- 2. GENE SETS: CUSTOM MODULES + REQUESTED HALLMARK / GO ==================
message("\n=== assembling gene sets ===")

# --- collaborator's custom modules (mouse symbols, verbatim) ---
custom_sets <- list(
  Module1_MetabolicStress_Hypoxia = c(
    "Epas1","Hif1a","Atf3","Atf4","Atf5","Ddit3","Slc2a6","Slc16a3","Ldha",
    "Pdk1","Pdk3","Gpx8","Hmox1","Bnip3","Bnip3l","Irs2","Pim3","Adora2a",
    "Bhlhe40","Nr4a3"),
  Module2_MetabolicResilience_OxFitness = c(
    "Foxo3","Cpt1a","Cpt2","Ldhb","Acss1","Nqo1","Nme4","Clybl","Cox7a1","Gpd2",
    "Slc25a24","Rxra","Mgll","Ppargc1b","Tfam","Sirt3","Sod2","Prdx3","S1pr1",
    "Klf2","Tcf7"),
  Module3_MitoQualityControl = c(
    "Opa1","Mfn1","Mfn2","Dnm1l","Fis1","Mff","Yme1l1","Oma1","Pink1","Prkn",
    "Bnip3","Bnip3l","Fundc1","Sqstm1","Optn","Tbk1","Park7","Lonp1","Clpp",
    "Hspd1","Hspa9"),
  Module3A_MitoOrganization = c(
    "Opa1","Mfn1","Mfn2","Dnm1l","Fis1","Mff","Yme1l1","Oma1","Tfam","Lonp1",
    "Clpp","Hspd1","Hspa9"),
  Module3B_Mitophagy_Damage = c(
    "Pink1","Prkn","Bnip3","Bnip3l","Fundc1","Sqstm1","Optn","Tbk1","Park7"),
  Module4_Activation_Anabolic = c(
    "Myc","Mki67","Top2a","Ccnb1","Cdk1","E2f1","E2f2","E2f8","Hk2","Slc2a1",
    "Slc2a3","Pkm","Ldha","Il2ra","Ifng","Gzma","Gzmb","Ccr5","Batf","Eomes")
)

# --- requested MSigDB sets: pull ALL mouse sets once, then select by gs_name ---
# (filtering on the stable gs_name column avoids msigdbr's version-specific
#  category/collection argument differences.)
msig <- msigdbr(species = "Mus musculus")
sym_col <- if ("gene_symbol" %in% colnames(msig)) "gene_symbol" else "gene_symbol"
all_sets <- split(msig[[sym_col]], msig$gs_name)
avail <- names(all_sets)

# label -> candidate exact names (first present wins) + regex fallback
pathway_spec <- list(
  HALLMARK_HYPOXIA                       = list(exact = "HALLMARK_HYPOXIA"),
  HALLMARK_GLYCOLYSIS                    = list(exact = "HALLMARK_GLYCOLYSIS"),
  HALLMARK_OXIDATIVE_PHOSPHORYLATION     = list(exact = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
  HALLMARK_FATTY_ACID_METABOLISM         = list(exact = "HALLMARK_FATTY_ACID_METABOLISM"),
  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = list(exact = "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"),
  HALLMARK_UNFOLDED_PROTEIN_RESPONSE     = list(exact = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
  HALLMARK_MTORC1_SIGNALING              = list(exact = "HALLMARK_MTORC1_SIGNALING"),
  HALLMARK_MYC_TARGETS_V1                = list(exact = "HALLMARK_MYC_TARGETS_V1"),
  HALLMARK_MYC_TARGETS_V2                = list(exact = "HALLMARK_MYC_TARGETS_V2"),
  HALLMARK_IL2_STAT5_SIGNALING           = list(exact = "HALLMARK_IL2_STAT5_SIGNALING"),
  HALLMARK_TNFA_SIGNALING_VIA_NFKB       = list(exact = "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  HALLMARK_APOPTOSIS                     = list(exact = "HALLMARK_APOPTOSIS"),
  GO_AUTOPHAGY = list(
    exact = c("GOBP_AUTOPHAGY", "GOBP_PROCESS_UTILIZING_AUTOPHAGIC_MECHANISM"),
    grep  = "^GOBP_(MACRO)?AUTOPHAGY$|AUTOPHAGIC_MECHANISM"),
  GO_MITOPHAGY = list(exact = "GOBP_MITOPHAGY", grep = "MITOPHAGY"),
  GO_MITOCHONDRIAL_BIOGENESIS = list(
    exact = c("GOBP_MITOCHONDRION_ORGANIZATION"),
    grep  = "MITOCHONDRI.*(BIOGENESIS|ORGANIZATION)"),
  GO_RIBOSOME_BIOGENESIS = list(exact = "GOBP_RIBOSOME_BIOGENESIS",
                                grep = "RIBOSOME_BIOGENESIS$")
)

resolve_set <- function(spec) {
  hit <- spec$exact[spec$exact %in% avail]
  if (length(hit)) return(hit[1])
  if (!is.null(spec$grep)) {
    g <- grep(spec$grep, avail, value = TRUE)
    if (length(g)) return(g[which.min(nchar(g))])  # prefer the simplest name
  }
  NA_character_
}
pathway_resolved <- vapply(pathway_spec, resolve_set, character(1))
for (nm in names(pathway_resolved)) {
  if (is.na(pathway_resolved[nm]))
    message("  [!] could not resolve ", nm, " -- skipping (check MSigDB name).")
  else if (pathway_resolved[nm] != nm)
    message("  ", nm, " -> using MSigDB set: ", pathway_resolved[nm])
}
pathway_resolved <- pathway_resolved[!is.na(pathway_resolved)]
pathway_sets <- setNames(all_sets[pathway_resolved], names(pathway_resolved))

# combined set list (custom modules first, then pathways)
gene_sets <- c(custom_sets, pathway_sets)

# membership report (how many genes of each set are in the object)
present <- rownames(seu)
membership <- data.frame(
  set        = names(gene_sets),
  n_total    = vapply(gene_sets, length, integer(1)),
  n_in_object = vapply(gene_sets, function(g) sum(g %in% present), integer(1))
)
membership$pct_found <- round(100 * membership$n_in_object / membership$n_total, 1)
write.csv(membership, file.path(dir_tab, "geneset_membership.csv"), row.names = FALSE)
print(membership)

# ---- 3. SCORE WITH UCell ====================================================
message("\n=== scoring with UCell (this can take a few minutes) ===")
seu <- AddModuleScore_UCell(seu, features = gene_sets, name = "_UCell",
                            ncores = 1, maxRank = 2000)
score_cols <- paste0(names(gene_sets), "_UCell")
score_cols <- score_cols[score_cols %in% colnames(seu@meta.data)]

# pretty labels for plotting (drop _UCell, tidy HALLMARK_/GO_)
pretty_label <- function(x) {
  x <- sub("_UCell$", "", x)
  x <- gsub("HALLMARK_", "H: ", x); x <- gsub("^GO_", "GO: ", x)
  gsub("_", " ", x)
}

qs_save(seu, file.path(saved_dir, "Mouse_CARTmiR29a_CD8_metabolic_scored.qs2"))
message("  saved scored object.")

# ---- 4. UMAP CONTEXT ========================================================
message("\n=== UMAP context ===")
p_umap <- DimPlot2(seu, group.by = "coarse_annot", cols = coarse_cols,
                   reduction = red_use, box = TRUE, label = TRUE,
                   theme = theme_umap_arrows()) +
  ggtitle("CD8 coarse subclusters")
ggsave(file.path(out_base, "umap_coarse_annot.png"), p_umap,
       width = 7, height = 6, dpi = 300, bg = "white")

p_umap_cond <- DimPlot2(seu, group.by = "coarse_annot", cols = coarse_cols,
                        reduction = red_use, split.by = cond_col, box = TRUE,
                        theme = theme_umap_arrows())
ggsave(file.path(out_base, "umap_coarse_annot_by_condition.png"), p_umap_cond,
       width = 15, height = 6, dpi = 300, bg = "white")

# ---- 5. FEATURE PLOTS (scCustomize), split by condition =====================
# Restrict to the two contrasted conditions so the split shows miR29a vs EV.
message("\n=== feature plots (miR29a vs EV) ===")
seu2 <- subset(seu, subset = condition %in% c(contrast_ctrl, contrast_test))
seu2$condition <- factor(as.character(seu2$condition),
                         levels = c(contrast_ctrl, contrast_test))

for (sc in score_cols) {
  p <- FeaturePlot_scCustom(seu2, features = sc, split.by = cond_col,
                            reduction = red_use, order = TRUE, pt.size = 0.2,
                            num_columns = 2, colors_use = viridis_plasma_dark_high) &
    theme(plot.title = element_text(size = 10))
  ggsave(file.path(dir_feat, paste0("feat_", sub("_UCell$", "", sc), ".png")),
         p, width = 11, height = 5.2, dpi = 300, bg = "white")
}

# ---- 6. VIOLINS (SeuratExtend), subcluster x condition (miR29a vs EV) =======
message("\n=== violins (miR29a vs EV per subcluster) ===")
cond_cols2 <- cond_cols[c(contrast_ctrl, contrast_test)]
for (sc in score_cols) {
  p <- VlnPlot2(seu2, features = sc, group.by = "coarse_annot",
                split.by = cond_col, cols = cond_cols2,
                stat.method = "wilcox.test", hide.points = TRUE) +
    ggtitle(pretty_label(sc)) + ylab("UCell score")
  ggsave(file.path(dir_vln, paste0("vln_", sub("_UCell$", "", sc), ".png")),
         p, width = 7.5, height = 5, dpi = 300, bg = "white")
}

# ---- 7. MEAN-SCORE + DELTA HEATMAPS =========================================
message("\n=== heatmaps ===")
meta <- seu@meta.data

# mean UCell per coarse group x condition (all three conditions)
mean_tbl <- meta %>%
  dplyr::select(dplyr::all_of(c("coarse_annot", cond_col, score_cols))) %>%
  tidyr::pivot_longer(dplyr::all_of(score_cols), names_to = "set", values_to = "score") %>%
  dplyr::group_by(coarse_annot, !!rlang::sym(cond_col), set) %>%
  dplyr::summarise(mean_score = mean(score), .groups = "drop")
write.csv(mean_tbl, file.path(dir_tab, "mean_scores_per_subcluster_condition.csv"),
          row.names = FALSE)

# matrix for the mean heatmap: rows = sets, cols = group|condition
mat_long <- mean_tbl %>%
  tidyr::unite("col", coarse_annot, !!rlang::sym(cond_col), sep = " | ") %>%
  tidyr::pivot_wider(names_from = col, values_from = mean_score) %>%
  tibble::column_to_rownames("set")
rownames(mat_long) <- pretty_label(rownames(mat_long))
pheatmap(as.matrix(mat_long), scale = "row", cluster_cols = FALSE, cluster_rows = TRUE,
         main = "CD8: mean UCell score (row-scaled) per subcluster x condition",
         color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(101),
         filename = file.path(dir_heat, "mean_score_heatmap.png"),
         width = 11, height = 9)

# ---- 8. STATS: per-cell Wilcoxon miR29a vs EV per subcluster per set ========
message("\n=== Wilcoxon (miR29a vs EV) per subcluster ===")
stat_rows <- list()
for (grp in coarse_levels) {
  for (sc in score_cols) {
    x <- meta[[sc]][meta$coarse_annot == grp & meta[[cond_col]] == contrast_test]
    y <- meta[[sc]][meta$coarse_annot == grp & meta[[cond_col]] == contrast_ctrl]
    if (length(x) < 10 || length(y) < 10) next
    wt <- suppressWarnings(wilcox.test(x, y))
    # rank-biserial effect size (robust, bounded -1..1)
    auc <- (wt$statistic) / (length(x) * length(y))
    rb  <- as.numeric(2 * auc - 1)
    stat_rows[[length(stat_rows) + 1]] <- data.frame(
      coarse_annot = grp, set = pretty_label(sc),
      n_miR29a = length(x), n_EV = length(y),
      median_miR29a = median(x), median_EV = median(y),
      delta_mean = mean(x) - mean(y),
      rank_biserial = round(rb, 3),
      p_value = wt$p.value)
  }
}
stat_df <- dplyr::bind_rows(stat_rows)
stat_df$p_adj_BH <- p.adjust(stat_df$p_value, method = "BH")
stat_df <- stat_df %>% dplyr::arrange(coarse_annot, dplyr::desc(abs(delta_mean)))
write.csv(stat_df, file.path(dir_tab, "wilcoxon_miR29a_vs_EV_per_subcluster.csv"),
          row.names = FALSE)
print(utils::head(stat_df, 20))

# pseudobulk (per replicate) means -- for biological-reproducibility assessment
pb <- meta %>%
  dplyr::select(dplyr::all_of(c("coarse_annot", cond_col, rep_col, score_cols))) %>%
  tidyr::pivot_longer(dplyr::all_of(score_cols), names_to = "set", values_to = "score") %>%
  dplyr::group_by(coarse_annot, !!rlang::sym(cond_col), !!rlang::sym(rep_col), set) %>%
  dplyr::summarise(mean_score = mean(score), n_cells = dplyr::n(), .groups = "drop")
pb$set <- pretty_label(pb$set)
write.csv(pb, file.path(dir_tab, "pseudobulk_per_sample.csv"), row.names = FALSE)

# ---- 9. DELTA HEATMAP (headline: miR29a - EV) ===============================
# rows = sets (in input order), cols = subclusters; fill = mean delta; * = BH<0.05
set_order  <- pretty_label(score_cols)
delta_mat  <- matrix(NA_real_, nrow = length(set_order), ncol = length(coarse_levels),
                     dimnames = list(set_order, coarse_levels))
star_mat   <- matrix("",      nrow = length(set_order), ncol = length(coarse_levels),
                     dimnames = list(set_order, coarse_levels))
for (i in seq_len(nrow(stat_df))) {
  r <- stat_df$set[i]; c <- as.character(stat_df$coarse_annot[i])
  delta_mat[r, c] <- stat_df$delta_mean[i]
  star_mat[r, c]  <- ifelse(stat_df$p_adj_BH[i] < 0.001, "***",
                            ifelse(stat_df$p_adj_BH[i] < 0.01,  "**",
                                   ifelse(stat_df$p_adj_BH[i] < 0.05,  "*", "")))
}
lim <- max(abs(delta_mat), na.rm = TRUE)
pheatmap(delta_mat, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = star_mat, number_color = "black", fontsize_number = 12,
         breaks = seq(-lim, lim, length.out = 101),
         color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100),
         main = "CD8: miR29a - EV mean UCell score (BH: * .05  ** .01  *** .001)",
         filename = file.path(dir_heat, "delta_miR29a_minus_EV_heatmap.png"),
         width = 7.5, height = 10)

message("\nDone. Outputs in: ", out_base)