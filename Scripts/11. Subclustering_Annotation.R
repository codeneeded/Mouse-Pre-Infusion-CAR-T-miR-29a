###############################################################################
# 12_annotation_support_CD8_CD4_mouse_CARTmiR29a.R
#
# ANNOTATION DECISION-SUPPORT for the CD8 (res 0.5) and CD4 (res 0.6)
# subclustered objects from script 11. This script does NOT assign final
# labels -- it produces evidence + suggestions for the wetlab biologist, who
# makes the final call. Specifically, per lineage it generates:
#
#   1. AUTOMATIC annotation        -- ProjecTILs projection onto the mouse LCMV
#                                     CD8/CD4 reference atlases: embeds the query
#                                     in reference space, predicts per-cell
#                                     functional.cluster, and emits the tutorial
#                                     diagnostics (projection UMAP, state
#                                     composition, marker radar). Reference-based
#                                     and EXPECTED to fit imperfectly for in-vitro
#                                     pre-infusion cells -- a cross-check, not
#                                     ground truth.
#   2/3. MANUAL-annotation aids     -- per-cluster QC table, curated marker
#                                     dot plot + avg-expression heatmap,
#                                     signature module-score heatmap, and the
#                                     per-cluster marker tables from script 11.
#   4. SUGGESTED cluster names      -- a worksheet CSV combining top signature,
#                                     ProjecTILs majority label, top markers,
#                                     and a tentative name (clearly a SUGGESTION).
#   5. SUGGESTED merges             -- cluster correlation heatmap + dendrogram,
#                                     plus targeted pairwise DE counts on highly
#                                     correlated pairs -> merge_candidates.csv.
#
# STRUCTURE: flat, same convention as script 11. Shared config + curated
# signatures at the top, then CD8 spelled out in numbered sections, then CD4.
# Run/inspect one section at a time; seu_cd8 / seu_cd4 stay in the environment.
#
# Inputs (from script 11):
#   saved_R_data/Mouse_CARTmiR29a_CD8_subclustered.qs2
#   saved_R_data/Mouse_CARTmiR29a_CD4_subclustered.qs2
#   Subclustering/<lineage>/Markers/top20_per_cluster_res<r>.csv
#
# Outputs:
#   Subclustering/<lineage>/annotation/  (QC, dotplot, heatmaps, dendrogram,
#     merge + name CSVs, and ProjecTILs projection/composition/radar/crosstab)
#   saved_R_data/Mouse_CARTmiR29a_<lineage>_annotated.qs2  (carries scores +
#     ProjecTILs labels + the chosen-resolution ident as annot_cluster; NO final
#     names -- those are added in a later step once the biologist decides)
###############################################################################

# ProjecTILs (+ its carmonalab deps) -- install once, in dependency order.
# STACAS is NOT on CRAN, so ProjecTILs' own install cannot auto-resolve it;
# the deps must be pulled from GitHub first, bottom-up: UCell -> scGate ->
# STACAS -> ProjecTILs. R_REMOTES_NO_ERRORS_FROM_WARNINGS stops a dependency
# build *warning* from aborting the install with a non-zero exit status.
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
if (!requireNamespace("remotes",    quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("UCell",      quietly = TRUE)) remotes::install_github("carmonalab/UCell")
if (!requireNamespace("scGate",     quietly = TRUE)) remotes::install_github("carmonalab/scGate")
if (!requireNamespace("STACAS",     quietly = TRUE)) remotes::install_github("carmonalab/STACAS")
if (!requireNamespace("ProjecTILs", quietly = TRUE)) remotes::install_github("carmonalab/ProjecTILs")

library(Seurat)
library(SeuratExtend)
library(scCustomize)
library(qs2)
library(dplyr); library(tidyr); library(tibble)
library(ggplot2); library(patchwork)
library(pheatmap)
library(ProjecTILs)

# ============================ Paths ==========================================
project_dir <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/10x_scRNAseq"
saved_dir   <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/saved_R_data"
sub_base    <- file.path(project_dir, "Subclustering")
# Annotation outputs live under each lineage's subcluster folder (mirrors
# script 11's Subclustering/<lineage>/ layout): Subclustering/<lineage>/annotation/

# ============================ Parameters (REVIEW) ============================
res_cd8 <- 0.5                              # committed CD8 resolution (script 11)
res_cd4 <- 0.6                              # committed CD4 resolution (script 11)

# Oddball clusters flagged during sweep review (possible doublet/ambient). The
# QC table will quantify; these are pre-flagged for the biologist's attention.
flag_clusters_cd8 <- c("7")                 # neuronal-ish genes on G2/M backbone
flag_clusters_cd4 <- c("12")                # Mctp2/Myo10/Il31ra oddball

# Merge-suggestion thresholds
merge_corr_threshold <- 0.90                # avg-expr Pearson r above which a
# pair is a merge CANDIDATE
de_padj <- 0.05; de_lfc <- 0.5              # "DE gene" definition for pair tests

# ---- FINALISED CD8 proposal (from manual markers + signatures + QC + ProjecTILs)
# cluster -> fine label and a coarse merged group. Drives proposed_annotation_CD8
# .csv. Edit freely; CD4 left empty until its manual review is done.
cd8_proposed_fine <- c(
  "0" = "Effector - cytotoxic",
  "1" = "Proliferating effector (S-phase)",
  "2" = "Proliferating effector (S/G2, histone)",
  "3" = "Proliferating effector (G2/M)",
  "4" = "Proliferating effector (mitotic, sterol-high)",
  "5" = "ISR / metabolic-stress effector (cycling)",
  "6" = "ISR / metabolic-stress effector (non-cycling)",
  "7" = "EXCLUDE - doublet / ambient",
  "8" = "Activated effector (cytokine+)",
  "9" = "Stem-like / memory-precursor (Tcf7+)")
cd8_proposed_group <- c(
  "0" = "Effector", "8" = "Activated effector",
  "1" = "Proliferating", "2" = "Proliferating", "3" = "Proliferating", "4" = "Proliferating",
  "5" = "ISR/stress", "6" = "ISR/stress",
  "9" = "Stem/precursor", "7" = "EXCLUDE")
cd4_proposed_fine  <- character(0)          # CD4 manual review pending
cd4_proposed_group <- character(0)

# ProjecTILs reference maps (mouse LCMV atlases; Andreatta & Carmona, eLife 2022).
# We PROJECT our query onto THEIR reference -- fully supported. (The Harmony
# limitation only bites if you try to BUILD a reference from a Harmony object.)
# Refs are LCMV-derived, so in-vitro pre-infusion cells fit imperfectly: treat
# the transferred labels as a CROSS-CHECK, not ground truth. Files auto-download
# into saved_dir on first run (URLs are the official figshare ndownloader links
# used in the ProjecTILs case studies).
projectils_ref_cd8 <- file.path(saved_dir, "ref_LCMV_Atlas_mouse_v1.rds")
projectils_ref_cd4 <- file.path(saved_dir, "ref_LCMV_CD4_mouse_release_v1.rds")
projectils_url_cd8 <- "https://ndownloader.figshare.com/files/23166794"   # CD8 LCMV atlas
projectils_url_cd4 <- "https://ndownloader.figshare.com/files/31057081"   # CD4 LCMV atlas

# Genes for the per-state radar plots (query vs reference), one panel per lineage.
genes4radar_cd8 <- c("Cd8a","Tcf7","Ccr7","Sell","Gzmb","Gzmk","Slamf6",
                     "Pdcd1","Havcr2","Tox","Cx3cr1","Mki67")
genes4radar_cd4 <- c("Cd4","Foxp3","Il2ra","Tcf7","Ccr7","Sell","Tbx21",
                     "Ifng","Cxcr5","Bcl6","Gzmb","Mki67")

# Per-sample projection: NULL projects all cells together (safe default). Only
# set to ".sample" if individual samples have >~1000 cells (ProjecTILs' batch
# handling needs the cells); CD8 is the minority compartment so NULL is wise.
projectils_split_cd8 <- NULL
projectils_split_cd4 <- NULL

# ============================ Curated signatures =============================
# Safe (no spaces/slashes) keys for column names; display names below.
cd8_signatures <- list(
  Naive_Tcm        = c("Tcf7","Lef1","Sell","Ccr7","Il7r","Bcl2","Pik3ip1"),
  Effector_cytotox = c("Gzmb","Gzmc","Prf1","Ifng","Nkg7","Klrg1","Cx3cr1"),
  Exhaustion       = c("Havcr2","Lag3","Pdcd1","Tox","Entpd1","Tigit","Ctla4"),
  Proliferation    = c("Mki67","Top2a","Birc5","Ube2c","Pcna","Cenpf","Plk1"),
  ISR_stress       = c("Ddit3","Atf4","Atf5","Trib3","Chac1","Nupr1","Slc7a11","Asns","Gpt2"),
  Activation       = c("Ccl3","Ccl4","Xcl1","Csf2","Egr2","Nr4a3","Ifng","Pdcd1")
)
cd8_name_lookup <- c(
  Naive_Tcm        = "Tcf7+ progenitor / stem-like",
  Effector_cytotox = "Effector (cytotoxic)",
  Exhaustion       = "Exhaustion-like",
  Proliferation    = "Proliferating (cell cycle)",
  ISR_stress       = "ISR / amino-acid stress",
  Activation       = "Activated / cytokine effector"
)

cd4_signatures <- list(
  Treg             = c("Foxp3","Ikzf2","Il2ra","Ctla4","Ikzf4","Tnfrsf18"),
  Th1              = c("Tbx21","Ifng","Cxcr3"),
  Naive_mem        = c("Tcf7","Sell","Ccr7","Il7r","S1pr1","Pik3ip1"),
  Proliferation    = c("Mki67","Top2a","Birc5","Ube2c","Pcna","Cenpf","Plk1"),
  ISR_stress       = c("Nupr1","Trib3","Ddit3","Atf5","Asns","Eif4ebp1"),
  AP1_early        = c("Jun","Fos","Egr1","Egr2","Egr3","Nr4a1","Nr4a3"),
  IFN_ISG          = c("Cxcl10","Ccl5","Ifit1","Ifit3","Isg15","Oasl2"),
  Th2_helper       = c("Il4","Il24","Ccr8","Maf","Rora")
)
cd4_name_lookup <- c(
  Treg             = "Treg (Foxp3+)",
  Th1              = "Th1 effector (Tbx21/Ifng)",
  Naive_mem        = "Quiescent / memory-like",
  Proliferation    = "Proliferating (cell cycle)",
  ISR_stress       = "ISR / amino-acid stress",
  AP1_early        = "Immediate-early activated (AP-1)",
  IFN_ISG          = "IFN / ISG responder",
  Th2_helper       = "Th2-like / Il4+ helper"
)


###############################################################################
###############################################################################
##                              C D 8                                        ##
###############################################################################
###############################################################################

cd8_dir <- file.path(sub_base, "CD8", "annotation"); dir.create(cd8_dir, recursive = TRUE, showWarnings = FALSE)
# Organised outputs: manual/ and automatic/ portions, each with a tables/
# subdir for CSVs (PNGs sit directly in the portion folder). Synthesis
# files (suggested_names, proposed_annotation, README) live at the top.
cd8_man     <- file.path(cd8_dir, "manual");        dir.create(file.path(cd8_man, "tables"), recursive = TRUE, showWarnings = FALSE)
cd8_man_tab <- file.path(cd8_man, "tables")
cd8_aut     <- file.path(cd8_dir, "automatic");     dir.create(file.path(cd8_aut, "tables"), recursive = TRUE, showWarnings = FALSE)
cd8_aut_tab <- file.path(cd8_aut, "tables")

# ######################## PORTION A : MANUAL ANNOTATION ######################
# ===== CD8 : 1. LOAD + SET IDENT FROM COMMITTED RESOLUTION + QC TABLE ========
message("\n=== CD8 : load + QC ===")
seu_cd8 <- qs_read(file.path(saved_dir, "Mouse_CARTmiR29a_CD8_subclustered.qs2"))
DefaultAssay(seu_cd8) <- "RNA"

# Derive the working clustering DIRECTLY from the committed-resolution column
# (robust to whatever $subclusters happened to be saved as).
res_col_cd8 <- paste0("RNA_snn_res.", res_cd8)
stopifnot(res_col_cd8 %in% colnames(seu_cd8@meta.data))
seu_cd8$annot_cluster <- factor(seu_cd8[[res_col_cd8]][, 1],
                                levels = sort(as.integer(levels(seu_cd8[[res_col_cd8]][, 1]))))
Idents(seu_cd8) <- "annot_cluster"
clusters_cd8 <- levels(seu_cd8$annot_cluster)

# percent.mt (compute fresh; mouse mito genes are "mt-")
if (!"percent_mt" %in% colnames(seu_cd8@meta.data))
  seu_cd8$percent_mt <- PercentageFeatureSet(seu_cd8, pattern = "^mt-")

qc_cd8 <- seu_cd8@meta.data %>%
  dplyr::group_by(annot_cluster) %>%
  dplyr::summarise(n_cells       = dplyr::n(),
                   med_nFeature  = median(nFeature_RNA),
                   med_nCount    = median(nCount_RNA),
                   med_percentmt = round(median(percent_mt), 2),
                   .groups = "drop") %>%
  dplyr::mutate(preflagged = annot_cluster %in% flag_clusters_cd8)
write.csv(qc_cd8, file.path(cd8_man_tab, "qc_per_cluster.csv"), row.names = FALSE)
print(qc_cd8)

# ===== CD8 : 2. CURATED MARKER DOT PLOT + AVG-EXPRESSION HEATMAP =============
message("\n=== CD8 : marker dotplot + heatmap ===")
panel_cd8 <- unique(unlist(cd8_signatures))
panel_cd8 <- intersect(panel_cd8, rownames(seu_cd8))

p_dot8 <- DotPlot(seu_cd8, features = panel_cd8, cluster.idents = FALSE) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic", size = 9)) +
  labs(title = "CD8: curated markers x cluster", x = NULL, y = "Cluster")
ggsave(file.path(cd8_man, "marker_dotplot.png"), p_dot8,
       width = 10, height = 11, dpi = 300, bg = "white")

avg8 <- AverageExpression(seu_cd8, assays = "RNA", layer = "data",
                          features = panel_cd8, group.by = "annot_cluster")$RNA
colnames(avg8) <- sub("^g", "", colnames(avg8))        # Seurat prefixes numeric idents with "g"
avg8 <- avg8[, order(as.integer(colnames(avg8))), drop = FALSE]
avg8_log <- log1p(as.matrix(avg8))
pheatmap(avg8_log, scale = "row", cluster_cols = FALSE, cluster_rows = TRUE,
         fontsize_row = 8, main = "CD8: curated markers (row-scaled log avg expr)",
         filename = file.path(cd8_man, "marker_avgexpr_heatmap.png"),
         width = 9, height = 11)

# ===== CD8 : 3. SIGNATURE MODULE SCORES + HEATMAP ===========================
message("\n=== CD8 : signature scores ===")
sig_cols_cd8 <- character(0)
for (sn in names(cd8_signatures)) {
  feats <- intersect(cd8_signatures[[sn]], rownames(seu_cd8))
  if (length(feats) < 2) { message("  skip signature ", sn, " (<2 genes present)"); next }
  seu_cd8 <- AddModuleScore(seu_cd8, features = list(feats),
                            name = paste0("sig_", sn, "_"), seed = 42)
  sig_cols_cd8 <- c(sig_cols_cd8, paste0("sig_", sn, "_1"))
}
sigmat_cd8 <- seu_cd8@meta.data %>%
  dplyr::group_by(annot_cluster) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(sig_cols_cd8), mean), .groups = "drop") %>%
  tibble::column_to_rownames("annot_cluster") %>% as.matrix()
colnames(sigmat_cd8) <- sub("^sig_", "", sub("_1$", "", colnames(sigmat_cd8)))
write.csv(round(as.data.frame(sigmat_cd8), 4),
          file.path(cd8_man_tab, "signature_scores_per_cluster.csv"))
pheatmap(t(sigmat_cd8), scale = "row", cluster_cols = FALSE, cluster_rows = TRUE,
         display_numbers = TRUE, number_format = "%.2f",
         main = "CD8: signature module scores (row-scaled)",
         filename = file.path(cd8_man, "signature_heatmap.png"),
         width = 9, height = 6)

# ===== CD8 : 4. CLUSTER CORRELATION + DENDROGRAM + MERGE CANDIDATES =========
message("\n=== CD8 : merge candidates ===")
hvg8  <- intersect(VariableFeatures(seu_cd8), rownames(seu_cd8))
avgH8 <- AverageExpression(seu_cd8, assays = "RNA", layer = "data",
                           features = hvg8, group.by = "annot_cluster")$RNA
colnames(avgH8) <- sub("^g", "", colnames(avgH8))
avgH8 <- avgH8[, order(as.integer(colnames(avgH8))), drop = FALSE]
cor8  <- cor(log1p(as.matrix(avgH8)), method = "pearson")
pheatmap(cor8, main = "CD8: cluster-cluster correlation (HVG log avg expr)",
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         filename = file.path(cd8_man, "cluster_correlation_heatmap.png"),
         width = 8, height = 7)
hc8 <- hclust(as.dist(1 - cor8), method = "average")
png(file.path(cd8_man, "cluster_dendrogram.png"), width = 1600, height = 1000, res = 200)
plot(hc8, main = "CD8: cluster dendrogram (1 - correlation)", xlab = "", sub = "")
dev.off()

# Targeted pairwise DE on highly-correlated pairs -> merge candidates
pairs8 <- which(upper.tri(cor8), arr.ind = TRUE)
cand8  <- data.frame()
for (k in seq_len(nrow(pairs8))) {
  i <- rownames(cor8)[pairs8[k, 1]]; j <- colnames(cor8)[pairs8[k, 2]]
  r <- cor8[pairs8[k, 1], pairs8[k, 2]]
  if (r < merge_corr_threshold) next
  dem <- FindMarkers(seu_cd8, ident.1 = i, ident.2 = j, only.pos = FALSE,
                     min.pct = 0.1, logfc.threshold = 0.1, verbose = FALSE)
  n_de <- sum(dem$p_val_adj < de_padj & abs(dem$avg_log2FC) > de_lfc, na.rm = TRUE)
  cand8 <- rbind(cand8, data.frame(cluster_a = i, cluster_b = j,
                                   correlation = round(r, 3), n_DE_genes = n_de))
}
if (nrow(cand8) > 0)
  cand8 <- cand8[order(-cand8$correlation, cand8$n_DE_genes), ]
write.csv(cand8, file.path(cd8_man_tab, "merge_candidates.csv"), row.names = FALSE)
message("  ", nrow(cand8), " merge-candidate pair(s) (corr >= ", merge_corr_threshold, ")")

# ###################### PORTION B : AUTOMATIC ANNOTATION #####################
# ===== CD8 : 5. AUTOMATIC ANNOTATION -- ProjecTILs projection ================
# Full projection of the CD8 query onto the mouse LCMV CD8 atlas. Run.ProjecTILs
# wraps make.projection + cellstate.predict: it embeds the query in the reference
# UMAP and predicts a per-cell functional.cluster. We transfer that label back to
# seu_cd8 by barcode (NA for cells scGate filters as non-CD8/low-quality), then
# emit the tutorial diagnostics: projection UMAP, predicted-state composition,
# and a query-vs-reference marker radar. LCMV-derived ref -> cross-check only.
message("\n=== CD8 : ProjecTILs projection ===")
options(timeout = 3000)
if (!file.exists(projectils_ref_cd8) || file.info(projectils_ref_cd8)$size < 1e7)
  download.file(projectils_url_cd8, projectils_ref_cd8, mode = "wb", method = "libcurl")
if (file.info(projectils_ref_cd8)$size < 1e7)
  stop("CD8 reference looks like a redirect stub (",
       file.info(projectils_ref_cd8)$size, " bytes). Delete it and re-download ",
       "using the ndownloader.figshare.com URL form.")
ref8 <- load.reference.map(projectils_ref_cd8)
ref8 <- UpdateSeuratObject(ref8)   # atlas is an old (v3-era) object: add missing slots

# project (filter.cells = TRUE by default runs scGate to drop non-T/low-quality)
query8 <- Run.ProjecTILs(seu_cd8, ref = ref8, split.by = projectils_split_cd8)

# transfer predicted label back to the full object (NA = filtered out)
seu_cd8$functional.cluster <- NA_character_
seu_cd8$functional.cluster[colnames(query8)] <- as.character(query8$functional.cluster)
n_lab8 <- sum(!is.na(seu_cd8$functional.cluster))
message(sprintf("  projected %d / %d cells (%.1f%%); rest filtered as non-T/low-quality.",
                n_lab8, ncol(seu_cd8), 100 * n_lab8 / ncol(seu_cd8)))

# diagnostic 1: query projected over the reference UMAP. Built manually with
# OPAQUE points (reference grey underneath, query coloured on top) -- this
# avoids the alpha-blend moiré (full-width banding) that plot.projection's
# default transparent rendering produced. Coords come from the projected
# query, which carries the reference UMAP embedding.
red_q8 <- grep("umap", Reductions(query8),  ignore.case = TRUE, value = TRUE)[1]
red_r8 <- grep("umap", Reductions(ref8), ignore.case = TRUE, value = TRUE)[1]
refemb8 <- as.data.frame(Embeddings(ref8, red_r8))[, 1:2]
qemb8   <- as.data.frame(Embeddings(query8, red_q8))[, 1:2]
colnames(refemb8) <- colnames(qemb8) <- c("UMAP_1", "UMAP_2")
qemb8$state <- factor(query8$functional.cluster)
qemb8 <- qemb8[!is.na(qemb8$state), ]
pal8 <- ref8@misc$atlas.palette        # reference's own state colours
p_proj8 <- ggplot() +
  geom_point(data = refemb8, aes(UMAP_1, UMAP_2),
             colour = "grey85", size = 0.3, stroke = 0) +
  geom_point(data = qemb8, aes(UMAP_1, UMAP_2, colour = state),
             size = 0.5, stroke = 0) +
  (if (!is.null(pal8)) scale_colour_manual(values = pal8, na.translate = FALSE)
   else scale_colour_hue(na.translate = FALSE)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  ggtitle("CD8: query projected on LCMV reference") +
  theme_classic() + theme(aspect.ratio = 1)
ggsave(file.path(cd8_aut, "projectils_projection_umap.png"), p_proj8,
       width = 8, height = 6, dpi = 300, bg = "white")

# diagnostic 2: predicted reference-state composition
p_comp8 <- plot.statepred.composition(ref8, query8, metric = "Percent") +
  ggtitle("CD8: predicted reference-state composition")
ggsave(file.path(cd8_aut, "projectils_state_composition.png"), p_comp8,
       width = 7, height = 5, dpi = 300, bg = "white")

# diagnostic 3: per-state marker radar (query in colour vs reference in black)
p_radar8 <- plot.states.radar(ref8, query = query8, genes4radar = genes4radar_cd8,
                              min.cells = 20, return = TRUE)
ggsave(file.path(cd8_aut, "projectils_marker_radar.png"), p_radar8,
       width = 11, height = 9, dpi = 300, bg = "white")

# our clusters x ProjecTILs label crosstab (row-normalised), NA labels dropped
lab8 <- factor(seu_cd8$functional.cluster)
keep8 <- !is.na(lab8)
ct8  <- as.matrix(table(seu_cd8$annot_cluster[keep8], droplevels(lab8[keep8])))
ct8_prop <- sweep(ct8, 1, pmax(rowSums(ct8), 1), "/")
write.csv(as.data.frame.matrix(ct8), file.path(cd8_aut_tab, "projectils_crosstab_counts.csv"))
write.csv(round(as.data.frame.matrix(ct8_prop), 3),
          file.path(cd8_aut_tab, "projectils_crosstab_prop.csv"))
pheatmap(ct8_prop, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         main = "CD8: ProjecTILs label fraction per cluster",
         filename = file.path(cd8_aut, "projectils_heatmap.png"),
         width = 9, height = 6)

# ===== CD8 : 6. SUGGESTED-NAMES WORKSHEET + SAVE ============================
message("\n=== CD8 : suggested names ===")
top20_cd8 <- tryCatch(
  read.csv(file.path(sub_base, "CD8", "Markers",
                     sprintf("top20_per_cluster_res%.2f.csv", res_cd8))),
  error = function(e) NULL)
top_markers_cd8 <- function(cl) {
  if (is.null(top20_cd8)) return("")
  g <- top20_cd8$gene[as.character(top20_cd8$cluster) == as.character(cl)]
  paste(head(g, 6), collapse = ", ")
}
rank_sig <- function(row) names(sort(row, decreasing = TRUE))
worksheet_cd8 <- lapply(clusters_cd8, function(cl) {
  sc <- sigmat_cd8[as.character(cl), ]
  ord <- rank_sig(sc)
  pl <- seu_cd8$functional.cluster[seu_cd8$annot_cluster == cl]
  pl <- pl[!is.na(pl)]
  pmaj <- if (length(pl)) names(sort(table(pl), decreasing = TRUE))[1] else NA
  pfrac <- if (length(pl)) round(max(table(pl)) / length(pl), 2) else NA
  data.frame(
    cluster            = cl,
    n_cells            = sum(seu_cd8$annot_cluster == cl),
    top_signature      = ord[1],
    second_signature   = ord[2],
    projectils_majority = pmaj,
    projectils_frac    = pfrac,
    top_markers        = top_markers_cd8(cl),
    suggested_name     = unname(cd8_name_lookup[ord[1]]),
    qc_flag            = ifelse(cl %in% flag_clusters_cd8, "REVIEW (possible doublet/ambient)", ""),
    stringsAsFactors   = FALSE)
}) %>% dplyr::bind_rows()
write.csv(worksheet_cd8, file.path(cd8_dir, "suggested_names.csv"), row.names = FALSE)
print(worksheet_cd8[, c("cluster","top_signature","projectils_majority","suggested_name","qc_flag")])

# Finalised proposal for the biologist (built on the worksheet; see params block).
if (length(cd8_proposed_fine)) {
  prop_cd8 <- worksheet_cd8[, c("cluster","n_cells","projectils_majority","projectils_frac","top_markers")]
  prop_cd8$proposed_state <- unname(cd8_proposed_fine[as.character(prop_cd8$cluster)])
  prop_cd8$merged_group   <- unname(cd8_proposed_group[as.character(prop_cd8$cluster)])
  prop_cd8$action         <- ifelse(grepl("EXCLUDE", prop_cd8$proposed_state), "exclude", "keep")
  prop_cd8 <- prop_cd8[, c("cluster","n_cells","proposed_state","merged_group","action",
                           "projectils_majority","projectils_frac","top_markers")]
  write.csv(prop_cd8, file.path(cd8_dir, "proposed_annotation_CD8.csv"), row.names = FALSE)
  print(prop_cd8[, c("cluster","proposed_state","merged_group","action")])
}

writeLines(c(
  "CD8 annotation support (script 12) -- SUGGESTIONS ONLY, biologist decides.",
  paste0("Resolution: ", res_cd8, "  (", length(clusters_cd8), " clusters)"),
  "Layout:",
  "  manual/    marker_dotplot, marker_avgexpr_heatmap, signature_heatmap,",
  "             cluster_correlation_heatmap, cluster_dendrogram (PNGs)",
  "  manual/tables/    qc_per_cluster, signature_scores_per_cluster, merge_candidates",
  "  automatic/ projectils_projection_umap, projectils_state_composition,",
  "             projectils_marker_radar, projectils_heatmap (PNGs)",
  "  automatic/tables/ projectils_crosstab_counts, projectils_crosstab_prop",
  "  (top)      suggested_names.csv, proposed_annotation_CD8.csv, README.txt",
  "Pre-flagged for review (possible doublet/ambient): cluster(s) ",
  paste(flag_clusters_cd8, collapse = ", "),
  "Note: ProjecTILs refs are LCMV/TIL-derived; in-vitro pre-infusion cells fit imperfectly."
), file.path(cd8_dir, "README.txt"))

qs_save(seu_cd8, file.path(saved_dir, "Mouse_CARTmiR29a_CD8_annotated.qs2"))
message("  saved CD8 annotated object (scores + ProjecTILs labels; NO final names).")


###############################################################################
###############################################################################
##                              C D 4                                        ##
###############################################################################
###############################################################################

cd4_dir <- file.path(sub_base, "CD4", "annotation"); dir.create(cd4_dir, recursive = TRUE, showWarnings = FALSE)
# Organised outputs: manual/ and automatic/ portions, each with a tables/
# subdir for CSVs (PNGs sit directly in the portion folder). Synthesis
# files (suggested_names, proposed_annotation, README) live at the top.
cd4_man     <- file.path(cd4_dir, "manual");        dir.create(file.path(cd4_man, "tables"), recursive = TRUE, showWarnings = FALSE)
cd4_man_tab <- file.path(cd4_man, "tables")
cd4_aut     <- file.path(cd4_dir, "automatic");     dir.create(file.path(cd4_aut, "tables"), recursive = TRUE, showWarnings = FALSE)
cd4_aut_tab <- file.path(cd4_aut, "tables")

# ######################## PORTION A : MANUAL ANNOTATION ######################
# ===== CD4 : 1. LOAD + SET IDENT FROM COMMITTED RESOLUTION + QC TABLE ========
message("\n=== CD4 : load + QC ===")
seu_cd4 <- qs_read(file.path(saved_dir, "Mouse_CARTmiR29a_CD4_subclustered.qs2"))
DefaultAssay(seu_cd4) <- "RNA"

res_col_cd4 <- paste0("RNA_snn_res.", res_cd4)
stopifnot(res_col_cd4 %in% colnames(seu_cd4@meta.data))
seu_cd4$annot_cluster <- factor(seu_cd4[[res_col_cd4]][, 1],
                                levels = sort(as.integer(levels(seu_cd4[[res_col_cd4]][, 1]))))
Idents(seu_cd4) <- "annot_cluster"
clusters_cd4 <- levels(seu_cd4$annot_cluster)

if (!"percent_mt" %in% colnames(seu_cd4@meta.data))
  seu_cd4$percent_mt <- PercentageFeatureSet(seu_cd4, pattern = "^mt-")

qc_cd4 <- seu_cd4@meta.data %>%
  dplyr::group_by(annot_cluster) %>%
  dplyr::summarise(n_cells       = dplyr::n(),
                   med_nFeature  = median(nFeature_RNA),
                   med_nCount    = median(nCount_RNA),
                   med_percentmt = round(median(percent_mt), 2),
                   .groups = "drop") %>%
  dplyr::mutate(preflagged = annot_cluster %in% flag_clusters_cd4)
write.csv(qc_cd4, file.path(cd4_man_tab, "qc_per_cluster.csv"), row.names = FALSE)
print(qc_cd4)

# ===== CD4 : 2. CURATED MARKER DOT PLOT + AVG-EXPRESSION HEATMAP =============
message("\n=== CD4 : marker dotplot + heatmap ===")
panel_cd4 <- unique(unlist(cd4_signatures))
panel_cd4 <- intersect(panel_cd4, rownames(seu_cd4))

p_dot4 <- DotPlot(seu_cd4, features = panel_cd4, cluster.idents = FALSE) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic", size = 9)) +
  labs(title = "CD4: curated markers x cluster", x = NULL, y = "Cluster")
ggsave(file.path(cd4_man, "marker_dotplot.png"), p_dot4,
       width = 11, height = 11, dpi = 300, bg = "white")

avg4 <- AverageExpression(seu_cd4, assays = "RNA", layer = "data",
                          features = panel_cd4, group.by = "annot_cluster")$RNA
colnames(avg4) <- sub("^g", "", colnames(avg4))        # Seurat prefixes numeric idents with "g"
avg4 <- avg4[, order(as.integer(colnames(avg4))), drop = FALSE]
avg4_log <- log1p(as.matrix(avg4))
pheatmap(avg4_log, scale = "row", cluster_cols = FALSE, cluster_rows = TRUE,
         fontsize_row = 8, main = "CD4: curated markers (row-scaled log avg expr)",
         filename = file.path(cd4_man, "marker_avgexpr_heatmap.png"),
         width = 10, height = 11)

# ===== CD4 : 3. SIGNATURE MODULE SCORES + HEATMAP ===========================
message("\n=== CD4 : signature scores ===")
sig_cols_cd4 <- character(0)
for (sn in names(cd4_signatures)) {
  feats <- intersect(cd4_signatures[[sn]], rownames(seu_cd4))
  if (length(feats) < 2) { message("  skip signature ", sn, " (<2 genes present)"); next }
  seu_cd4 <- AddModuleScore(seu_cd4, features = list(feats),
                            name = paste0("sig_", sn, "_"), seed = 42)
  sig_cols_cd4 <- c(sig_cols_cd4, paste0("sig_", sn, "_1"))
}
sigmat_cd4 <- seu_cd4@meta.data %>%
  dplyr::group_by(annot_cluster) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(sig_cols_cd4), mean), .groups = "drop") %>%
  tibble::column_to_rownames("annot_cluster") %>% as.matrix()
colnames(sigmat_cd4) <- sub("^sig_", "", sub("_1$", "", colnames(sigmat_cd4)))
write.csv(round(as.data.frame(sigmat_cd4), 4),
          file.path(cd4_man_tab, "signature_scores_per_cluster.csv"))
pheatmap(t(sigmat_cd4), scale = "row", cluster_cols = FALSE, cluster_rows = TRUE,
         display_numbers = TRUE, number_format = "%.2f",
         main = "CD4: signature module scores (row-scaled)",
         filename = file.path(cd4_man, "signature_heatmap.png"),
         width = 10, height = 7)

# ===== CD4 : 4. CLUSTER CORRELATION + DENDROGRAM + MERGE CANDIDATES =========
message("\n=== CD4 : merge candidates ===")
hvg4  <- intersect(VariableFeatures(seu_cd4), rownames(seu_cd4))
avgH4 <- AverageExpression(seu_cd4, assays = "RNA", layer = "data",
                           features = hvg4, group.by = "annot_cluster")$RNA
colnames(avgH4) <- sub("^g", "", colnames(avgH4))
avgH4 <- avgH4[, order(as.integer(colnames(avgH4))), drop = FALSE]
cor4  <- cor(log1p(as.matrix(avgH4)), method = "pearson")
pheatmap(cor4, main = "CD4: cluster-cluster correlation (HVG log avg expr)",
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         filename = file.path(cd4_man, "cluster_correlation_heatmap.png"),
         width = 9, height = 8)
hc4 <- hclust(as.dist(1 - cor4), method = "average")
png(file.path(cd4_man, "cluster_dendrogram.png"), width = 1700, height = 1000, res = 200)
plot(hc4, main = "CD4: cluster dendrogram (1 - correlation)", xlab = "", sub = "")
dev.off()

pairs4 <- which(upper.tri(cor4), arr.ind = TRUE)
cand4  <- data.frame()
for (k in seq_len(nrow(pairs4))) {
  i <- rownames(cor4)[pairs4[k, 1]]; j <- colnames(cor4)[pairs4[k, 2]]
  r <- cor4[pairs4[k, 1], pairs4[k, 2]]
  if (r < merge_corr_threshold) next
  dem <- FindMarkers(seu_cd4, ident.1 = i, ident.2 = j, only.pos = FALSE,
                     min.pct = 0.1, logfc.threshold = 0.1, verbose = FALSE)
  n_de <- sum(dem$p_val_adj < de_padj & abs(dem$avg_log2FC) > de_lfc, na.rm = TRUE)
  cand4 <- rbind(cand4, data.frame(cluster_a = i, cluster_b = j,
                                   correlation = round(r, 3), n_DE_genes = n_de))
}
if (nrow(cand4) > 0)
  cand4 <- cand4[order(-cand4$correlation, cand4$n_DE_genes), ]
write.csv(cand4, file.path(cd4_man_tab, "merge_candidates.csv"), row.names = FALSE)
message("  ", nrow(cand4), " merge-candidate pair(s) (corr >= ", merge_corr_threshold, ")")

# ###################### PORTION B : AUTOMATIC ANNOTATION #####################
# ===== CD4 : 5. AUTOMATIC ANNOTATION -- ProjecTILs projection ================
# Full projection of the CD4 query onto the mouse LCMV CD4 atlas. Run.ProjecTILs
# wraps make.projection + cellstate.predict: it embeds the query in the reference
# UMAP and predicts a per-cell functional.cluster. We transfer that label back to
# seu_cd4 by barcode (NA for cells scGate filters as non-CD4/low-quality), then
# emit the tutorial diagnostics: projection UMAP, predicted-state composition,
# and a query-vs-reference marker radar. LCMV-derived ref -> cross-check only.
message("\n=== CD4 : ProjecTILs projection ===")
options(timeout = 3000)
if (!file.exists(projectils_ref_cd4) || file.info(projectils_ref_cd4)$size < 1e7)
  download.file(projectils_url_cd4, projectils_ref_cd4, mode = "wb", method = "libcurl")
if (file.info(projectils_ref_cd4)$size < 1e7)
  stop("CD4 reference looks like a redirect stub (",
       file.info(projectils_ref_cd4)$size, " bytes). Delete it and re-download ",
       "using the ndownloader.figshare.com URL form.")
ref4 <- load.reference.map(projectils_ref_cd4)
ref4 <- UpdateSeuratObject(ref4)   # atlas is an old (v3-era) object: add missing slots

# project (filter.cells = TRUE by default runs scGate to drop non-T/low-quality)
query4 <- Run.ProjecTILs(seu_cd4, ref = ref4, split.by = projectils_split_cd4)

# transfer predicted label back to the full object (NA = filtered out)
seu_cd4$functional.cluster <- NA_character_
seu_cd4$functional.cluster[colnames(query4)] <- as.character(query4$functional.cluster)
n_lab4 <- sum(!is.na(seu_cd4$functional.cluster))
message(sprintf("  projected %d / %d cells (%.1f%%); rest filtered as non-T/low-quality.",
                n_lab4, ncol(seu_cd4), 100 * n_lab4 / ncol(seu_cd4)))

# diagnostic 1: query projected over the reference UMAP. Built manually with
# OPAQUE points (reference grey underneath, query coloured on top) -- this
# avoids the alpha-blend moiré (full-width banding) that plot.projection's
# default transparent rendering produced. Coords come from the projected
# query, which carries the reference UMAP embedding.
red_q4 <- grep("umap", Reductions(query4),  ignore.case = TRUE, value = TRUE)[1]
red_r4 <- grep("umap", Reductions(ref4), ignore.case = TRUE, value = TRUE)[1]
refemb4 <- as.data.frame(Embeddings(ref4, red_r4))[, 1:2]
qemb4   <- as.data.frame(Embeddings(query4, red_q4))[, 1:2]
colnames(refemb4) <- colnames(qemb4) <- c("UMAP_1", "UMAP_2")
qemb4$state <- factor(query4$functional.cluster)
qemb4 <- qemb4[!is.na(qemb4$state), ]
pal4 <- ref4@misc$atlas.palette        # reference's own state colours
p_proj4 <- ggplot() +
  geom_point(data = refemb4, aes(UMAP_1, UMAP_2),
             colour = "grey85", size = 0.3, stroke = 0) +
  geom_point(data = qemb4, aes(UMAP_1, UMAP_2, colour = state),
             size = 0.5, stroke = 0) +
  (if (!is.null(pal4)) scale_colour_manual(values = pal4, na.translate = FALSE)
   else scale_colour_hue(na.translate = FALSE)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  ggtitle("CD4: query projected on LCMV reference") +
  theme_classic() + theme(aspect.ratio = 1)
ggsave(file.path(cd4_aut, "projectils_projection_umap.png"), p_proj4,
       width = 8, height = 6, dpi = 300, bg = "white")

# diagnostic 2: predicted reference-state composition
p_comp4 <- plot.statepred.composition(ref4, query4, metric = "Percent") +
  ggtitle("CD4: predicted reference-state composition")
ggsave(file.path(cd4_aut, "projectils_state_composition.png"), p_comp4,
       width = 7, height = 5, dpi = 300, bg = "white")

# diagnostic 3: per-state marker radar (query in colour vs reference in black)
p_radar4 <- plot.states.radar(ref4, query = query4, genes4radar = genes4radar_cd4,
                              min.cells = 20, return = TRUE)
ggsave(file.path(cd4_aut, "projectils_marker_radar.png"), p_radar4,
       width = 11, height = 9, dpi = 300, bg = "white")

# our clusters x ProjecTILs label crosstab (row-normalised), NA labels dropped
lab4 <- factor(seu_cd4$functional.cluster)
keep4 <- !is.na(lab4)
ct4  <- as.matrix(table(seu_cd4$annot_cluster[keep4], droplevels(lab4[keep4])))
ct4_prop <- sweep(ct4, 1, pmax(rowSums(ct4), 1), "/")
write.csv(as.data.frame.matrix(ct4), file.path(cd4_aut_tab, "projectils_crosstab_counts.csv"))
write.csv(round(as.data.frame.matrix(ct4_prop), 3),
          file.path(cd4_aut_tab, "projectils_crosstab_prop.csv"))
pheatmap(ct4_prop, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         main = "CD4: ProjecTILs label fraction per cluster",
         filename = file.path(cd4_aut, "projectils_heatmap.png"),
         width = 9, height = 6)

# ===== CD4 : 6. SUGGESTED-NAMES WORKSHEET + SAVE ============================
message("\n=== CD4 : suggested names ===")
top20_cd4 <- tryCatch(
  read.csv(file.path(sub_base, "CD4", "Markers",
                     sprintf("top20_per_cluster_res%.2f.csv", res_cd4))),
  error = function(e) NULL)
top_markers_cd4 <- function(cl) {
  if (is.null(top20_cd4)) return("")
  g <- top20_cd4$gene[as.character(top20_cd4$cluster) == as.character(cl)]
  paste(head(g, 6), collapse = ", ")
}
worksheet_cd4 <- lapply(clusters_cd4, function(cl) {
  sc <- sigmat_cd4[as.character(cl), ]
  ord <- names(sort(sc, decreasing = TRUE))
  pl <- seu_cd4$functional.cluster[seu_cd4$annot_cluster == cl]
  pl <- pl[!is.na(pl)]
  pmaj <- if (length(pl)) names(sort(table(pl), decreasing = TRUE))[1] else NA
  pfrac <- if (length(pl)) round(max(table(pl)) / length(pl), 2) else NA
  data.frame(
    cluster            = cl,
    n_cells            = sum(seu_cd4$annot_cluster == cl),
    top_signature      = ord[1],
    second_signature   = ord[2],
    projectils_majority = pmaj,
    projectils_frac    = pfrac,
    top_markers        = top_markers_cd4(cl),
    suggested_name     = unname(cd4_name_lookup[ord[1]]),
    qc_flag            = ifelse(cl %in% flag_clusters_cd4, "REVIEW (possible doublet/ambient)", ""),
    stringsAsFactors   = FALSE)
}) %>% dplyr::bind_rows()
write.csv(worksheet_cd4, file.path(cd4_dir, "suggested_names.csv"), row.names = FALSE)
print(worksheet_cd4[, c("cluster","top_signature","projectils_majority","suggested_name","qc_flag")])

writeLines(c(
  "CD4 annotation support (script 12) -- SUGGESTIONS ONLY, biologist decides.",
  paste0("Resolution: ", res_cd4, "  (", length(clusters_cd4), " clusters)"),
  "Layout:",
  "  manual/    marker_dotplot, marker_avgexpr_heatmap, signature_heatmap,",
  "             cluster_correlation_heatmap, cluster_dendrogram (PNGs)",
  "  manual/tables/    qc_per_cluster, signature_scores_per_cluster, merge_candidates",
  "  automatic/ projectils_projection_umap, projectils_state_composition,",
  "             projectils_marker_radar, projectils_heatmap (PNGs)",
  "  automatic/tables/ projectils_crosstab_counts, projectils_crosstab_prop",
  "  (top)      suggested_names.csv, README.txt",
  "Pre-flagged for review (possible doublet/ambient): cluster(s) ",
  paste(flag_clusters_cd4, collapse = ", "),
  "Note: ProjecTILs refs are LCMV/TIL-derived; in-vitro pre-infusion cells fit imperfectly."
), file.path(cd4_dir, "README.txt"))

qs_save(seu_cd4, file.path(saved_dir, "Mouse_CARTmiR29a_CD4_annotated.qs2"))
message("  saved CD4 annotated object (scores + ProjecTILs labels; NO final names).")


# ============================================================================
message("\nDone. Annotation-support package in: ", file.path(sub_base, "<lineage>", "annotation"))
message("Hand suggested_names.csv + merge_candidates.csv + the heatmaps to the")
message("biologist. Final labels are theirs to set; a later script will apply the")
message("agreed name map to annot_cluster and write the definitive labels.")
###############################################################################

