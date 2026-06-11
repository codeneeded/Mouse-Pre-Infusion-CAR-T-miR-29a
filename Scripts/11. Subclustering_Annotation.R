###############################################################################
# 12_annotation_support_CD8_CD4_mouse_CARTmiR29a.R
#
# ANNOTATION DECISION-SUPPORT for the CD8 (res 0.5) and CD4 (res 0.6)
# subclustered objects from script 11. This script does NOT assign final
# labels -- it produces evidence + suggestions for the wetlab biologist, who
# makes the final call. Specifically, per lineage it generates:
#
#   1. AUTOMATIC annotation        -- ProjecTILs projection onto mouse CD8/CD4
#                                     reference atlases (optional; guarded so
#                                     the rest runs without it). Reference-based
#                                     and EXPECTED to fit imperfectly for
#                                     in-vitro pre-infusion cells -- a cross-
#                                     check, not ground truth.
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
#   Annotation/<lineage>/  (QC, dotplot, heatmaps, dendrogram, merge + name CSVs)
#   saved_R_data/Mouse_CARTmiR29a_<lineage>_annotated.qs2  (carries scores +
#     ProjecTILs labels + the chosen-resolution ident as annot_cluster; NO final
#     names -- those are added in a later step once the biologist decides)
###############################################################################

library(Seurat)
library(SeuratExtend)
library(scCustomize)
library(qs2)
library(dplyr); library(tidyr); library(tibble)
library(ggplot2); library(patchwork)
library(pheatmap)
# ProjecTILs is optional -- loaded inside the guarded section if installed.

# ============================ Paths ==========================================
project_dir <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/10x_scRNAseq"
saved_dir   <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/saved_R_data"
sub_base    <- file.path(project_dir, "Subclustering")
annot_base  <- file.path(project_dir, "Annotation")
dir.create(annot_base, recursive = TRUE, showWarnings = FALSE)

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

# ProjecTILs reference maps (mouse). Download from the Carmona lab / SPICA
# (https://spica.unil.ch/refs) or Zenodo and point these at the .rds files.
# Leave as NA / nonexistent paths to SKIP auto-annotation gracefully.
projectils_ref_cd8 <- file.path(saved_dir, "ref_LCMV_CD8_mouse_v1.rds")
projectils_ref_cd4 <- file.path(saved_dir, "ref_LCMV_CD4_mouse_v1.rds")

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

cd8_dir <- file.path(annot_base, "CD8"); dir.create(cd8_dir, recursive = TRUE, showWarnings = FALSE)

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
write.csv(qc_cd8, file.path(cd8_dir, "qc_per_cluster.csv"), row.names = FALSE)
print(qc_cd8)

# ===== CD8 : 2. AUTOMATIC ANNOTATION -- ProjecTILs (guarded/optional) ========
message("\n=== CD8 : ProjecTILs (auto) ===")
seu_cd8$functional.cluster <- NA_character_
if (requireNamespace("ProjecTILs", quietly = TRUE) && file.exists(projectils_ref_cd8)) {
  ref8 <- ProjecTILs::load.reference.map(projectils_ref_cd8)
  seu_cd8 <- ProjecTILs::ProjecTILs.classifier(query = seu_cd8, ref = ref8,
                                               split.by = NULL)
  # cluster x predicted-label crosstab (row-normalised)
  ct8 <- as.matrix(table(seu_cd8$annot_cluster,
                         factor(seu_cd8$functional.cluster)))
  ct8_prop <- sweep(ct8, 1, pmax(rowSums(ct8), 1), "/")
  write.csv(as.data.frame.matrix(ct8), file.path(cd8_dir, "projectils_crosstab_counts.csv"))
  write.csv(round(as.data.frame.matrix(ct8_prop), 3),
            file.path(cd8_dir, "projectils_crosstab_prop.csv"))
  pheatmap(ct8_prop, cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.2f",
           main = "CD8: ProjecTILs label fraction per cluster",
           filename = file.path(cd8_dir, "projectils_heatmap.png"),
           width = 9, height = 6)
} else {
  message("  ProjecTILs unavailable or reference missing -- skipping auto-annotation.")
  message("  (install carmonalab/ProjecTILs and set projectils_ref_cd8 to enable.)")
}

# ===== CD8 : 3. CURATED MARKER DOT PLOT + AVG-EXPRESSION HEATMAP =============
message("\n=== CD8 : marker dotplot + heatmap ===")
panel_cd8 <- unique(unlist(cd8_signatures))
panel_cd8 <- intersect(panel_cd8, rownames(seu_cd8))

p_dot8 <- DotPlot(seu_cd8, features = panel_cd8, cluster.idents = TRUE) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic", size = 9)) +
  labs(title = "CD8: curated markers x cluster", x = NULL, y = "Cluster")
ggsave(file.path(cd8_dir, "marker_dotplot.png"), p_dot8,
       width = 10, height = 11, dpi = 300, bg = "white")

avg8 <- AverageExpression(seu_cd8, assays = "RNA", layer = "data",
                          features = panel_cd8, group.by = "annot_cluster")$RNA
avg8_log <- log1p(as.matrix(avg8))
pheatmap(avg8_log, scale = "row", cluster_cols = TRUE, cluster_rows = TRUE,
         fontsize_row = 8, main = "CD8: curated markers (row-scaled log avg expr)",
         filename = file.path(cd8_dir, "marker_avgexpr_heatmap.png"),
         width = 9, height = 11)

# ===== CD8 : 4. SIGNATURE MODULE SCORES + HEATMAP ===========================
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
          file.path(cd8_dir, "signature_scores_per_cluster.csv"))
pheatmap(t(sigmat_cd8), scale = "row", cluster_cols = TRUE, cluster_rows = TRUE,
         display_numbers = TRUE, number_format = "%.2f",
         main = "CD8: signature module scores (row-scaled)",
         filename = file.path(cd8_dir, "signature_heatmap.png"),
         width = 9, height = 6)

# ===== CD8 : 5. CLUSTER CORRELATION + DENDROGRAM + MERGE CANDIDATES =========
message("\n=== CD8 : merge candidates ===")
hvg8  <- intersect(VariableFeatures(seu_cd8), rownames(seu_cd8))
avgH8 <- AverageExpression(seu_cd8, assays = "RNA", layer = "data",
                           features = hvg8, group.by = "annot_cluster")$RNA
cor8  <- cor(log1p(as.matrix(avgH8)), method = "pearson")
pheatmap(cor8, main = "CD8: cluster-cluster correlation (HVG log avg expr)",
         display_numbers = TRUE, number_format = "%.2f",
         filename = file.path(cd8_dir, "cluster_correlation_heatmap.png"),
         width = 8, height = 7)
hc8 <- hclust(as.dist(1 - cor8), method = "average")
png(file.path(cd8_dir, "cluster_dendrogram.png"), width = 1600, height = 1000, res = 200)
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
write.csv(cand8, file.path(cd8_dir, "merge_candidates.csv"), row.names = FALSE)
message("  ", nrow(cand8), " merge-candidate pair(s) (corr >= ", merge_corr_threshold, ")")

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

writeLines(c(
  "CD8 annotation support (script 12).",
  paste0("Resolution: ", res_cd8, "  (", length(clusters_cd8), " clusters)"),
  "Files: qc_per_cluster, projectils_* (if run), marker_dotplot, marker_avgexpr_heatmap,",
  "       signature_heatmap, signature_scores_per_cluster, cluster_correlation_heatmap,",
  "       cluster_dendrogram, merge_candidates, suggested_names.",
  "Pre-flagged for review (possible doublet/ambient): cluster(s) " ,
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

cd4_dir <- file.path(annot_base, "CD4"); dir.create(cd4_dir, recursive = TRUE, showWarnings = FALSE)

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
write.csv(qc_cd4, file.path(cd4_dir, "qc_per_cluster.csv"), row.names = FALSE)
print(qc_cd4)

# ===== CD4 : 2. AUTOMATIC ANNOTATION -- ProjecTILs (guarded/optional) ========
message("\n=== CD4 : ProjecTILs (auto) ===")
seu_cd4$functional.cluster <- NA_character_
if (requireNamespace("ProjecTILs", quietly = TRUE) && file.exists(projectils_ref_cd4)) {
  ref4 <- ProjecTILs::load.reference.map(projectils_ref_cd4)
  seu_cd4 <- ProjecTILs::ProjecTILs.classifier(query = seu_cd4, ref = ref4,
                                               split.by = NULL)
  ct4 <- as.matrix(table(seu_cd4$annot_cluster,
                         factor(seu_cd4$functional.cluster)))
  ct4_prop <- sweep(ct4, 1, pmax(rowSums(ct4), 1), "/")
  write.csv(as.data.frame.matrix(ct4), file.path(cd4_dir, "projectils_crosstab_counts.csv"))
  write.csv(round(as.data.frame.matrix(ct4_prop), 3),
            file.path(cd4_dir, "projectils_crosstab_prop.csv"))
  pheatmap(ct4_prop, cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.2f",
           main = "CD4: ProjecTILs label fraction per cluster",
           filename = file.path(cd4_dir, "projectils_heatmap.png"),
           width = 9, height = 6)
} else {
  message("  ProjecTILs unavailable or reference missing -- skipping auto-annotation.")
  message("  (install carmonalab/ProjecTILs and set projectils_ref_cd4 to enable.)")
}

# ===== CD4 : 3. CURATED MARKER DOT PLOT + AVG-EXPRESSION HEATMAP =============
message("\n=== CD4 : marker dotplot + heatmap ===")
panel_cd4 <- unique(unlist(cd4_signatures))
panel_cd4 <- intersect(panel_cd4, rownames(seu_cd4))

p_dot4 <- DotPlot(seu_cd4, features = panel_cd4, cluster.idents = TRUE) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic", size = 9)) +
  labs(title = "CD4: curated markers x cluster", x = NULL, y = "Cluster")
ggsave(file.path(cd4_dir, "marker_dotplot.png"), p_dot4,
       width = 11, height = 11, dpi = 300, bg = "white")

avg4 <- AverageExpression(seu_cd4, assays = "RNA", layer = "data",
                          features = panel_cd4, group.by = "annot_cluster")$RNA
avg4_log <- log1p(as.matrix(avg4))
pheatmap(avg4_log, scale = "row", cluster_cols = TRUE, cluster_rows = TRUE,
         fontsize_row = 8, main = "CD4: curated markers (row-scaled log avg expr)",
         filename = file.path(cd4_dir, "marker_avgexpr_heatmap.png"),
         width = 10, height = 11)

# ===== CD4 : 4. SIGNATURE MODULE SCORES + HEATMAP ===========================
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
          file.path(cd4_dir, "signature_scores_per_cluster.csv"))
pheatmap(t(sigmat_cd4), scale = "row", cluster_cols = TRUE, cluster_rows = TRUE,
         display_numbers = TRUE, number_format = "%.2f",
         main = "CD4: signature module scores (row-scaled)",
         filename = file.path(cd4_dir, "signature_heatmap.png"),
         width = 10, height = 7)

# ===== CD4 : 5. CLUSTER CORRELATION + DENDROGRAM + MERGE CANDIDATES =========
message("\n=== CD4 : merge candidates ===")
hvg4  <- intersect(VariableFeatures(seu_cd4), rownames(seu_cd4))
avgH4 <- AverageExpression(seu_cd4, assays = "RNA", layer = "data",
                           features = hvg4, group.by = "annot_cluster")$RNA
cor4  <- cor(log1p(as.matrix(avgH4)), method = "pearson")
pheatmap(cor4, main = "CD4: cluster-cluster correlation (HVG log avg expr)",
         display_numbers = TRUE, number_format = "%.2f",
         filename = file.path(cd4_dir, "cluster_correlation_heatmap.png"),
         width = 9, height = 8)
hc4 <- hclust(as.dist(1 - cor4), method = "average")
png(file.path(cd4_dir, "cluster_dendrogram.png"), width = 1700, height = 1000, res = 200)
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
write.csv(cand4, file.path(cd4_dir, "merge_candidates.csv"), row.names = FALSE)
message("  ", nrow(cand4), " merge-candidate pair(s) (corr >= ", merge_corr_threshold, ")")

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
  "CD4 annotation support (script 12).",
  paste0("Resolution: ", res_cd4, "  (", length(clusters_cd4), " clusters)"),
  "Files: qc_per_cluster, projectils_* (if run), marker_dotplot, marker_avgexpr_heatmap,",
  "       signature_heatmap, signature_scores_per_cluster, cluster_correlation_heatmap,",
  "       cluster_dendrogram, merge_candidates, suggested_names.",
  "Pre-flagged for review (possible doublet/ambient): cluster(s) ",
  paste(flag_clusters_cd4, collapse = ", "),
  "Note: ProjecTILs refs are LCMV/TIL-derived; in-vitro pre-infusion cells fit imperfectly."
), file.path(cd4_dir, "README.txt"))

qs_save(seu_cd4, file.path(saved_dir, "Mouse_CARTmiR29a_CD4_annotated.qs2"))
message("  saved CD4 annotated object (scores + ProjecTILs labels; NO final names).")


# ============================================================================
message("\nDone. Annotation-support package in: ", annot_base)
message("Hand suggested_names.csv + merge_candidates.csv + the heatmaps to the")
message("biologist. Final labels are theirs to set; a later script will apply the")
message("agreed name map to annot_cluster and write the definitive labels.")
###############################################################################