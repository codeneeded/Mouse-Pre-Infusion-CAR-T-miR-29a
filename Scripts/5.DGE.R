###############################################################################
# 05 -- Differential Expression (Pseudobulk DESeq2 + MAST)                    #
# Mouse pre-infusion CAR-T scRNA-seq, miR-29a project                         #
###############################################################################
# Conditions and controls:
#   EV     = empty vector (no insert)         -- secondary control
#   Scr    = scramble miRNA (matched backbone) -- PRIMARY control (specificity)
#   miR29a = miR-29a construct
#
# Contrasts:
#   miR29a vs Scr  -- PRIMARY  : miR-29a-specific effect (isolates sequence)
#   miR29a vs EV   -- secondary: confirm direction with EV control
#   EV vs Scr      -- QC       : should be minimal if Scr is well-behaved
#
# Two DE methods run side-by-side:
#
#   1. PSEUDOBULK (DESeq2)  -- output: Differential_Expression/DGE_pseudobulk/
#      Aggregate counts per (condition x replicate), DESeq2 with
#      design ~ replicate + condition, lfcShrink (ashr) for LFC reliability.
#      CALIBRATED at low n (information-sharing across genes via empirical
#      Bayes dispersion fit). Defensible to a careful reviewer. This is the
#      inferentially-honest hit list and should be the manuscript primary.
#
#   2. MAST (FindMarkers)   -- output: Differential_Expression/DGE_MAST/
#      Single-cell-level DE with nCount_RNA as latent variable (OPIS pattern;
#      field convention). Inflates significance at n=2 reps because cells
#      from the same mouse are non-independent; treat the p-values as
#      EXPLORATORY. Useful as a broader screen and for comparability with
#      OPIS-style analyses.
#
# Method comparison output:
#   - DGE_summary_counts.csv / .png : DE counts per compartment x contrast x method
#   - DGE_method_overlap.csv        : sig-gene overlap PB vs MAST per group/contrast
#     (high overlap on strong hits = methods agree on what's real;
#      MAST-only = likely inflated; PB-only = rarer, calibration disagreement
#      on a real signal.)
#
# Compartments (both methods):
#   - by LINEAGE (CD4, CD8)        -- wetlab-requested (CD8 must be separate
#                                      from CD4 due to 2.7:1 skew)
#   - by CLUSTER (excluding 8, 13) -- secondary
#
# Output CSV format matches the OPIS pathway-script pattern:
#   first column = gene symbol (used as rownames on read);
#   p_val_adj, avg_log2FC, p_val columns. Same format for both methods.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(DESeq2)
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(stringr)
})

# ============================ Paths ==========================================
project_dir <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a"
saved_dir   <- file.path(project_dir, "saved_R_data")
de_base     <- file.path(project_dir, "Differential_Expression")

# Pseudobulk (PRIMARY -- calibrated at low n)
dir_PB_lin <- file.path(de_base, "DGE_pseudobulk", "by_lineage")
dir_PB_cl  <- file.path(de_base, "DGE_pseudobulk", "by_cluster")
# MAST (SECONDARY / exploratory -- field-default, but inflates significance
#       at n=2 reps; we run it for comparability and as a broader screen)
dir_MAST_lin <- file.path(de_base, "DGE_MAST", "by_lineage")
dir_MAST_cl  <- file.path(de_base, "DGE_MAST", "by_cluster")

for (d in c(dir_PB_lin, dir_PB_cl, dir_MAST_lin, dir_MAST_cl))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ============================ Load ===========================================
obj <- qs_read(file.path(saved_dir, "Mouse_CARTmiR29a_PreAnnotation.qs2"))
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")     # idempotent; safe even if already joined

# ensure log-normalized data layer exists (MAST needs it)
if (!"data" %in% Layers(obj[["RNA"]])) {
  message("NormalizeData(RNA) -- needed for MAST"); obj <- NormalizeData(obj, verbose = FALSE)
}

# All clusters and lineages are tested. The `exclude_cluster` flag remains in
# metadata as wetlab provenance but is NOT used to drop clusters from DE --
# the inferential output stands or falls on its own merits per cluster.
obj$clusters          <- droplevels(factor(obj$clusters))
obj$tentative_lineage <- droplevels(factor(obj$tentative_lineage))
obj$tentative_state   <- droplevels(factor(obj$tentative_state))

obj$condition <- factor(obj$condition, levels = c("EV","Scr","miR29a"))
obj$replicate <- factor(obj$replicate)

# pseudobulk sample key: condition-replicate
# (Seurat rewrites underscores in identity values to dashes, so we use dashes
# directly to avoid the round-trip mangling; condition and replicate values
# contain no dashes themselves.)
obj$.sample <- paste(obj$condition, obj$replicate, sep = "-")

# ============================ Contrasts ======================================
# c(factor_col, ident.1, ident.2) -- log2FC sign is ident.1 - ident.2
contrasts_list <- list(
  miR29a_vs_Scr = c("condition", "miR29a", "Scr"),   # PRIMARY
  miR29a_vs_EV  = c("condition", "miR29a", "EV"),    # secondary
  EV_vs_Scr     = c("condition", "EV",     "Scr")    # QC
)

# ============================ Helpers ========================================
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# Build cluster -> safe stem map: "NN_State_Label" (e.g. "00_Activated_intermediate_CD4")
# Zero-padded so directories/files sort numerically (10, 11, 12 don't land before 2).
cluster_label_map <- tapply(as.character(obj$tentative_state),
                            obj$clusters, function(x) x[1])
cluster_stem <- function(cl_id) {
  lbl <- cluster_label_map[as.character(cl_id)]
  if (is.na(lbl) || nchar(lbl) == 0) return(paste0("cluster_", cl_id))
  safe_name(paste0(sprintf("%02d", as.integer(as.character(cl_id))), "_", lbl))
}
# pretty version for plot labels (underscores -> spaces)
cluster_pretty <- function(cl_id) gsub("_", " ", cluster_stem(cl_id))

# Aggregate a Seurat subset into pseudobulk counts + sample-level metadata
make_pb <- function(seu_subset, min_cells_per_sample = 10) {
  n_tbl <- seu_subset@meta.data %>%
    dplyr::count(condition, replicate, name = "n_cells") %>%
    dplyr::filter(n_cells >= min_cells_per_sample)
  if (nrow(n_tbl) < 4) return(NULL)   # need >=4 pseudobulks for any contrast
  
  agg <- AggregateExpression(seu_subset, assays = "RNA", slot = "counts",
                             group.by = ".sample", return.seurat = FALSE)$RNA
  colnames(agg) <- sub("^g", "", colnames(agg))  # defensive (numeric-like levels)
  
  parts <- strsplit(colnames(agg), "-", fixed = TRUE)
  cm <- data.frame(
    sample_id = colnames(agg),
    condition = sapply(parts, `[`, 1),
    replicate = sapply(parts, `[`, 2),
    row.names = colnames(agg),
    stringsAsFactors = FALSE
  )
  cm$n_cells <- n_tbl$n_cells[
    match(paste(cm$condition, cm$replicate),
          paste(n_tbl$condition, n_tbl$replicate))
  ]
  keep <- !is.na(cm$n_cells) & cm$n_cells >= min_cells_per_sample
  if (sum(keep) < 4) return(NULL)
  agg <- agg[, keep, drop = FALSE]
  cm  <- cm[keep, , drop = FALSE]
  cm$condition <- factor(cm$condition, levels = c("EV","Scr","miR29a"))
  cm$replicate <- factor(cm$replicate)
  list(counts = as.matrix(agg), col_meta = cm)
}

# Run DESeq2 on one contrast within a pseudobulk; write CSV in OPIS-compatible format
run_de <- function(pb_counts, pb_meta, contrast, out_csv,
                   min_count = 10, min_samples = 4,
                   shrink_type = "ashr") {
  
  if (!all(contrast[2:3] %in% as.character(pb_meta$condition))) {
    message("    -> skip: contrast levels not both present"); return(invisible(NULL))
  }
  # restrict to just the two contrast groups
  keep <- pb_meta$condition %in% contrast[2:3]
  pb_meta <- droplevels(pb_meta[keep, , drop = FALSE])
  pb_counts <- pb_counts[, rownames(pb_meta), drop = FALSE]
  
  # gene filter: >= min_count in >= min_samples
  ok_g <- rowSums(pb_counts >= min_count) >= min_samples
  if (sum(ok_g) < 100) { message("    -> skip: <100 genes pass count filter"); return(invisible(NULL)) }
  pb_counts <- pb_counts[ok_g, , drop = FALSE]
  
  # design: include replicate if >1 level remain
  n_rep <- nlevels(droplevels(pb_meta$replicate))
  design_formula <- if (n_rep > 1) ~ replicate + condition else ~ condition
  
  dds <- DESeqDataSetFromMatrix(countData = pb_counts,
                                colData   = pb_meta,
                                design    = design_formula)
  dds <- DESeq(dds, quiet = TRUE)
  
  # LFC shrinkage (ashr handles arbitrary contrasts; fallback to unshrunken)
  res_obj <- tryCatch(
    lfcShrink(dds, contrast = contrast, type = shrink_type, quiet = TRUE),
    error = function(e) {
      message("    -> LFC shrinkage failed (", e$message, "); unshrunken results")
      results(dds, contrast = contrast)
    }
  )
  
  res_df <- as.data.frame(res_obj) %>%
    dplyr::rename(p_val_adj  = padj,
                  avg_log2FC = log2FoldChange,
                  p_val      = pvalue) %>%
    dplyr::arrange(p_val_adj)
  
  write.csv(res_df, out_csv)   # rownames = gene symbols (OPIS-compatible)
  invisible(res_df)
}

# ============================ Per-lineage DE (pseudobulk, primary) ===========
message("\n=== Pseudobulk DE: per LINEAGE ===")
for (lin in levels(obj$tentative_lineage)) {
  message("\n[Lineage: ", lin, "]")
  sub <- subset(obj, subset = tentative_lineage == lin)
  pb  <- make_pb(sub, min_cells_per_sample = 10)
  if (is.null(pb)) { message("  insufficient pseudobulks; skipping"); next }
  
  dir.create(file.path(dir_PB_lin, lin), showWarnings = FALSE, recursive = TRUE)
  write.csv(pb$col_meta,
            file.path(dir_PB_lin, lin, "pseudobulk_samples.csv"),
            row.names = FALSE)
  
  for (cn in names(contrasts_list)) {
    message("  -> ", cn)
    out_csv <- file.path(dir_PB_lin, lin,
                         paste0(safe_name(lin), "_", cn, ".csv"))
    run_de(pb$counts, pb$col_meta, contrasts_list[[cn]], out_csv)
  }
}

# ============================ Per-cluster DE (pseudobulk, secondary) =========
message("\n=== Pseudobulk DE: per CLUSTER ===")
for (cl in levels(obj$clusters)) {
  message("\n[Cluster ", cl, "]")
  sub <- subset(obj, subset = clusters == cl)
  pb  <- make_pb(sub, min_cells_per_sample = 10)
  if (is.null(pb)) { message("  insufficient pseudobulks; skipping"); next }
  
  stem <- cluster_stem(cl)
  cl_dir <- file.path(dir_PB_cl, stem)
  dir.create(cl_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(pb$col_meta,
            file.path(cl_dir, "pseudobulk_samples.csv"),
            row.names = FALSE)
  
  for (cn in names(contrasts_list)) {
    message("  -> ", cn)
    out_csv <- file.path(cl_dir,
                         paste0(stem, "_", cn, ".csv"))
    run_de(pb$counts, pb$col_meta, contrasts_list[[cn]], out_csv)
  }
}

# ============================ MAST helper ====================================
# Single-cell DE via FindMarkers + MAST, with nCount_RNA as latent variable
# (the OPIS pattern). Run for comparability and as a broader exploratory
# screen -- the inferentially-honest hit list comes from pseudobulk above.
run_mast <- function(seu_subset, ident1, ident2, out_csv,
                     group_col = "condition",
                     latent    = c("nCount_RNA"),
                     min_cells = 10) {
  grp <- seu_subset@meta.data[[group_col]]
  tab <- table(grp)
  if (!all(c(ident1, ident2) %in% names(tab)) ||
      tab[ident1] < min_cells || tab[ident2] < min_cells) {
    message("    -> skip: <", min_cells, " cells in one or both groups")
    return(invisible(NULL))
  }
  de <- tryCatch(
    FindMarkers(seu_subset,
                ident.1     = ident1,
                ident.2     = ident2,
                group.by    = group_col,
                test.use    = "MAST",
                latent.vars = latent,
                verbose     = FALSE),
    error = function(e) { message("    -> MAST failed: ", e$message); NULL }
  )
  if (is.null(de) || nrow(de) == 0) return(invisible(NULL))
  de <- de[order(de$p_val_adj), ]
  write.csv(de, out_csv)   # rownames = gene symbols (OPIS-compatible)
  invisible(de)
}

# ============================ Per-lineage DE (MAST) ==========================
message("\n=== MAST DE: per LINEAGE ===")
for (lin in levels(obj$tentative_lineage)) {
  message("\n[Lineage: ", lin, "]")
  sub <- subset(obj, subset = tentative_lineage == lin)
  dir.create(file.path(dir_MAST_lin, lin), showWarnings = FALSE, recursive = TRUE)
  
  for (cn in names(contrasts_list)) {
    contrast <- contrasts_list[[cn]]
    message("  -> ", cn)
    out_csv <- file.path(dir_MAST_lin, lin,
                         paste0(safe_name(lin), "_", cn, ".csv"))
    run_mast(sub, ident1 = contrast[2], ident2 = contrast[3], out_csv = out_csv)
  }
}

# ============================ Per-cluster DE (MAST) ==========================
message("\n=== MAST DE: per CLUSTER ===")
for (cl in levels(obj$clusters)) {
  message("\n[Cluster ", cl, "]")
  sub <- subset(obj, subset = clusters == cl)
  stem <- cluster_stem(cl)
  cl_dir <- file.path(dir_MAST_cl, stem)
  dir.create(cl_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (cn in names(contrasts_list)) {
    contrast <- contrasts_list[[cn]]
    message("  -> ", cn)
    out_csv <- file.path(cl_dir,
                         paste0(stem, "_", cn, ".csv"))
    run_mast(sub, ident1 = contrast[2], ident2 = contrast[3], out_csv = out_csv)
  }
}

# ============================ Summary: pseudobulk vs MAST ====================
count_sig <- function(csv, thr = 0.05) {
  df <- tryCatch(read.csv(csv, row.names = 1), error = function(e) NULL)
  if (is.null(df) || !"p_val_adj" %in% colnames(df)) return(NA_integer_)
  sum(df$p_val_adj < thr, na.rm = TRUE)
}

contrast_pat <- "_(miR29a_vs_Scr|miR29a_vs_EV|EV_vs_Scr)\\.csv$"

# walk a method's directory tree and collect DE counts per compartment/contrast
gather_counts <- function(method_root, method_label) {
  rows <- list()
  lin_root <- file.path(method_root, "by_lineage")
  cl_root  <- file.path(method_root, "by_cluster")
  for (lin_dir in list.dirs(lin_root, recursive = FALSE)) {
    for (csv in list.files(lin_dir, pattern = contrast_pat, full.names = TRUE)) {
      cn <- sub(paste0(".*", contrast_pat), "\\1", csv)
      rows[[length(rows)+1]] <- data.frame(
        group    = paste0("Lineage: ", basename(lin_dir)),
        contrast = cn, n_DE = count_sig(csv), method = method_label
      )
    }
  }
  for (cl_dir in list.dirs(cl_root, recursive = FALSE)) {
    for (csv in list.files(cl_dir, pattern = contrast_pat, full.names = TRUE)) {
      cn <- sub(paste0(".*", contrast_pat), "\\1", csv)
      rows[[length(rows)+1]] <- data.frame(
        group    = gsub("_", " ", basename(cl_dir)),
        contrast = cn, n_DE = count_sig(csv), method = method_label
      )
    }
  }
  if (length(rows) == 0) NULL else do.call(rbind, rows)
}

summary_df <- rbind(
  gather_counts(file.path(de_base, "DGE_pseudobulk"), "Pseudobulk"),
  gather_counts(file.path(de_base, "DGE_MAST"),       "MAST")
)

if (!is.null(summary_df) && nrow(summary_df) > 0) {
  summary_df$contrast <- factor(summary_df$contrast,
                                levels = c("miR29a_vs_Scr","miR29a_vs_EV","EV_vs_Scr"))
  summary_df$method   <- factor(summary_df$method, levels = c("Pseudobulk","MAST"))
  write.csv(summary_df, file.path(de_base, "DGE_summary_counts.csv"),
            row.names = FALSE)
  
  # side-by-side bars, faceted by method (MAST counts will dwarf pseudobulk
  # -- expected; pseudobulk is calibrated and MAST is significance-inflated
  # at n=2 reps. Compare TRENDS across compartments, not absolute counts.)
  p_sum <- ggplot(summary_df, aes(x = group, y = n_DE, fill = contrast)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ method, ncol = 1, scales = "free_y") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "DE gene counts: pseudobulk vs MAST (padj < 0.05)",
         x = NULL, y = "# DE genes")
  ggsave(file.path(de_base, "DGE_summary_counts.png"),
         p_sum, width = 14, height = 10, dpi = 300, bg = "white")
}

# ============================ Method-overlap table ===========================
# For each (compartment x contrast), how many sig genes overlap between
# pseudobulk and MAST? This is the headline diagnostic: high overlap on the
# strong hits = methods agree on what's real; MAST-only = likely inflated;
# PB-only = rarer (suggests calibration disagreement on a real signal).
sig_genes_from <- function(csv, thr = 0.05) {
  df <- tryCatch(read.csv(csv, row.names = 1), error = function(e) NULL)
  if (is.null(df) || !"p_val_adj" %in% colnames(df)) return(character(0))
  rownames(df)[!is.na(df$p_val_adj) & df$p_val_adj < thr]
}

build_overlap <- function() {
  rows <- list()
  for (lin in levels(obj$tentative_lineage)) {
    for (cn in names(contrasts_list)) {
      pb_csv   <- file.path(dir_PB_lin,   lin, paste0(safe_name(lin), "_", cn, ".csv"))
      mast_csv <- file.path(dir_MAST_lin, lin, paste0(safe_name(lin), "_", cn, ".csv"))
      pb   <- sig_genes_from(pb_csv); ma <- sig_genes_from(mast_csv)
      rows[[length(rows)+1]] <- data.frame(
        group = paste0("Lineage: ", lin), contrast = cn,
        n_pseudobulk = length(pb), n_MAST = length(ma),
        n_overlap    = length(intersect(pb, ma)),
        n_PB_only    = length(setdiff(pb, ma)),
        n_MAST_only  = length(setdiff(ma, pb))
      )
    }
  }
  for (cl in levels(obj$clusters)) {
    for (cn in names(contrasts_list)) {
      stem <- cluster_stem(cl)
      pb_csv   <- file.path(dir_PB_cl,   stem, paste0(stem, "_", cn, ".csv"))
      mast_csv <- file.path(dir_MAST_cl, stem, paste0(stem, "_", cn, ".csv"))
      pb   <- sig_genes_from(pb_csv); ma <- sig_genes_from(mast_csv)
      rows[[length(rows)+1]] <- data.frame(
        group = cluster_pretty(cl), contrast = cn,
        n_pseudobulk = length(pb), n_MAST = length(ma),
        n_overlap    = length(intersect(pb, ma)),
        n_PB_only    = length(setdiff(pb, ma)),
        n_MAST_only  = length(setdiff(ma, pb))
      )
    }
  }
  do.call(rbind, rows)
}

overlap_df <- build_overlap()
write.csv(overlap_df, file.path(de_base, "DGE_method_overlap.csv"), row.names = FALSE)

message("\nDone. See ", de_base)
###############################################################################