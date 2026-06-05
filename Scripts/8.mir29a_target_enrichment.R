###############################################################################
# 08 -- miR-29a Target Enrichment in DE Results                               #
# Mouse pre-infusion CAR-T scRNA-seq, miR-29a project                         #
###############################################################################
# Two complementary tests against the conserved miR-29-3p TargetScan list
# (Resources/miR29a_targetscan_conserved.csv), per DE contrast/compartment/
# method:
#
#   (1) Spearman correlation -- TargetScan cumulative weighted score vs
#       experimental avg_log2FC. The HEADLINE figure: more-negative
#       TargetScan score (stronger predicted target) should correspond to
#       more-negative experimental log2FC (more strongly repressed) when
#       miR-29a is functional. Each contrast yields a scatter with regression
#       line and Spearman rho/p-value.
#
#   (2) Hypergeometric test -- are significantly DOWN-regulated genes
#       enriched for miR-29a targets relative to the tested gene universe?
#       Run against the top-N target list (default 200) and against the full
#       conserved list, with the universe defined per contrast (rows in
#       that DE CSV).
#
# Runs on BOTH script 05 outputs:
#   DGE_pseudobulk  (primary; calibrated)
#   DGE_MAST        (secondary; broader; significance inflated at n=2)
#
# Writes per-contrast scatter plots, per-contrast summary CSVs, and a
# top-level summary table aggregating Spearman rho / hypergeometric p-values
# across every compartment x contrast x method.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
})

# ============================ Paths ==========================================
project_dir   <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a"
resources_dir <- file.path(project_dir, "Resources")
mir29a_csv    <- file.path(resources_dir, "miR29a_targetscan_conserved.csv")

de_root       <- file.path(project_dir, "Differential_Expression")
out_base      <- file.path(project_dir, "miR29a_Target_Enrichment")
plot_dir      <- file.path(out_base, "Plots")
table_dir     <- file.path(out_base, "Tables")
for (d in c(plot_dir, table_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ============================ Load target list ===============================
targets <- read.csv(mir29a_csv, stringsAsFactors = FALSE)
stopifnot(all(c("gene", "cum_weighted_context_score") %in% colnames(targets)))
targets <- targets %>% dplyr::arrange(cum_weighted_context_score)   # most neg first

top_n              <- 200            # top-N target set for hypergeometric
target_top_set     <- head(targets$gene, top_n)
target_full_set    <- targets$gene
message("Target list: ", nrow(targets), " total; top-", top_n,
        " (score <= ", round(targets$cum_weighted_context_score[top_n], 3), ")")

# ============================ Discover DE CSVs ===============================
# Scan both methods' lineage and cluster trees for *_<contrast>.csv files.
contrast_pat <- "_(miR29a_vs_Scr|miR29a_vs_EV|EV_vs_Scr)\\.csv$"
methods      <- c("DGE_pseudobulk", "DGE_MAST")

de_index <- list()
for (mth in methods) {
  for (level in c("by_lineage", "by_cluster")) {
    root <- file.path(de_root, mth, level)
    if (!dir.exists(root)) next
    for (sub in list.dirs(root, recursive = FALSE)) {
      for (csv in list.files(sub, pattern = contrast_pat, full.names = TRUE)) {
        cn <- sub(paste0(".*", contrast_pat), "\\1", csv)
        de_index[[length(de_index)+1]] <- data.frame(
          method    = mth,
          level     = level,                  # by_lineage / by_cluster
          group     = basename(sub),          # e.g. "CD8" or "02_Proliferative_CD8_effector"
          contrast  = cn,
          csv       = csv,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
de_index <- dplyr::bind_rows(de_index)
if (nrow(de_index) == 0) stop("No DE CSVs found under ", de_root)
message("Found ", nrow(de_index), " DE CSVs")

# ============================ Helpers ========================================
safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# Hypergeometric enrichment of sig-down genes for target set
hyper_enrich <- function(sig_down, universe, target_set) {
  K <- sum(universe %in% target_set)        # successes in pop
  N <- length(universe)                     # pop size
  n <- length(sig_down)                     # sample size
  k <- sum(sig_down %in% target_set)        # successes in sample
  if (n == 0 || K == 0)
    return(list(k = k, n = n, K = K, N = N,
                p_value = NA_real_, odds = NA_real_, expected = NA_real_))
  p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)     # P(>= k)
  expected <- n * K / N
  odds <- (k / max(n - k, 1)) / (K / max(N - K, 1))
  list(k = k, n = n, K = K, N = N, p_value = p, odds = odds, expected = expected)
}

# scatter: TargetScan score vs experimental LFC for one DE result
scatter_score_vs_lfc <- function(merged, title, out_png,
                                 label_top = 15) {
  if (nrow(merged) < 5) return(invisible(NULL))
  rho <- suppressWarnings(cor.test(merged$cum_weighted_context_score,
                                   merged$avg_log2FC, method = "spearman",
                                   exact = FALSE))
  # label the most strongly repressed targets (most neg LFC) plus the strongest
  # predicted (most neg target score)
  to_label <- unique(c(
    merged %>% arrange(avg_log2FC)                  %>% head(label_top) %>% pull(gene),
    merged %>% arrange(cum_weighted_context_score)  %>% head(label_top) %>% pull(gene)
  ))
  p <- ggplot(merged,
              aes(x = cum_weighted_context_score, y = avg_log2FC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(alpha = 0.45, size = 1.3, color = "steelblue4") +
    geom_smooth(method = "lm", color = "firebrick", se = TRUE, linewidth = 0.6) +
    ggrepel::geom_text_repel(
      data = subset(merged, gene %in% to_label),
      aes(label = gene), size = 3, max.overlaps = 25, min.segment.length = 0
    ) +
    labs(
      title = title,
      subtitle = sprintf("n = %d  |  Spearman \u03c1 = %.3f  |  p = %.2e",
                         nrow(merged), rho$estimate, rho$p.value),
      x = "TargetScan cumulative weighted context++ score (more neg = stronger predicted target)",
      y = "Experimental avg log2FC (more neg = more repressed)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(out_png, p, width = 10, height = 8, dpi = 300, bg = "white")
  list(rho = unname(rho$estimate), p = rho$p.value, n = nrow(merged))
}

# ============================ Main loop ======================================
summary_rows <- list()

for (i in seq_len(nrow(de_index))) {
  rec <- de_index[i, ]
  de  <- tryCatch(read.csv(rec$csv, row.names = 1, check.names = FALSE),
                  error = function(e) NULL)
  if (is.null(de) || !"avg_log2FC" %in% colnames(de) ||
      !"p_val_adj" %in% colnames(de)) next
  
  # universe = tested genes (rows of DE result)
  de$gene <- rownames(de)
  universe <- de$gene
  
  # significantly DOWN-regulated genes
  sig_down <- de %>%
    dplyr::filter(!is.na(p_val_adj), p_val_adj < 0.05, avg_log2FC < 0) %>%
    dplyr::pull(gene)
  sig_up <- de %>%
    dplyr::filter(!is.na(p_val_adj), p_val_adj < 0.05, avg_log2FC > 0) %>%
    dplyr::pull(gene)
  
  # hypergeometric -- target top-N
  ht_top  <- hyper_enrich(sig_down, universe, target_top_set)
  ht_full <- hyper_enrich(sig_down, universe, target_full_set)
  # hypergeometric for UP genes too (negative control -- shouldn't be enriched)
  ht_top_up  <- hyper_enrich(sig_up, universe, target_top_set)
  
  # spearman: merge DE with target weights
  merged <- de %>%
    dplyr::inner_join(targets, by = "gene") %>%
    dplyr::filter(!is.na(avg_log2FC))
  
  # scatter plot
  group_safe <- safe(rec$group)
  out_png <- file.path(plot_dir,
                       paste0(rec$method, "_", rec$level, "_", group_safe,
                              "_", rec$contrast, "_scatter.png"))
  title <- paste0(rec$method, " | ", gsub("_", " ", rec$group),
                  " | ", gsub("_", " ", rec$contrast))
  sp <- scatter_score_vs_lfc(merged, title, out_png)
  
  summary_rows[[length(summary_rows)+1]] <- data.frame(
    method                 = rec$method,
    level                  = rec$level,
    group                  = rec$group,
    contrast               = rec$contrast,
    n_tested_universe      = length(universe),
    n_in_target_top        = ht_top$K,
    n_in_target_full       = ht_full$K,
    n_sig_down             = length(sig_down),
    n_sig_up               = length(sig_up),
    
    spearman_rho           = if (!is.null(sp)) sp$rho else NA_real_,
    spearman_p             = if (!is.null(sp)) sp$p   else NA_real_,
    spearman_n             = if (!is.null(sp)) sp$n   else NA_integer_,
    
    hyper_top_overlap      = ht_top$k,
    hyper_top_expected     = ht_top$expected,
    hyper_top_odds         = ht_top$odds,
    hyper_top_p            = ht_top$p_value,
    hyper_full_overlap     = ht_full$k,
    hyper_full_p           = ht_full$p_value,
    
    hyper_top_overlap_UP   = ht_top_up$k,
    hyper_top_p_UP         = ht_top_up$p_value,
    
    stringsAsFactors = FALSE
  )
}

summary_df <- dplyr::bind_rows(summary_rows)
# BH-adjust enrichment p-values within each (method, level) family
summary_df <- summary_df %>%
  dplyr::group_by(method, level) %>%
  dplyr::mutate(hyper_top_p_BH  = p.adjust(hyper_top_p,  method = "BH"),
                hyper_full_p_BH = p.adjust(hyper_full_p, method = "BH"),
                spearman_p_BH   = p.adjust(spearman_p,   method = "BH")) %>%
  dplyr::ungroup()

write.csv(summary_df, file.path(table_dir, "miR29a_enrichment_summary.csv"),
          row.names = FALSE)

# ============================ Headline summary plot ==========================
# Bar plot: Spearman rho per (compartment, contrast) for pseudobulk lineage-level
# results -- the manuscript headline.
if (nrow(summary_df) > 0) {
  headline <- summary_df %>%
    dplyr::filter(method == "DGE_pseudobulk", level == "by_lineage") %>%
    dplyr::mutate(contrast = factor(contrast,
                                    levels = c("miR29a_vs_Scr",
                                               "miR29a_vs_EV",
                                               "EV_vs_Scr")))
  if (nrow(headline) > 0) {
    p_head <- ggplot(headline,
                     aes(x = group, y = spearman_rho, fill = contrast)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_hline(yintercept = 0, color = "grey50") +
      labs(title = "miR-29a target enrichment in DE (pseudobulk, lineage)",
           subtitle = "Spearman correlation: TargetScan score vs experimental log2FC",
           x = NULL,
           y = "Spearman \u03c1 (positive = targets repressed in ident.1)") +
      theme_minimal(base_size = 12)
    ggsave(file.path(plot_dir, "headline_spearman_lineage_pseudobulk.png"),
           p_head, width = 10, height = 6, dpi = 300, bg = "white")
  }
}

message("\nDone. See ", out_base)
message("Headline file: Tables/miR29a_enrichment_summary.csv")
message("  Look at: spearman_rho > 0 with spearman_p_BH < 0.05")
message("  AND     hyper_top_p_BH < 0.05 with overlap > expected")
###############################################################################