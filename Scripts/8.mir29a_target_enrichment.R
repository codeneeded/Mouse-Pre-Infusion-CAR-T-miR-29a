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
                                 label_top = 15, sig_thresh = 0.05) {
  if (nrow(merged) < 5) return(invisible(NULL))
  rho <- suppressWarnings(cor.test(merged$cum_weighted_context_score,
                                   merged$avg_log2FC, method = "spearman",
                                   exact = FALSE))
  
  # classify each gene by significance + direction
  merged$sig_class <- dplyr::case_when(
    is.na(merged$p_val_adj)                                  ~ "ns",
    merged$p_val_adj < sig_thresh & merged$avg_log2FC < 0    ~ "sig_down",
    merged$p_val_adj < sig_thresh & merged$avg_log2FC > 0    ~ "sig_up",
    TRUE                                                      ~ "ns"
  )
  merged$sig_class <- factor(merged$sig_class,
                             levels = c("ns", "sig_up", "sig_down"))
  
  # label only the strongest validated repression: sig-down genes ranked by
  # combined evidence (most negative target score + most negative LFC).
  sig_down <- merged %>% dplyr::filter(sig_class == "sig_down")
  to_label <- if (nrow(sig_down) > 0) {
    sig_down %>%
      dplyr::mutate(combined = cum_weighted_context_score + avg_log2FC) %>%
      dplyr::arrange(combined) %>%
      head(label_top) %>%
      dplyr::pull(gene)
  } else character(0)
  
  cols <- c("ns" = "grey75", "sig_up" = "firebrick3", "sig_down" = "steelblue3")
  lbls <- c("ns" = "n.s.", "sig_up" = "Sig. up", "sig_down" = "Sig. down")
  
  p <- ggplot(merged,
              aes(x = cum_weighted_context_score, y = avg_log2FC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(aes(color = sig_class), alpha = 0.6, size = 1.3) +
    # regression on ALL points (one line, not per sig class)
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 0.6) +
    scale_color_manual(values = cols, labels = lbls, drop = FALSE) +
    ggrepel::geom_text_repel(
      data = subset(merged, gene %in% to_label),
      aes(label = gene), size = 3, max.overlaps = 25,
      min.segment.length = 0, color = "black"
    ) +
    labs(
      title = title,
      subtitle = sprintf("n = %d  |  Spearman \u03c1 = %.3f  |  p = %.2e",
                         nrow(merged), rho$estimate, rho$p.value),
      x = "TargetScan cumulative weighted context++ score (more neg = stronger predicted target)",
      y = "Experimental avg log2FC (more neg = more repressed)",
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title       = element_text(face = "bold"),
          legend.position  = "top")
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
  
  # scatter plot -- organized by method x level to keep the per-contrast folder
  # navigable rather than flat.
  group_safe <- safe(rec$group)
  scatter_subdir <- file.path(plot_dir, "Per_contrast",
                              sub("DGE_", "", rec$method), rec$level)
  dir.create(scatter_subdir, showWarnings = FALSE, recursive = TRUE)
  out_png <- file.path(scatter_subdir,
                       paste0(group_safe, "_", rec$contrast, "_scatter.png"))
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

# ============================ Summary plots =================================
# Four manuscript-grade summary plots per DE method:
#   1. Spearman ρ forest plot (the statistical claim)
#   2. Hypergeometric enrichment dot plot (top-200 robustness check)
#   3. Direction asymmetry scatter (miRNA-repressor fingerprint)
#   4. Canonical miR-29a target heatmap (the biological claim)

# Pretty compartment labels and a stable ordering: lineages first, then clusters
# in numeric order. Consistent across all summary plots and both methods.
make_compartment_label <- function(group, level) {
  ifelse(level == "by_lineage",
         paste0(group, "  (lineage)"),
         gsub("_+", " ", group))
}
build_compartment_order <- function(df) {
  lin <- df %>% dplyr::filter(level == "by_lineage") %>% dplyr::distinct(group)
  lin_order <- c("CD4", "CD8", "Non-T", "Innate-like")
  lin <- lin$group[order(match(lin$group, lin_order))]
  lin_labs <- paste0(lin, "  (lineage)")
  
  cl <- df %>% dplyr::filter(level == "by_cluster") %>% dplyr::distinct(group)
  cl_num <- suppressWarnings(as.numeric(sub("^(\\d+).*", "\\1", cl$group)))
  cl <- cl$group[order(cl_num)]
  cl_labs <- gsub("_+", " ", cl)
  rev(c(lin_labs, cl_labs))  # reverse so first row plots at TOP
}

contrast_levels <- c("miR29a_vs_Scr", "miR29a_vs_EV", "EV_vs_Scr")
contrast_cols   <- c("miR29a_vs_Scr" = "#1F3A6B",   # deep navy -- PRIMARY
                     "miR29a_vs_EV"  = "#7BAFD4",   # mid-blue  -- secondary
                     "EV_vs_Scr"     = "grey65")    # grey      -- control
contrast_labs   <- c("miR29a_vs_Scr" = "miR29a vs Scr  (primary)",
                     "miR29a_vs_EV"  = "miR29a vs EV  (secondary)",
                     "EV_vs_Scr"     = "EV vs Scr  (control)")

# helper to read one DE CSV for the heatmap panel
read_de_panel <- function(method, group, level, contrast = "miR29a_vs_Scr") {
  fname <- if (level == "by_lineage")
    paste0(safe(group), "_", contrast, ".csv")
  else
    paste0(group, "_", contrast, ".csv")
  csv <- file.path(de_root, method, level, group, fname)
  if (!file.exists(csv)) return(NULL)
  tryCatch(read.csv(csv, row.names = 1, check.names = FALSE),
           error = function(e) NULL)
}

# canonical miR-29 target panel -- T cell relevant.
# 'Bach2_Foxo3' = quiescence/memory TFs; 'Tbx21_Eomes' = effector T-box TFs;
# 'Tet2_Tet3' / 'Dnmt3a_Dnmt3b' / 'Tdg' = DNA (de)methylation; 'Atad2b_Mycn'
# etc. are strong data-driven hits already labelled on the per-contrast scatters.
canonical_genes <- c(
  "Tet2", "Tet3", "Dnmt3a", "Dnmt3b", "Tdg",            # methylation axis
  "Tbx21", "Eomes",                                      # T-box effector TFs
  "Bach2", "Foxo3",                                      # memory / quiescence
  "Icos",                                                # costimulation
  "Atad2b", "Mycn", "Pdgfa", "Efna5", "Gas7", "Dab2ip"   # data-driven hits
)

# generate the four plots for one method ----------------------------------
make_summary_plots <- function(mth, sub_df, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  comp_order <- build_compartment_order(sub_df)
  d <- sub_df %>%
    dplyr::mutate(
      compartment = make_compartment_label(group, level),
      compartment = factor(compartment, levels = comp_order),
      contrast    = factor(contrast, levels = contrast_levels)
    )
  
  base_theme <- theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(color = "grey35", size = 10),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "white", color = NA),
          plot.margin   = margin(t = 12, r = 14, b = 12, l = 12))
  
  # ----- Plot 1: Spearman forest --------------------------------------------
  forest_d <- d %>%
    dplyr::filter(!is.na(spearman_rho), !is.na(spearman_n), spearman_n > 3) %>%
    dplyr::mutate(
      z_rho = atanh(pmin(pmax(spearman_rho, -0.999), 0.999)),
      z_se  = 1 / sqrt(spearman_n - 3),
      ci_lo = tanh(z_rho - 1.96 * z_se),
      ci_hi = tanh(z_rho + 1.96 * z_se),
      sig_score = -log10(pmax(spearman_p_BH, 1e-15))
    )
  
  if (nrow(forest_d) > 0) {
    p1 <- ggplot(forest_d, aes(x = spearman_rho, y = compartment, color = contrast)) +
      annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
               alpha = 0.025, fill = "red") +
      annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
               alpha = 0.025, fill = "blue") +
      geom_vline(xintercept = 0, color = "grey35", linewidth = 0.4) +
      geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                     position = position_dodge(width = 0.7),
                     height = 0, linewidth = 0.55, alpha = 0.75) +
      geom_point(aes(size = sig_score),
                 position = position_dodge(width = 0.7)) +
      scale_color_manual(values = contrast_cols, labels = contrast_labs,
                         drop = FALSE) +
      scale_size_continuous(range = c(1.6, 5.5), name = expression(-log[10]~p[BH]),
                            breaks = c(2, 5, 10)) +
      labs(title    = sprintf("Spearman \u03c1: TargetScan score vs experimental log2FC  (%s)",
                              gsub("DGE_", "", mth)),
           subtitle = "Each compartment shows three contrasts with 95% CI. EV vs Scr should sit at \u03c1 = 0.",
           x = expression(Spearman~rho),
           y = NULL, color = NULL) +
      base_theme +
      theme(panel.grid.major.y = element_blank()) +
      guides(color = guide_legend(order = 1, override.aes = list(size = 3.5)),
             size  = guide_legend(order = 2))
    
    ggsave(file.path(out_dir, "1_spearman_forest.png"),
           p1, width = 11, height = 9.5, dpi = 300, bg = "white")
  }
  
  # ----- Plot 2: Hypergeometric enrichment dot plot -------------------------
  dot_d <- d %>%
    dplyr::filter(!is.na(hyper_top_p_BH), n_sig_down > 0) %>%
    dplyr::mutate(
      sig_score = -log10(pmax(hyper_top_p_BH, 1e-15)),
      log_odds  = log2(pmin(hyper_top_odds + 1, 64))   # cap for visual range
    )
  
  if (nrow(dot_d) > 0) {
    p2 <- ggplot(dot_d, aes(x = contrast, y = compartment)) +
      geom_point(aes(size = sig_score, fill = log_odds),
                 shape = 21, color = "grey25", stroke = 0.4) +
      scale_size_continuous(range = c(2, 12),
                            name  = expression(-log[10]~p[BH])) +
      scale_fill_distiller(palette = "YlOrRd", direction = 1,
                           name = expression(log[2](odds+1))) +
      scale_x_discrete(labels = c(miR29a_vs_Scr = "miR29a vs Scr",
                                  miR29a_vs_EV  = "miR29a vs EV",
                                  EV_vs_Scr     = "EV vs Scr")) +
      labs(title    = sprintf("Hypergeometric enrichment of sig-down genes in top-200 TargetScan targets  (%s)",
                              gsub("DGE_", "", mth)),
           subtitle = "Larger and warmer = stronger enrichment. Empty rows: too few sig-down genes to test.",
           x = NULL, y = NULL) +
      base_theme +
      theme(axis.text.x      = element_text(angle = 25, hjust = 1),
            panel.grid.major = element_line(color = "grey94", linewidth = 0.3))
    
    ggsave(file.path(out_dir, "2_hyper_enrichment_dotplot.png"),
           p2, width = 9, height = 9.5, dpi = 300, bg = "white")
  }
  
  # ----- Plot 3: Direction asymmetry ----------------------------------------
  asym_d <- d %>%
    dplyr::filter(n_sig_down + n_sig_up > 0) %>%
    dplyr::mutate(level_lab = ifelse(level == "by_lineage", "Lineage", "Cluster"))
  
  if (nrow(asym_d) > 0) {
    max_n <- max(asym_d$n_sig_down, asym_d$n_sig_up, 2)
    p3 <- ggplot(asym_d, aes(x = n_sig_up + 1, y = n_sig_down + 1)) +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "grey45") +
      geom_point(aes(color = contrast, shape = level_lab),
                 size = 3.2, alpha = 0.85, stroke = 0.4) +
      scale_x_log10() + scale_y_log10() +
      scale_color_manual(values = contrast_cols, labels = contrast_labs) +
      scale_shape_manual(values = c(Lineage = 16, Cluster = 17), name = NULL) +
      annotate("text", x = max_n * 0.6, y = 2, hjust = 0.5,
               label = "more sig-UP than sig-DOWN",
               color = "grey55", size = 3.3) +
      annotate("text", x = 2, y = max_n * 0.6, hjust = 0,
               label = "more sig-DOWN \u2192 miRNA repressor fingerprint",
               color = "grey20", size = 3.3, fontface = "bold") +
      labs(title    = sprintf("Direction asymmetry: sig-DOWN vs sig-UP per compartment  (%s)",
                              gsub("DGE_", "", mth)),
           subtitle = "Below the diagonal = repression dominates. Control (EV vs Scr) should sit on the diagonal.",
           x = "n sig-UP  (+1, log scale)",
           y = "n sig-DOWN  (+1, log scale)",
           color = NULL) +
      base_theme +
      guides(color = guide_legend(order = 1, override.aes = list(size = 3.5)),
             shape = guide_legend(order = 2, override.aes = list(size = 3.5)))
    
    ggsave(file.path(out_dir, "3_direction_asymmetry.png"),
           p3, width = 9, height = 7.5, dpi = 300, bg = "white")
  }
  
  # ----- Plot 4: Canonical target heatmap -----------------------------------
  comp_list <- d %>% dplyr::distinct(group, level)
  hm_d <- do.call(rbind, lapply(seq_len(nrow(comp_list)), function(i) {
    de <- read_de_panel(mth, comp_list$group[i], comp_list$level[i],
                        contrast = "miR29a_vs_Scr")
    if (is.null(de) || !"avg_log2FC" %in% colnames(de)) return(NULL)
    hits <- de[rownames(de) %in% canonical_genes,
               c("p_val_adj", "avg_log2FC"), drop = FALSE]
    if (nrow(hits) == 0) return(NULL)
    hits$gene  <- rownames(hits)
    hits$group <- comp_list$group[i]
    hits$level <- comp_list$level[i]
    hits
  }))
  
  if (!is.null(hm_d) && nrow(hm_d) > 0) {
    hm_d <- hm_d %>%
      dplyr::mutate(
        compartment = make_compartment_label(group, level),
        compartment = factor(compartment, levels = rev(comp_order)),
        gene        = factor(gene, levels = canonical_genes),
        lfc_clip    = pmin(pmax(avg_log2FC, -1.5), 1.5),
        sig_lab     = dplyr::case_when(
          is.na(p_val_adj)   ~ "",
          p_val_adj < 0.001  ~ "***",
          p_val_adj < 0.01   ~ "**",
          p_val_adj < 0.05   ~ "*",
          TRUE               ~ ""
        )
      )
    
    p4 <- ggplot(hm_d, aes(x = compartment, y = gene, fill = lfc_clip)) +
      geom_tile(color = "white", linewidth = 0.6) +
      geom_text(aes(label = sig_lab), size = 3.5,
                color = "black", fontface = "bold") +
      scale_fill_distiller(palette = "RdBu", limits = c(-1.5, 1.5),
                           name = expression(avg~log[2]~FC),
                           na.value = "grey90") +
      scale_y_discrete(limits = rev) +
      labs(title    = sprintf("Canonical miR-29a targets: log2FC (miR29a vs Scr)  (%s)",
                              gsub("DGE_", "", mth)),
           subtitle = "Blue = repressed in miR-29a. * p_adj<0.05  ** p_adj<0.01  *** p_adj<0.001  (LFC clipped at \u00b11.5)",
           x = NULL, y = NULL) +
      base_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(face = "italic"),
            panel.grid  = element_blank(),
            plot.margin = margin(t = 12, r = 14, b = 60, l = 12))
    
    ggsave(file.path(out_dir, "4_canonical_targets_heatmap.png"),
           p4, width = 13, height = 8.5, dpi = 300, bg = "white")
  }
}

# Run for each method --------------------------------------------------------
if (nrow(summary_df) > 0) {
  for (mth in unique(summary_df$method)) {
    sub <- summary_df %>% dplyr::filter(method == mth)
    if (nrow(sub) == 0) next
    out_dir <- file.path(plot_dir,
                         paste0("Summary_", sub("DGE_", "", mth)))
    message("Summary plots for ", mth, " -> ", out_dir)
    make_summary_plots(mth, sub, out_dir)
  }
}

message("\nDone. See ", out_base)
message("Per-contrast scatters:  Plots/Per_contrast/<method>/<level>/")
message("Summary plots:          Plots/Summary_<method>/  (4 plots each)")
message("Headline table:         Tables/miR29a_enrichment_summary.csv")
message("  Read out: spearman_rho > 0 with spearman_p_BH < 0.05")
message("       AND  hyper_top_p_BH < 0.05 with overlap > expected")
###############################################################################