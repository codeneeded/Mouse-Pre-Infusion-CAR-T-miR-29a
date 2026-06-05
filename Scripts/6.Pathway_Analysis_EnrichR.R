###############################################################################
# 06 -- Pathway & miRNA-target enrichment (EnrichR)                           #
# Mouse pre-infusion CAR-T scRNA-seq, miR-29a project                         #
###############################################################################
# Reads DGE CSVs produced by script 05 and runs enrichment in three groups:
#   1. TF / regulators           (TRRUST, ChEA, JASPAR PWMs)
#   2. Pathways / ontologies     (KEGG/Wiki MOUSE; GO BP, Reactome,
#                                 MSigDB Hallmark, BioPlanet, Panther HUMAN)
#   3. miRNA target interactions (miRTarBase, TargetScan)
#
# miRNA logic: miR-29a represses its targets, so when miR-29a is OVER-
# expressed (miR29a vs control), its validated/predicted targets should appear
# DOWN-regulated. We therefore enrich miRNA databases on the DOWN-regulated
# gene list SEPARATELY from UP, and flag any "mir-29" hits explicitly.
#
# Gene symbols: DGE input is mouse Title case (e.g. Tbx21). Most EnrichR
# libraries are human-centric, so we pass toupper() symbols to TF /
# pathway-human / miRNA libraries (mouse->human ortholog approximation by
# symbol, which works for most T-cell genes). Mouse-case symbols are passed
# to the mouse-specific KEGG / Wiki libraries. For paper-grade orthology,
# swap toupper() for babelgene::orthologs() or biomaRt.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(enrichR); library(openxlsx); library(ggplot2)
})

# ============================ Paths ==========================================
project_dir <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a"

# Which DE method's results to feed into pathway analysis?
# Options (must match script 05's output dirs):
#   "DGE_pseudobulk"  -- CALIBRATED, recommended primary (manuscript)
#   "DGE_MAST"        -- exploratory; broader gene lists, inflated significance
de_method   <- "DGE_MAST"

de_base     <- file.path(project_dir, "Differential_Expression", de_method)
pa_base     <- file.path(project_dir, "Pathway_Analysis_EnrichR", de_method)
dir.create(pa_base, recursive = TRUE, showWarnings = FALSE)

# Discover input DGE directories from script 05's output structure
input_dirs <- list()
# per-lineage
for (d in list.dirs(file.path(de_base, "by_lineage"), recursive = FALSE)) {
  input_dirs[[paste0("Lineage_", basename(d))]] <- d
}
# per-cluster
for (d in list.dirs(file.path(de_base, "by_cluster"), recursive = FALSE)) {
  input_dirs[[basename(d)]] <- d
}

# ============================ Databases ======================================
tf_databases <- c(
  "TRRUST_Transcription_Factors_2019",
  "ChEA_2022",
  "TRANSFAC_and_JASPAR_PWMs"
)
pathway_databases_human <- c(
  "GO_Biological_Process_2023",
  "MSigDB_Hallmark_2020",
  "Reactome_2022",
  "BioPlanet_2019",
  "Panther_2016"
)
pathway_databases_mouse <- c(
  "KEGG_2021_Mouse",
  "WikiPathways_2024_Mouse"
)
mirna_databases <- c(
  "miRTarBase_2017",          # experimentally validated miRNA-target pairs
  "TargetScan_microRNA_2017"  # seed-based predicted targets
)
# (verify availability on your system with listEnrichrDbs() if newer versions exist)

# ============================ Helpers ========================================
normalize_combined <- function(df) {
  if ("Combined Score" %in% colnames(df))
    df <- df %>% rename(Combined.Score = `Combined Score`)
  df
}

# collect top-N significant terms per database, tagging Database name
collect_top <- function(enrich_list, top_n = 10) {
  out <- list()
  for (db in names(enrich_list)) {
    r <- enrich_list[[db]]
    if (!is.data.frame(r) || nrow(r) == 0) next
    r <- normalize_combined(r)
    if (!all(c("Combined.Score","Adjusted.P.value") %in% colnames(r))) next
    sig <- r %>% filter(Adjusted.P.value < 0.05) %>%
      arrange(desc(Combined.Score)) %>% slice_head(n = top_n) %>%
      mutate(Database = db)
    if (nrow(sig) > 0) out[[db]] <- sig
  }
  if (length(out) == 0) return(NULL)
  bind_rows(out)
}

plot_enrichment <- function(df, title, out_png, top_n = 20) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  df <- df %>% arrange(desc(Combined.Score)) %>% slice_head(n = top_n)
  p <- ggplot(df, aes(x = reorder(Term, Combined.Score),
                      y = Combined.Score, fill = Database)) +
    geom_bar(stat = "identity") +
    scale_y_log10() +
    coord_flip() +
    labs(title = title, x = NULL, y = "log10(Combined Score)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10, color = "black"),
          axis.title  = element_text(size = 13),
          plot.title  = element_text(size = 14, hjust = 0.5))
  ggsave(out_png, plot = p, width = 12, height = 10, dpi = 300, bg = "white")
}

write_xlsx_groups <- function(enrich_groups, out_xlsx) {
  wb <- createWorkbook()
  for (grp_name in names(enrich_groups)) {
    enr <- enrich_groups[[grp_name]]
    if (length(enr) == 0) next
    for (db in names(enr)) {
      sheet <- substr(paste0(grp_name, "_", db), 1, 31)
      addWorksheet(wb, sheet)
      writeData(wb, sheet, enr[[db]])
    }
  }
  saveWorkbook(wb, out_xlsx, overwrite = TRUE)
}

# Flag any "mir-29" terms in miRNA enrichment results
flag_mir29 <- function(enr_list, set_label) {
  rows <- list()
  for (db in names(enr_list)) {
    r <- enr_list[[db]]
    if (!is.data.frame(r) || nrow(r) == 0) next
    hits <- r[grepl("mir.?29|MIR29", r$Term, ignore.case = TRUE), ]
    if (nrow(hits) > 0) {
      hits$Database <- db
      hits$gene_set <- set_label
      rows[[db]] <- hits
    }
  }
  if (length(rows) == 0) return(NULL)
  bind_rows(rows)
}

# ============================ Main enrichment ================================
run_enrichment <- function(dge_df, label, out_dir, sig_thresh = 0.05) {
  if (!all(c("p_val_adj","avg_log2FC") %in% colnames(dge_df))) {
    message("  - ", label, ": missing p_val_adj/avg_log2FC; skipping")
    return(invisible(NULL))
  }
  sig <- dge_df %>% filter(!is.na(p_val_adj), p_val_adj < sig_thresh)
  if (nrow(sig) == 0) { message("  - ", label, ": no sig genes"); return(invisible(NULL)) }
  
  genes_all  <- rownames(sig)
  genes_up   <- rownames(sig %>% filter(avg_log2FC > 0))
  genes_down <- rownames(sig %>% filter(avg_log2FC < 0))
  
  message("  - ", label, ": ", length(genes_all), " sig (",
          length(genes_up), " up / ", length(genes_down), " down)")
  
  # display title: turn 00_Activated_intermediate_CD4_miR29a_vs_Scr into
  # "00 Activated intermediate CD4 - miR29a vs Scr" (keep _vs_ as " vs ")
  pretty <- gsub("_vs_", " vs ", label, fixed = TRUE)
  pretty <- gsub("_",      " ",    pretty,  fixed = TRUE)
  
  csv_dir  <- file.path(out_dir, "CSVs")
  plot_dir <- file.path(out_dir, "Plots")
  dir.create(csv_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- TF + pathway enrichment (on full sig list) -------------------------
  enr_tf      <- enrichr(toupper(genes_all), tf_databases)
  enr_path_h  <- enrichr(toupper(genes_all), pathway_databases_human)
  enr_path_m  <- enrichr(genes_all,          pathway_databases_mouse)  # mouse symbols
  
  # ---- miRNA target enrichment (up/down separately) -----------------------
  enr_mirna_down <- if (length(genes_down) > 0)
    enrichr(toupper(genes_down), mirna_databases) else list()
  enr_mirna_up   <- if (length(genes_up)   > 0)
    enrichr(toupper(genes_up),   mirna_databases) else list()
  
  # ---- write workbook with all sheets -------------------------------------
  write_xlsx_groups(
    list(TF       = enr_tf,
         PathwayH = enr_path_h,
         PathwayM = enr_path_m,
         miRNAdn  = enr_mirna_down,
         miRNAup  = enr_mirna_up),
    file.path(csv_dir, paste0(label, "_Enrichment.xlsx"))
  )
  
  # ---- plots --------------------------------------------------------------
  plot_enrichment(collect_top(enr_tf),
                  paste("Top TFs --", pretty),
                  file.path(plot_dir, paste0(label, "_TF.png")))
  plot_enrichment(collect_top(c(enr_path_h, enr_path_m)),
                  paste("Top pathways --", pretty),
                  file.path(plot_dir, paste0(label, "_Pathways.png")))
  plot_enrichment(collect_top(enr_mirna_down),
                  paste("Top miRNAs (DOWN-reg genes) --", pretty),
                  file.path(plot_dir, paste0(label, "_miRNA_DOWN.png")))
  plot_enrichment(collect_top(enr_mirna_up),
                  paste("Top miRNAs (UP-reg genes) --", pretty),
                  file.path(plot_dir, paste0(label, "_miRNA_UP.png")))
  
  # ---- highlight: any miR-29* hits? --------------------------------------
  mir29_hits <- bind_rows(
    flag_mir29(enr_mirna_down, "DOWN"),
    flag_mir29(enr_mirna_up,   "UP")
  )
  if (!is.null(mir29_hits) && nrow(mir29_hits) > 0) {
    write.csv(mir29_hits,
              file.path(csv_dir, paste0(label, "_miR29_hits.csv")),
              row.names = FALSE)
    message("    >>> miR-29 hits: ", nrow(mir29_hits),
            " (see ", label, "_miR29_hits.csv)")
  }
}

# ============================ Main loop ======================================
for (group_name in names(input_dirs)) {
  in_dir  <- input_dirs[[group_name]]
  out_dir <- file.path(pa_base, group_name)
  message("\n=== ", group_name, " ===")
  
  csvs <- list.files(in_dir,
                     pattern = "_(miR29a_vs_Scr|miR29a_vs_EV|EV_vs_Scr)\\.csv$",
                     full.names = TRUE)
  if (length(csvs) == 0) { message("  no DE CSVs found"); next }
  
  for (csv in csvs) {
    label <- tools::file_path_sans_ext(basename(csv))
    dge <- tryCatch(read.csv(csv, row.names = 1, check.names = FALSE),
                    error = function(e) { message("  read fail: ", e$message); NULL })
    if (is.null(dge)) next
    run_enrichment(dge, label, out_dir)
  }
}

message("\nDone. See ", pa_base)
###############################################################################

