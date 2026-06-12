###############################################################################
# 07b -- Responder vs Non-Responder module scoring (CD4 / CD8 subsets)        #
#        + hypoxia / ROS / stress readout on CD4 & CD8 subsets and clusters    #
# Mouse pre-infusion CAR-T scRNA-seq, miR-29a project                         #
#-----------------------------------------------------------------------------#
# STANDALONE EXTENSION to script 07. It does NOT recompute the base modules;   #
# it loads script 07's checkpoint (Mouse_CARTmiR29a_WithModuleScores.qs2),     #
# adds the human CAR-T Responder/Non-Responder signatures (converted human ->  #
# mouse), and renders the score FeaturePlots on the CD4 and CD8 lineage        #
# SUBSETS (subset, NOT subclustered -- existing umap.harmony coords reused).   #
#                                                                              #
# Wetlab request (NK_AI email):                                                #
#   * add module scoring from the human Responder vs Non-Responder gene sets   #
#     (CAR-T_Responder_v_Non-Responder_Meta.xlsx); Fraietta 2018 is the        #
#     strongest reference and is treated as the headline set.                  #
#   * plot the module-score FeaturePlots on CD4 and CD8 subsets separately.    #
#   * report whether hypoxia / ROS / stress modules light up on the CD8        #
#     subset, the CD4 subset, or specific clusters.                            #
#                                                                              #
# Human -> mouse conversion: babelgene orthologs (offline, multi-DB support),  #
#   with a validated title-case fallback for unmapped symbols that ONLY keeps  #
#   a candidate if it actually exists in the object (no invented genes).       #
#                                                                              #
# Requires once:  install.packages("babelgene")                                #
###############################################################################

library(Seurat); library(qs2); library(readxl); library(dplyr); library(tidyr)
library(ggplot2); library(patchwork); library(viridis)
library(SeuratExtend); library(scCustomize); library(babelgene)
# ============================ Paths (reused from script 07) ==================
project_dir   <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/10x_scRNAseq"
saved_dir     <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/saved_R_data"
resources_dir <- file.path(project_dir, "Resources")

# NEW input: the human Responder vs Non-Responder workbook.
rvnr_xlsx     <- file.path(resources_dir, "CAR-T_Responder_v_Non-Responder_Meta.xlsx")

# NEW outputs, kept separate from script 07's Module_Scores tree.
out_base      <- file.path(project_dir, "Module_Scores", "Responder_NonResponder")
plot_dir      <- file.path(out_base, "Plots")
umap_cd4_dir  <- file.path(plot_dir, "UMAP_CD4")
umap_cd8_dir  <- file.path(plot_dir, "UMAP_CD8")
vln_cd4_dir   <- file.path(plot_dir, "VlnByCluster_CD4")
vln_cd8_dir   <- file.path(plot_dir, "VlnByCluster_CD8")
stress_dir    <- file.path(plot_dir, "Hypoxia_ROS_Stress")
data_dir      <- file.path(out_base, "Tables")
for (d in c(plot_dir, umap_cd4_dir, umap_cd8_dir, vln_cd4_dir, vln_cd8_dir,
            stress_dir, data_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ============================ Load script-07 checkpoint ======================
# We deliberately load the WithModuleScores object so the hypoxia/ROS/stress
# scores computed in script 07 are already present and can be re-read here.
chk <- file.path(saved_dir, "Mouse_CARTmiR29a_WithModuleScores.qs2")
if (!file.exists(chk))
  stop("Missing ", chk,
       "\nRun script 07 first -- this extension reads its module scores.")
obj <- qs_read(chk)

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")
if (!"data" %in% Layers(obj[["RNA"]]))
  obj <- NormalizeData(obj, verbose = FALSE)

obj$clusters          <- droplevels(factor(obj$clusters))
obj$tentative_lineage <- droplevels(factor(obj$tentative_lineage))
obj$tentative_state   <- droplevels(factor(obj$tentative_state))
obj$condition         <- factor(obj$condition, levels = c("EV", "Scr", "miR29a"))

umap_reduction <- "umap.harmony"
if (!umap_reduction %in% Reductions(obj))
  umap_reduction <- Reductions(obj)[grepl("umap", Reductions(obj),
                                          ignore.case = TRUE)][1]

# ============================ Resolve CD4 / CD8 lineages =====================
# "subset, not subclustered" -- we select cells by lineage label and reuse the
# existing embedding. Match labels flexibly so naming variants are handled.
lin_levels <- levels(obj$tentative_lineage)
message("tentative_lineage levels present: ", paste(lin_levels, collapse = ", "))
cd4_levels <- grep("CD4", lin_levels, ignore.case = TRUE, value = TRUE)
cd8_levels <- grep("CD8", lin_levels, ignore.case = TRUE, value = TRUE)
if (length(cd4_levels) == 0 || length(cd8_levels) == 0)
  stop("Could not find CD4/CD8 lineage labels in tentative_lineage.\n",
       "Levels seen: ", paste(lin_levels, collapse = ", "),
       "\nEdit cd4_levels / cd8_levels to match your labels.")
message("CD4 lineage labels: ", paste(cd4_levels, collapse = ", "))
message("CD8 lineage labels: ", paste(cd8_levels, collapse = ", "))

# ============================ Read R-vs-NR workbook ==========================
# Long sheet is the clean single source (set_name | gene_symbol | Group | Study)
long <- read_excel(rvnr_xlsx, sheet = "ALL _by Publication")
stopifnot(all(c("set_name", "gene_symbol", "Group") %in% colnames(long)))
long$gene_symbol <- trimws(as.character(long$gene_symbol))
long <- long[!is.na(long$gene_symbol) & long$gene_symbol != "", ]
# Relabel the ICI2025 set_name (this is the Shorer 2025 dataset) for clarity.
long$set_name <- gsub("ICI2025", "Shorer2025", long$set_name)

human_sets <- split(unique(long$gene_symbol), long$set_name)  # named list

# De-duplicated combined R / NR (wide sheet).
comb <- read_excel(rvnr_xlsx, sheet = "ALL_combined, Dups removed")
human_sets$Combined_R  <- unique(trimws(na.omit(as.character(comb[[1]]))))
human_sets$Combined_NR <- unique(trimws(na.omit(as.character(comb[[2]]))))

# Headline sets get composite grids + condition tests; the rest still get
# individual subset FeaturePlots (Fraietta = strongest reference per wetlab).
headline_sets <- c("Fraietta2018_R", "Fraietta2018_NR",
                   "Combined_R",     "Combined_NR")

# ============================ Human -> mouse conversion ======================
present_mouse <- rownames(obj[["RNA"]])

title_case_gene <- function(x) {
  # human ALL-CAPS symbol -> mouse Title-case (CXCL13 -> Cxcl13). Symbols with
  # separators (HLA-DQA1) are passed through but will simply fail to match.
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
}

convert_to_mouse <- function(genes, set_label) {
  genes <- unique(genes[!is.na(genes) & genes != ""])
  # also offer a de-dotted form (R made.names artefacts: HLA.DQB1 -> HLA-DQB1)
  genes_query <- unique(c(genes, gsub("\\.", "-", genes)))
  
  orth <- tryCatch(
    babelgene::orthologs(genes = genes_query, species = "mouse", human = TRUE),
    error = function(e) { message("  babelgene error: ", conditionMessage(e)); NULL })
  
  mapped <- data.frame(human = character(0), mouse = character(0),
                       support_n = numeric(0), method = character(0),
                       stringsAsFactors = FALSE)
  if (!is.null(orth) && nrow(orth) > 0) {
    # best mouse ortholog per human symbol (highest support_n)
    best <- orth %>%
      dplyr::group_by(human_symbol) %>%
      dplyr::slice_max(order_by = support_n, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    mapped <- data.frame(human = best$human_symbol, mouse = best$symbol,
                         support_n = best$support_n, method = "babelgene",
                         stringsAsFactors = FALSE)
  }
  
  mapped_human <- toupper(mapped$human)
  unmapped <- genes[!toupper(genes) %in% mapped_human]
  
  # validated title-case fallback -- keep ONLY if present in the object
  if (length(unmapped) > 0) {
    cand <- title_case_gene(unmapped)
    keep <- which(cand %in% present_mouse)
    if (length(keep) > 0)
      mapped <- rbind(mapped, data.frame(
        human = unmapped[keep], mouse = cand[keep], support_n = NA_real_,
        method = "title_case_validated", stringsAsFactors = FALSE))
  }
  
  mapped$set_name   <- set_label
  mapped$in_object  <- mapped$mouse %in% present_mouse
  mapped
}

message("\nConverting human -> mouse for ", length(human_sets), " gene sets...")
conv_all <- dplyr::bind_rows(lapply(names(human_sets), function(s)
  convert_to_mouse(human_sets[[s]], s)))

# Conversion / presence audit (one row per human gene per set).
write.csv(conv_all, file.path(data_dir, "RvNR_human_to_mouse_audit.csv"),
          row.names = FALSE)

# per-set conversion summary
conv_summary <- conv_all %>%
  dplyr::group_by(set_name) %>%
  dplyr::summarise(
    n_human_input   = length(unique(human_sets[[dplyr::first(set_name)]])),
    n_mapped        = dplyr::n(),
    n_in_object     = sum(in_object),
    n_babelgene     = sum(method == "babelgene"),
    n_titlecase_fb  = sum(method == "title_case_validated"),
    .groups = "drop")
write.csv(conv_summary, file.path(data_dir, "RvNR_conversion_summary.csv"),
          row.names = FALSE)
message("Conversion audit written: RvNR_human_to_mouse_audit.csv (+ summary)")

# Build mouse module lists (genes present in object, >=3 to score).
module_lists <- lapply(split(conv_all[conv_all$in_object, ], 
                             conv_all$set_name[conv_all$in_object]),
                       function(df) unique(df$mouse))
module_lists <- module_lists[lengths(module_lists) >= 3]
if (length(module_lists) == 0)
  stop("No R/NR set retained >=3 mouse genes present in the object.")

modules_used <- do.call(rbind, lapply(names(module_lists), function(m)
  data.frame(module = m, gene = module_lists[[m]], stringsAsFactors = FALSE)))
write.csv(modules_used, file.path(data_dir, "RvNR_module_genes_used.csv"),
          row.names = FALSE)
message("Mouse module gene lists written: RvNR_module_genes_used.csv")

# ============================ AddModuleScore (full object) ===================
# Scores are computed once on the FULL object so control-gene binning spans the
# whole landscape; we only SUBSET afterwards for plotting.
message("Computing R/NR module scores: ", length(module_lists), " modules")
for (m in names(module_lists)) {
  obj <- AddModuleScore(obj, features = list(module_lists[[m]]),
                        name = paste0(m, "__"), assay = "RNA",
                        ctrl = max(10, min(100, length(module_lists[[m]]) * 5)))
  raw <- paste0(m, "__1"); new <- paste0("score_", m)
  obj@meta.data[[new]] <- obj@meta.data[[raw]]
  obj@meta.data[[raw]] <- NULL
  # display mirror column (clean plot titles -- same trick as script 07)
  obj@meta.data[[m]] <- obj@meta.data[[new]]
}

# Net responder score (R - NR) for headline datasets: a single signed readout
# of "responder-like vs non-responder-like" per cell.
make_net <- function(prefix) {
  rcol <- paste0("score_", prefix, "_R"); ncol <- paste0("score_", prefix, "_NR")
  if (all(c(rcol, ncol) %in% colnames(obj@meta.data))) {
    nm <- paste0("net_", prefix)
    obj@meta.data[[nm]] <<- obj@meta.data[[rcol]] - obj@meta.data[[ncol]]
    nm
  } else NULL
}
net_cols <- Filter(Negate(is.null), lapply(c("Fraietta2018", "Combined"), make_net))
net_cols <- unlist(net_cols)

# ============================ Build CD4 / CD8 subsets ========================
# subset() retains the umap.harmony reduction for the kept cells -> we plot on
# the ORIGINAL coordinates (no re-clustering, no re-embedding).
cd4_cells <- colnames(obj)[obj$tentative_lineage %in% cd4_levels]
cd8_cells <- colnames(obj)[obj$tentative_lineage %in% cd8_levels]
obj_cd4 <- subset(obj, cells = cd4_cells)
obj_cd8 <- subset(obj, cells = cd8_cells)
for (o in c("obj_cd4", "obj_cd8")) {
  tmp <- get(o)
  tmp$tentative_state <- droplevels(factor(tmp$tentative_state))
  tmp$clusters        <- droplevels(factor(tmp$clusters))
  assign(o, tmp)
}
message("CD4 subset: ", ncol(obj_cd4), " cells | CD8 subset: ", ncol(obj_cd8))

# ============================ Plot helpers (script-07 style) =================
pal  <- viridis(n = 10, option = "A")               # magma, matches script 07
safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

plot_score_umap <- function(object, feat, ttl_suffix) {
  FeaturePlot_scCustom(object, reduction = umap_reduction, features = feat,
                       colors_use = pal, order = TRUE) +
    ggtitle(paste0(feat, "  (", ttl_suffix, ")"))
}

plot_score_vln_cluster <- function(object, feat) {
  VlnPlot2(object, features = feat, group.by = "tentative_state",
           cols = "default", show.mean = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          plot.margin = margin(t = 5, r = 5, b = 40, l = 40, unit = "pt"))
}

# signed net score -> diverging palette centred at zero
plot_net_umap <- function(object, feat, ttl_suffix) {
  FeaturePlot(object, reduction = umap_reduction, features = feat, order = TRUE) +
    scale_color_gradient2(low = "#1B4965", mid = "grey92", high = "#F50057",
                          midpoint = 0) +
    ggtitle(paste0(feat, "  (", ttl_suffix, ")"))
}

# ============================ Per-module subset FeaturePlots =================
message("Saving R/NR FeaturePlots on CD4 and CD8 subsets...")
for (m in names(module_lists)) {
  ggsave(file.path(umap_cd4_dir, paste0(safe(m), ".png")),
         plot_score_umap(obj_cd4, m, "CD4 subset"),
         width = 8, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(umap_cd8_dir, paste0(safe(m), ".png")),
         plot_score_umap(obj_cd8, m, "CD8 subset"),
         width = 8, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(vln_cd4_dir, paste0(safe(m), ".png")),
         plot_score_vln_cluster(obj_cd4, m),
         width = 14, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(vln_cd8_dir, paste0(safe(m), ".png")),
         plot_score_vln_cluster(obj_cd8, m),
         width = 14, height = 8, dpi = 300, bg = "white")
}

# net responder score plots (headline)
for (nc in net_cols) {
  ggsave(file.path(umap_cd4_dir, paste0(safe(nc), ".png")),
         plot_net_umap(obj_cd4, nc, "CD4 subset"),
         width = 8, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(umap_cd8_dir, paste0(safe(nc), ".png")),
         plot_net_umap(obj_cd8, nc, "CD8 subset"),
         width = 8, height = 7, dpi = 300, bg = "white")
}

# ============================ Composite headline grids =======================
hl_present <- intersect(headline_sets, names(module_lists))
if (length(hl_present) > 0) {
  g4 <- wrap_plots(lapply(hl_present, function(m) plot_score_umap(obj_cd4, m, "CD4")),
                   ncol = 2)
  g8 <- wrap_plots(lapply(hl_present, function(m) plot_score_umap(obj_cd8, m, "CD8")),
                   ncol = 2)
  ggsave(file.path(plot_dir, "UMAP_headline_RvNR_CD4.png"), g4,
         width = 16, height = 7 * ceiling(length(hl_present) / 2),
         dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(file.path(plot_dir, "UMAP_headline_RvNR_CD8.png"), g8,
         width = 16, height = 7 * ceiling(length(hl_present) / 2),
         dpi = 300, bg = "white", limitsize = FALSE)
}

# ============================ Condition test within subsets ==================
# Exploratory two-sided Wilcoxon: does miR-29a shift the responder/NR signature
# vs each control, within CD4 and within CD8 (overall + per cluster)?
headline_score_cols <- c(paste0("score_", hl_present), net_cols)
run_condition_tests <- function(object, lineage_label) {
  md <- object@meta.data
  res <- list()
  for (sc in intersect(headline_score_cols, colnames(md))) {
    levels_to_test <- c("OVERALL", levels(droplevels(factor(md$tentative_state))))
    for (lv in levels_to_test) {
      sub <- if (lv == "OVERALL") md else md[md$tentative_state == lv, ]
      for (ctrl in c("EV", "Scr")) {
        x <- sub[[sc]][sub$condition == "miR29a"]
        y <- sub[[sc]][sub$condition == ctrl]
        if (length(x) > 10 && length(y) > 10) {
          w <- wilcox.test(x, y)               # two-sided, exploratory
          res[[length(res) + 1]] <- data.frame(
            lineage = lineage_label, module = sc, level = lv,
            group = paste0("miR29a_vs_", ctrl),
            median_miR29a = median(x), median_ctrl = median(y),
            delta = median(x) - median(y),
            direction = ifelse(median(x) > median(y), "up_in_miR29a", "down_in_miR29a"),
            W = unname(w$statistic), p_value = w$p.value,
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  dplyr::bind_rows(res)
}
cond_df <- dplyr::bind_rows(run_condition_tests(obj_cd4, "CD4"),
                            run_condition_tests(obj_cd8, "CD8"))
if (nrow(cond_df) > 0) {
  cond_df <- cond_df %>%
    dplyr::group_by(lineage, module) %>%
    dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup()
  write.csv(cond_df,
            file.path(data_dir, "RvNR_module_condition_wilcoxon_CD4_CD8.csv"),
            row.names = FALSE)
  message("Condition tests written: RvNR_module_condition_wilcoxon_CD4_CD8.csv")
}

###############################################################################
# ============================ Hypoxia / ROS / Stress readout ================
# Answers the wetlab question: do these modules pop up on the CD8 subset, the
# CD4 subset, or specific clusters? Uses scores ALREADY in the object from 07.
###############################################################################
existing_scores <- grep("^score_", colnames(obj@meta.data), value = TRUE)
existing_scores <- setdiff(existing_scores, paste0("score_", names(module_lists)))
stress_pat <- "hypox|oxidative|reactive.?oxygen|(^|_)ros($|_)|stress|unfold|heat.?shock"
stress_scores <- existing_scores[grepl(stress_pat, existing_scores, ignore.case = TRUE)]

if (length(stress_scores) == 0) {
  message("\n[Hypoxia/ROS/Stress] No matching module found in the object.")
  message("  score_ columns available: ",
          paste(sub("^score_", "", existing_scores), collapse = ", "))
  message("  -> If these modules should exist, add them to ",
          "Module_Gene_Lists.xlsx and re-run script 07.")
} else {
  message("\n[Hypoxia/ROS/Stress] modules detected: ",
          paste(sub("^score_", "", stress_scores), collapse = ", "))
  
  # per-cluster, per-lineage summary using a global z so values are comparable
  full_md <- obj@meta.data
  summ <- list()
  for (sc in stress_scores) {
    z_full <- as.numeric(scale(full_md[[sc]]))            # z across ALL cells
    full_md[[paste0(sc, "_z")]] <- z_full
    for (lin_lab in c("CD4", "CD8")) {
      lin_set <- if (lin_lab == "CD4") cd4_levels else cd8_levels
      keep <- full_md$tentative_lineage %in% lin_set
      sub  <- full_md[keep, ]
      for (cl in levels(droplevels(factor(sub$tentative_state)))) {
        v  <- sub[[sc]][sub$tentative_state == cl]
        vz <- sub[[paste0(sc, "_z")]][sub$tentative_state == cl]
        summ[[length(summ) + 1]] <- data.frame(
          module = sub("^score_", "", sc), lineage = lin_lab, cluster = cl,
          n_cells = length(v), mean_score = mean(v), median_score = median(v),
          mean_z = mean(vz), pct_z_gt0 = mean(vz > 0) * 100,
          elevated = mean(vz) > 0.25, stringsAsFactors = FALSE)
      }
      # whole-lineage row
      v  <- sub[[sc]]; vz <- sub[[paste0(sc, "_z")]]
      summ[[length(summ) + 1]] <- data.frame(
        module = sub("^score_", "", sc), lineage = lin_lab, cluster = "ALL",
        n_cells = length(v), mean_score = mean(v), median_score = median(v),
        mean_z = mean(vz), pct_z_gt0 = mean(vz > 0) * 100,
        elevated = mean(vz) > 0.25, stringsAsFactors = FALSE)
    }
  }
  stress_summary <- dplyr::bind_rows(summ)
  write.csv(stress_summary,
            file.path(data_dir, "hypoxia_ROS_stress_subset_cluster_summary.csv"),
            row.names = FALSE)
  
  # CD8-vs-CD4 lineage-level Wilcoxon per module (which lineage is higher)
  lin_cmp <- lapply(stress_scores, function(sc) {
    x <- obj@meta.data[[sc]][obj$tentative_lineage %in% cd8_levels]
    y <- obj@meta.data[[sc]][obj$tentative_lineage %in% cd4_levels]
    w <- wilcox.test(x, y)
    data.frame(module = sub("^score_", "", sc),
               median_CD8 = median(x), median_CD4 = median(y),
               delta_CD8_minus_CD4 = median(x) - median(y),
               higher_in = ifelse(median(x) > median(y), "CD8", "CD4"),
               p_value = w$p.value, stringsAsFactors = FALSE)
  }) %>% dplyr::bind_rows()
  write.csv(lin_cmp,
            file.path(data_dir, "hypoxia_ROS_stress_CD8_vs_CD4.csv"),
            row.names = FALSE)
  
  # FeaturePlots + per-cluster violins on each subset
  for (sc in stress_scores) {
    feat <- sub("^score_", "", sc)
    if (!feat %in% colnames(obj@meta.data)) obj@meta.data[[feat]] <- obj@meta.data[[sc]]
    obj_cd4@meta.data[[feat]] <- obj@meta.data[colnames(obj_cd4), feat]
    obj_cd8@meta.data[[feat]] <- obj@meta.data[colnames(obj_cd8), feat]
    ggsave(file.path(stress_dir, paste0(safe(feat), "_CD4_umap.png")),
           plot_score_umap(obj_cd4, feat, "CD4 subset"),
           width = 8, height = 7, dpi = 300, bg = "white")
    ggsave(file.path(stress_dir, paste0(safe(feat), "_CD8_umap.png")),
           plot_score_umap(obj_cd8, feat, "CD8 subset"),
           width = 8, height = 7, dpi = 300, bg = "white")
    ggsave(file.path(stress_dir, paste0(safe(feat), "_CD4_vln.png")),
           plot_score_vln_cluster(obj_cd4, feat),
           width = 14, height = 8, dpi = 300, bg = "white")
    ggsave(file.path(stress_dir, paste0(safe(feat), "_CD8_vln.png")),
           plot_score_vln_cluster(obj_cd8, feat),
           width = 14, height = 8, dpi = 300, bg = "white")
  }
  
  # ---- human-readable console answer ----
  message("\n================= HYPOXIA / ROS / STRESS READOUT =================")
  for (mod in unique(stress_summary$module)) {
    lc <- lin_cmp[lin_cmp$module == mod, ]
    message(sprintf("* %s: higher overall in %s (median %.3f vs %.3f, p=%.2e)",
                    mod, lc$higher_in, lc$median_CD8, lc$median_CD4, lc$p_value))
    top <- stress_summary %>%
      dplyr::filter(module == mod, cluster != "ALL") %>%
      dplyr::arrange(dplyr::desc(mean_z)) %>% head(3)
    for (i in seq_len(nrow(top)))
      message(sprintf("    top cluster: %s/%s  mean_z=%.2f  (%d cells)%s",
                      top$lineage[i], top$cluster[i], top$mean_z[i], top$n_cells[i],
                      ifelse(top$elevated[i], "  [elevated]", "")))
  }
  message("==================================================================")
  message("Tables: hypoxia_ROS_stress_subset_cluster_summary.csv, ",
          "hypoxia_ROS_stress_CD8_vs_CD4.csv")
}

# ============================ Re-save checkpoint =============================
# New checkpoint; does NOT overwrite script 07's WithModuleScores.qs2.
qs_save(obj, file = file.path(saved_dir,
                              "Mouse_CARTmiR29a_WithModuleScores_RvNR.qs2"))

# ============================ README =========================================
writeLines(c(
  "07b -- Responder vs Non-Responder module scores (CD4/CD8 subsets)",
  "",
  "Inputs:",
  "  saved_R_data/Mouse_CARTmiR29a_WithModuleScores.qs2  (script 07 output)",
  "  Resources/CAR-T_Responder_v_Non-Responder_Meta.xlsx",
  "",
  "Human->mouse conversion: babelgene orthologs + validated title-case fallback.",
  "Scores computed on the FULL object, then plotted on CD4 / CD8 subsets",
  "(subset only -- existing umap.harmony coords reused, no re-clustering).",
  "",
  "Plots/UMAP_CD4, Plots/UMAP_CD8      per-set score FeaturePlots + net_* scores",
  "Plots/VlnByCluster_CD4 / _CD8       per-cluster violins within each lineage",
  "Plots/UMAP_headline_RvNR_CD4/CD8    Fraietta + Combined composite grids",
  "Plots/Hypoxia_ROS_Stress            hypoxia/ROS/stress FeaturePlots + violins",
  "",
  "Tables/RvNR_human_to_mouse_audit.csv         per-gene mapping + method + in_object",
  "Tables/RvNR_conversion_summary.csv           per-set mapping counts",
  "Tables/RvNR_module_genes_used.csv            mouse genes actually scored",
  "Tables/RvNR_module_condition_wilcoxon_CD4_CD8.csv  miR29a vs EV/Scr per subset",
  "Tables/hypoxia_ROS_stress_subset_cluster_summary.csv",
  "Tables/hypoxia_ROS_stress_CD8_vs_CD4.csv"),
  con = file.path(out_base, "README.txt"))

message("\nDone. Outputs under: ", out_base)
###############################################################################

