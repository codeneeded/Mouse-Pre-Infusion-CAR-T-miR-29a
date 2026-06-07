## =====================================================================
## 05_differential_targeting.R  —  STEP 5 / AIM 3: how targeting drives differentiation
##
## Primary lens = the Eclipsebio DESeq2 differential TSVs (IP-vs-input
## normalized -> controls for transcript abundance, unlike step-4 raw
## convergence). 6 pairwise contrasts. We focus on:
##   - what targeting is GAINED/LOST per contrast (esp. into exhaustion)
##   - differentiation MASTER REGULATORS that change targeting
##   - TARGET SWITCHING: what the top switchers (miR-210/miR-155) bind per state
##
## Abundant-RNA background (mt-*/Rpl*/Rps*/Mrp*/Actb) is FILTERED for gene-level
## views (per step-4 interpretation).
##
## Output (flat): eCLIP/05_differential/ D1..D5 PNGs + tables/ + README.
## =====================================================================

## ---- libraries (plain loads) -------------------------------------------
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(ggrepel)
library(qs2)

## =====================================================================
## CONFIG (self-contained)
## =====================================================================
DATA_DIR <- "/media/akshay-iyer/Elements/data_from_hpc/QN-0000916_miR-eCLIP_No-Gel/files"
ROOT_OUT <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/eCLIP"
OBJ_DIR  <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/saved_R_data"

STATES      <- c("Tn","Teff","Tmem","Tex")
STATE_LABEL <- c(Tn="Naive", Teff="Effector", Tmem="Memory", Tex="Exhausted")
state_factor <- function(x) factor(x, levels = STATES)
STATE_COLS  <- c(Tn="#3B6FB6", Teff="#E08A3C", Tmem="#4FA168", Tex="#C0413B")

## contrasts present on disk; log2FC sign is (first state vs second)
CONTRASTS <- list(c("Teff","Tn"), c("Tmem","Tn"), c("Tex","Tn"),
                  c("Tex","Tmem"), c("Tex","Teff"), c("Tmem","Teff"))
f_mt_compare <- function(a,b) {                          # files exist in one orientation only
  f1 <- file.path(DATA_DIR, sprintf("%s_vs_%s.miRNA_target_sample_comparisons.final.tsv", a, b))
  f2 <- file.path(DATA_DIR, sprintf("%s_vs_%s.miRNA_target_sample_comparisons.final.tsv", b, a))
  if (file.exists(f1)) list(path=f1, flip=FALSE) else list(path=f2, flip=TRUE)
}

theme_pub <- function(base_size = 13) theme_minimal(base_size) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = .25, colour = "grey88"),
        axis.line = element_line(linewidth = .4, colour = "grey20"),
        axis.ticks = element_line(linewidth = .3, colour = "grey40"),
        plot.title = element_text(face = "bold", size = rel(1.05)),
        plot.subtitle = element_text(colour = "grey35", size = rel(.85)),
        strip.text = element_text(face = "bold"), legend.key.size = unit(.8,"lines"),
        plot.margin = margin(10,14,10,10))

OUT <- file.path(ROOT_OUT, "05_differential"); TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
save_png <- function(p, name, w = 8, h = 5.5)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
csv <- function(x, f) fwrite(x, file.path(TAB, f))
save_obj <- function(x, n) { qs2::qs_save(x, file.path(OBJ_DIR, paste0(n, ".qs2"))); message("saved ", n) }
read_obj <- function(n) qs2::qs_read(file.path(OBJ_DIR, paste0(n, ".qs2")))

MASTER_TF <- c("Tcf7","Lef1","Myb","Bach2","Tbx21","Eomes","Prdm1","Tox",
               "Nr4a1","Nr4a2","Nr4a3","Irf4","Batf","Zeb2","Id2","Id3")
## abundant-RNA background filter
is_background <- function(g) grepl("^(mt-|Rpl|Rps|Mrpl|Mrps|Rn[0-9]|Gm[0-9]+$)", g) | g %in% c("Actb","Malat1","Neat1")

## =====================================================================
## 1.  Load + clean the differential miRNA-target TSVs (input-normalized)
## =====================================================================
clean_cmp <- function(dt, a, b, flip) {
  dt <- data.table::copy(dt)                       # explicit copy -> no shallow-copy warnings
  setnames(dt, c("Gene Name","Gene ID","miRNA"), c("gene","gene_id","miRNA"), skip_absent = TRUE)
  setnames(dt, grep("Fold Change", names(dt), value = TRUE), "log2FC")
  setnames(dt, grep("adjusted p", names(dt), value = TRUE), "neglog10_padj")
  setnames(dt, grep("^-log10\\(p-value", names(dt), value = TRUE), "neglog10_p")
  setnames(dt, c("Chromosome","Start","End","Strand","Feature"),
           c("chr","start","end","strand","feature"), skip_absent = TRUE)
  dt <- dt[gene != "Gene Name"]
  for (cc in c("log2FC","neglog10_padj","neglog10_p","start","end"))
    if (cc %in% names(dt)) dt[[cc]] <- as.numeric(dt[[cc]])
  if (flip) dt[, log2FC := -log2FC]                       # reorient so log2FC = a vs b
  dt[, `:=`(contrast = paste0(a, "_vs_", b), grpA = a, grpB = b)]
  dt[]
}
dtar <- rbindlist(lapply(CONTRASTS, function(cc) {
  fi <- f_mt_compare(cc[1], cc[2])
  clean_cmp(fread(fi$path, sep = "\t"), cc[1], cc[2], fi$flip)
}), fill = TRUE)
dtar[, padj := 10^(-neglog10_padj)]
dtar[, sig := !is.na(neglog10_padj) & neglog10_padj > -log10(0.05) & abs(log2FC) >= 1]
dtar[, background := is_background(gene)]
csv(dtar[sig == TRUE][order(contrast, -abs(log2FC))], "significant_differential_targeting.csv")

## =====================================================================
## 2.  D1 — volcanoes per contrast (master TFs labelled; background dimmed)
## =====================================================================
labv <- dtar[sig == TRUE & !background & gene %in% MASTER_TF]
d1 <- ggplot(dtar[background == FALSE], aes(log2FC, neglog10_padj)) +
  geom_point(aes(colour = sig), size = .7, alpha = .55) +
  scale_colour_manual(values = c(`TRUE` = "#C0413B", `FALSE` = "grey78"), guide = "none") +
  geom_vline(xintercept = c(-1, 1), lty = 2, colour = "grey60") +
  geom_hline(yintercept = -log10(0.05), lty = 2, colour = "grey60") +
  geom_text_repel(data = labv, aes(label = paste0(gene, " / ", sub("mmu-","",miRNA))),
                  size = 2.4, max.overlaps = 16) +
  facet_wrap(~ contrast, nrow = 2, scales = "free") +
  labs(title = "Differential miRNA targeting across states (Eclipsebio DESeq2)",
       subtitle = "log2FC > 0 = stronger targeting in the first state; abundant-RNA background removed",
       x = "log2 fold change (IP-normalized targeting)", y = "-log10 adjusted p") + theme_pub(11)
save_png(d1, "D1_differential_volcano", 12, 8)

## =====================================================================
## 3.  D2 — gained/lost targeting per contrast (directionality)
## =====================================================================
dir_sum <- dtar[sig == TRUE & !background,
                .(gained = sum(log2FC > 0), lost = sum(log2FC < 0)), by = contrast]
dl <- melt(dir_sum, id.vars = "contrast", variable.name = "direction", value.name = "n")
csv(dir_sum, "gained_lost_targeting_per_contrast.csv")
d2 <- ggplot(dl, aes(contrast, ifelse(direction=="lost", -n, n), fill = direction)) +
  geom_col(width = .7, colour = "grey25", linewidth = .2) +
  geom_hline(yintercept = 0, colour = "grey40") +
  scale_fill_manual(values = c(gained = "#C0413B", lost = "#3B6FB6"),
                    labels = c("gained in 1st state","lost in 1st state")) +
  coord_flip() +
  labs(title = "Differential targeting: interactions gained vs lost",
       subtitle = "significant sites (padj<0.05, |log2FC|>=1), background removed",
       x = NULL, y = "# interactions (left = lost, right = gained)", fill = NULL) + theme_pub()
save_png(d2, "D2_gained_lost", 9, 5.5)

## =====================================================================
## 4.  D3 — master-regulator differential targeting heatmap
## =====================================================================
tf <- dtar[gene %in% MASTER_TF & sig == TRUE]
csv(tf[order(gene, contrast)], "masterTF_differential_targeting.csv")
if (nrow(tf)) {
  tf[, pair := paste0(gene, " | ", sub("mmu-","",miRNA))]
  d3 <- ggplot(tf, aes(contrast, pair, fill = log2FC)) +
    geom_tile(colour = "white", linewidth = .4) +
    geom_text(aes(label = sprintf("%.1f", log2FC)), size = 2.6) +
    scale_fill_gradient2(low = "#3B6FB6", mid = "grey95", high = "#C0413B", midpoint = 0) +
    labs(title = "Differential targeting of master regulators",
         subtitle = "log2FC > 0 = stronger in first state of the contrast",
         x = NULL, y = NULL, fill = "log2FC") +
    theme_pub(11) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_png(d3, "D3_masterTF_differential", 9, max(4, .4*nrow(tf)))
}

## =====================================================================
## 5.  D4 — target SWITCHING: what top switchers bind per state
## (uses step-4 interactions_by_state; falls back to a note if missing)
## =====================================================================
switchers <- c("mmu-miR-210-3p","mmu-miR-155-5p","mmu-miR-23a-3p","mmu-miR-150-5p")
if (file.exists(file.path(OBJ_DIR, "interactions_by_state.qs2"))) {
  inter <- as.data.table(read_obj("interactions_by_state"))
  inter[, state := state_factor(as.character(state))]
  sw <- inter[miRNA %in% switchers & !is_background(gene)]
  ## keep, per miRNA, the genes that are among its top targets in any state
  topg <- sw[, .(cov = sum(coverage)), by = .(miRNA, gene)][order(miRNA, -cov)][, head(.SD, 12), by = miRNA]
  swp <- sw[gene %in% topg$gene]
  csv(swp[order(miRNA, gene, state)], "top_switcher_targets_by_state.csv")
  if (nrow(swp)) {
    d4 <- ggplot(swp, aes(state, gene, fill = log2(coverage + 1))) +
      geom_tile(colour = "white", linewidth = .3) +
      facet_wrap(~ sub("mmu-","",miRNA), scales = "free_y", nrow = 1) +
      scale_fill_viridis_c(option = "A", end = .9) +
      scale_x_discrete(labels = STATE_LABEL[STATES]) +
      labs(title = "Target switching: top switcher miRNAs bind different genes per state",
           subtitle = "miR-210 / miR-155 / miR-23a / miR-150 (background removed)",
           x = NULL, y = NULL, fill = "log2 cov") +
      theme_pub(10) + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                            axis.text.y = element_text(size = 7))
    save_png(d4, "D4_target_switching_detail", 13, 7)
  }
} else message("NOTE: run step 4 first for D4 (interactions_by_state.qs2 missing)")

## =====================================================================
## 6.  D5 — exhaustion focus: targeting gained INTO exhaustion
## =====================================================================
exh <- dtar[contrast %in% c("Tex_vs_Tmem","Tex_vs_Teff","Tex_vs_Tn") &
              sig == TRUE & !background & log2FC > 0]
exh_genes <- exh[, .(n = uniqueN(contrast), miRNAs = paste(sort(unique(sub("mmu-","",miRNA))), collapse=","),
                     max_lfc = max(log2FC)), by = gene][order(-n, -max_lfc)]
csv(exh_genes, "exhaustion_gained_targeting_genes.csv")
csv(exh[order(-log2FC)], "exhaustion_gained_targeting_sites.csv")
top_exh <- exh[order(-log2FC)][1:min(25, .N)]
if (nrow(top_exh)) {
  top_exh[, lab := paste0(gene, " / ", sub("mmu-","",miRNA))]
  d5 <- ggplot(top_exh, aes(reorder(lab, log2FC), log2FC, fill = contrast)) +
    geom_col(width = .72, colour = "grey25", linewidth = .2) + coord_flip() +
    scale_fill_brewer(palette = "Reds") +
    labs(title = "Targeting gained into exhaustion",
         subtitle = "miR->target interactions stronger in Tex (vs Tmem/Teff/Tn), background removed",
         x = NULL, y = "log2 fold change", fill = "contrast") + theme_pub(11)
  save_png(d5, "D5_exhaustion_gained_targeting", 9.5, 7)
}

## =====================================================================
## 7.  Save objects + README
## =====================================================================
save_obj(dtar, "differential_targeting")
ntot <- dtar[sig == TRUE & !background, .N]
sig_by <- dtar[sig == TRUE & !background, .(n = .N,
                                            gained = sum(log2FC>0), lost = sum(log2FC<0)), by = contrast]
tf_pairs <- if (nrow(tf)) tf[, .(pair = paste0(gene,"/",sub("mmu-","",miRNA)), contrast, log2FC)] else NULL
L <- c(
  "=====================================================================",
  " STEP 5 / AIM 3 — differential targeting across CD8 states",
  "=====================================================================",
  "",
  "LENS: Eclipsebio DESeq2 differential TSVs (IP-vs-input normalized), which",
  "control for transcript abundance — the clean complement to step-4 raw",
  "convergence. Abundant-RNA background (mt-/Rpl/Rps/Mrp/Actb/Malat1/Neat1)",
  "is removed from gene-level views. log2FC sign = first state vs second.",
  "",
  sprintf("Significant differential interactions (padj<0.05,|log2FC|>=1, non-bg): %d", ntot),
  "",
  "--- GAINED / LOST PER CONTRAST -------------------------------------",
  apply(sig_by, 1, function(r) sprintf("    %-12s n=%-4s gained=%-4s lost=%s",
                                       r[["contrast"]], r[["n"]], r[["gained"]], r[["lost"]])),
  "",
  "--- MASTER-REGULATOR DIFFERENTIAL TARGETING ------------------------",
  if (!is.null(tf_pairs) && nrow(tf_pairs))
    apply(tf_pairs, 1, function(r) sprintf("    %-22s %-12s log2FC %s",
                                           r[["pair"]], r[["contrast"]], round(as.numeric(r[["log2FC"]]),2)))
  else "    (no master-TF interactions cleared threshold in the differential TSVs)",
  "",
  "--- FIGURES --------------------------------------------------------",
  "  D1 volcanoes per contrast (master TFs labelled)",
  "  D2 gained vs lost targeting per contrast (directionality)",
  "  D3 master-regulator differential targeting heatmap",
  "  D4 target switching: top switchers' genes per state (miR-210/-155/-23a/-150)",
  "  D5 targeting gained INTO exhaustion (Tex vs others)",
  "",
  "--- CAVEATS --------------------------------------------------------",
  "  * n=2/state (Eclipsebio DESeq2): treat p-values as exploratory; weight by",
  "    effect size + biological coherence.",
  "  * The differential TSVs are the abundance-controlled view; cross-check",
  "    candidate targets against step-4 reproducible interactions.",
  "  * Tex has lowest capture (step-1 QC): 'lost into exhaustion' can be partly",
  "    technical — interpret losses more cautiously than gains.",
  "",
  "--- CSVs TO UPLOAD (eCLIP/05_differential/tables) ------------------",
  "  significant_differential_targeting.csv, gained_lost_targeting_per_contrast.csv,",
  "  masterTF_differential_targeting.csv, top_switcher_targets_by_state.csv,",
  "  exhaustion_gained_targeting_genes.csv",
  "",
  "NEXT (step 6, aim 4): therapeutic candidates — miRNAs targeting exhaustion",
  "  drivers / effector programs, ranked; miR-29 family revisited at relaxed",
  "  support against the CAR-T DE targets.",
  "=====================================================================")
writeLines(L, file.path(OUT, "README_05_differential.txt"))
message("05_differential_targeting COMPLETE -> ", OUT)

