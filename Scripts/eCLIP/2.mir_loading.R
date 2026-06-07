## =====================================================================
## 02_miR_loading.R  —  STEP 2 / AIM 1: per-miRNA engaged abundance per state
##
## Question: which miRNAs are abundantly ENGAGED (loaded + actively targeting)
## at each CD8 state, and how does loading shift across differentiation?
##
## Quantity: chimeric-cluster COVERAGE (the per-miRNA read support; the BAM
## carries no miRNA identity, so coverage IS the quantification — see step-1 QC).
## Universe: miRNAs in the REPRODUCIBLE set of >=1 state (engaged, both-rep).
## Normalization: DESeq2 size factors on per-library summed coverage -> removes
## the naive IP7 depth spike. n=2/state, so differential loading is EXPLORATORY.
##
## Output (flat): eCLIP/02_miR_loading/ M1..M5 PNGs + tables/ + README.
## Objects (qs2) -> saved_R_data/.
## =====================================================================

## ---- libraries (plain loads) -------------------------------------------
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(ggrepel)
library(DESeq2)
library(qs2)

## =====================================================================
## CONFIG (self-contained — edit paths)
## =====================================================================
DATA_DIR <- "/media/akshay-iyer/Elements/data_from_hpc/QN-0000916_miR-eCLIP_No-Gel/files"
ROOT_OUT <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/eCLIP"
OBJ_DIR  <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/saved_R_data"

STATES      <- c("Tn","Teff","Tmem","Tex")
STATE_LABEL <- c(Tn="Naive", Teff="Effector", Tmem="Memory", Tex="Exhausted")
state_factor <- function(x) factor(x, levels = STATES)
STATE_COLS  <- c(Tn="#3B6FB6", Teff="#E08A3C", Tmem="#4FA168", Tex="#C0413B")

SAMPLES <- data.table(
  ip    = c("IP1_Teff_1_RR2","IP3_Teff_2_RR2","IP2_Tex_1_RR2","IP4_Tex_2_RR2",
            "IP5_Tmem_1_RR2","IP6_Tmem_2_RR2","IP7_Tn_1_2_RR2","IP8_Tn_3_RR2"),
  state = c("Teff","Teff","Tex","Tex","Tmem","Tmem","Tn","Tn"),
  rep   = c(1L,2L,1L,2L,1L,2L,1L,2L))
SAMPLES[, state := state_factor(state)]

f_persample_clusters <- function(ip) file.path(DATA_DIR, paste0(ip, ".mir_targets_clusters.bed"))
f_repro_clusters     <- function(st) file.path(DATA_DIR, paste0(st, ".mir_targets_reproducible_clusters.bed"))
read_cluster_bed <- function(path) {            # 6-col, no header, mm10
  dt <- fread(path, header = FALSE, sep = "\t")
  setnames(dt, 1:6, c("chr","start","end","miRNA","coverage","strand"))
  dt[, start := start + 1L]; dt[]
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

OUT <- file.path(ROOT_OUT, "02_miR_loading"); TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)
save_png <- function(p, name, w = 8, h = 5.5)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
csv <- function(x, f) fwrite(x, file.path(TAB, f))
save_obj <- function(x, n) { qs2::qs_save(x, file.path(OBJ_DIR, paste0(n, ".qs2"))); message("saved ", n) }

## curated immunology miRNAs to track across the differentiation axis
CURATED <- c("mmu-miR-155-5p","mmu-miR-150-5p","mmu-miR-146a-5p","mmu-miR-21a-5p",
             "mmu-miR-29a-3p","mmu-miR-29b-3p","mmu-miR-29c-3p",
             "mmu-miR-17-5p","mmu-miR-19b-3p","mmu-miR-92a-3p",
             "mmu-let-7a-5p","mmu-let-7b-5p","mmu-let-7c-5p",
             "mmu-miR-31-5p","mmu-miR-210-3p","mmu-miR-142a-3p","mmu-miR-142a-5p",
             "mmu-miR-101a-3p","mmu-miR-23a-3p","mmu-miR-16-5p","mmu-miR-9-5p")

## =====================================================================
## 1.  Load layers (self-contained, from raw files)
## =====================================================================
clusters_persample <- rbindlist(lapply(SAMPLES$ip, function(ip) {
  dt <- read_cluster_bed(f_persample_clusters(ip)); dt[, ip := ip]; dt }))
clusters_persample <- merge(clusters_persample, SAMPLES[, .(ip, state, rep)], by = "ip")

clusters_repro <- rbindlist(lapply(STATES, function(st) {
  dt <- read_cluster_bed(f_repro_clusters(st)); dt[, state := st]; dt }))
clusters_repro[, state := state_factor(state)]

## engaged universe = miRNAs reproducible in >=1 state
engaged <- sort(unique(clusters_repro$miRNA))
message(length(engaged), " engaged (reproducible) miRNAs")

## =====================================================================
## 2.  Per-library coverage count matrix (engaged universe) + DESeq2 norm
## =====================================================================
ps <- clusters_persample[miRNA %in% engaged, .(cov = sum(coverage)), by = .(ip, miRNA)]
cmat <- dcast(ps, miRNA ~ ip, value.var = "cov", fill = 0)
M <- as.matrix(cmat[, -1]); rownames(M) <- cmat$miRNA
M <- round(M)                                          # integer pseudo-counts (coverage ~ depth)
M <- M[, SAMPLES$ip]                                   # order columns to sample sheet
M <- M[rowSums(M >= 1) >= 2, , drop = FALSE]           # light prefilter

coldata <- data.frame(state = SAMPLES$state, rep = SAMPLES$rep, row.names = SAMPLES$ip)
dds <- DESeqDataSetFromMatrix(M, coldata, design = ~ state)
dds <- DESeq(dds)
normM <- counts(dds, normalized = TRUE)                # depth-normalized coverage
logN  <- log2(normM + 1)

## abundance per state = mean log2-normalized coverage; + within-state relative %
ab_long <- as.data.table(reshape2::melt(logN, varnames = c("miRNA","ip"), value.name = "logN"))
ab_long <- merge(ab_long, SAMPLES[, .(ip, state)], by = "ip")
abundance <- ab_long[, .(mean_logN = mean(logN)), by = .(miRNA, state)]
## within-state relative abundance from normalized linear scale
lin <- as.data.table(reshape2::melt(normM, varnames=c("miRNA","ip"), value.name="norm"))
lin <- merge(lin, SAMPLES[, .(ip, state)], by="ip")
relab <- lin[, .(mean_norm = mean(norm)), by=.(miRNA,state)][
  , rel_pct := 100*mean_norm/sum(mean_norm), by=state]
abundance <- merge(abundance, relab[, .(miRNA,state,mean_norm,rel_pct)], by=c("miRNA","state"))
abundance[, state := state_factor(state)]
csv(dcast(abundance, miRNA ~ state, value.var = "mean_logN"), "miRNA_abundance_by_state_log2norm.csv")
csv(abundance, "miRNA_abundance_by_state_long.csv")

## =====================================================================
## 3.  M1 — top engaged miRNAs per state (within-state relative %)
## =====================================================================
topN <- 15
top_state <- abundance[order(-rel_pct), head(.SD, topN), by = state]
csv(top_state, "top_miRNAs_per_state.csv")
top_state[, lab := factor(paste(miRNA, state, sep="@@"),
                          levels = rev(paste(miRNA, state, sep="@@")))]
m1 <- ggplot(top_state, aes(lab, rel_pct, fill = state)) +
  geom_col(width = .74, colour = "grey25", linewidth = .2) +
  coord_flip() +
  scale_x_discrete(labels = function(z) sub("@@.*$", "", z)) +
  scale_fill_manual(values = STATE_COLS, guide = "none") +
  facet_wrap(~ state, scales = "free_y", labeller = labeller(state = STATE_LABEL)) +
  labs(title = "Top engaged miRNAs per CD8 state",
       subtitle = "within-state relative abundance (DESeq2-normalized chimeric coverage)",
       x = NULL, y = "% of state engaged signal") + theme_pub(11)
save_png(m1, "M1_top_miRNAs_per_state", 11, 8)

## =====================================================================
## 4.  M2 — state x miRNA heatmap (z-scored, most state-variable engaged miR)
## =====================================================================
wide <- dcast(abundance, miRNA ~ state, value.var = "mean_logN")
WM <- as.matrix(wide[, ..STATES]); rownames(WM) <- wide$miRNA
vary <- order(apply(WM, 1, function(r) max(r) - min(r)), decreasing = TRUE)
WMv <- WM[head(vary, 45), , drop = FALSE]
Z <- t(scale(t(WMv)))                                  # row z-score
ord <- hclust(dist(Z))$order
zl <- as.data.table(reshape2::melt(Z[ord, ], varnames = c("miRNA","state"), value.name = "z"))
zl[, miRNA := factor(miRNA, levels = rownames(Z)[ord])]
zl[, state := factor(state, levels = STATES)]
m2 <- ggplot(zl, aes(state, miRNA, fill = z)) +
  geom_tile(colour = "white", linewidth = .3) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "grey96", high = "#C0413B", midpoint = 0) +
  scale_x_discrete(labels = STATE_LABEL[STATES]) +
  labs(title = "State-variable engaged miRNAs", subtitle = "row z-score of log2-normalized coverage",
       x = NULL, y = NULL, fill = "z") +
  theme_pub(10) + theme(axis.text.y = element_text(size = 7), panel.grid = element_blank())
save_png(m2, "M2_state_miRNA_heatmap", 6.5, 9)

## =====================================================================
## 5.  M3 — state-specificity (tau index) + peak state
## =====================================================================
tau_of <- function(x) { x <- pmax(x, 0); if (max(x) <= 0) return(NA_real_); sum(1 - x/max(x))/(length(x)-1) }
spec <- abundance[, .(tau = tau_of(mean_norm),
                      peak_state = STATES[which.max(mean_norm)],
                      max_norm = max(mean_norm)), by = miRNA]
spec[, peak_state := state_factor(peak_state)]
csv(spec[order(-tau)], "miRNA_state_specificity.csv")
## show specific + reasonably abundant miRNAs
sp_show <- spec[max_norm >= median(spec$max_norm)][order(-tau)][1:25]
sp_show[, miRNA := factor(miRNA, levels = rev(miRNA))]
m3 <- ggplot(sp_show, aes(miRNA, tau, fill = peak_state)) +
  geom_col(width = .74, colour = "grey25", linewidth = .2) + coord_flip() +
  scale_fill_manual(values = STATE_COLS, labels = STATE_LABEL) +
  labs(title = "Most state-specific engaged miRNAs",
       subtitle = "tau = 1 fully state-restricted; bar colour = state of peak loading",
       x = NULL, y = "specificity (tau)", fill = "peaks in") + theme_pub(11)
save_png(m3, "M3_state_specificity", 8.5, 7)

## =====================================================================
## 6.  M4 — differential loading across states (DESeq2; n=2 EXPLORATORY)
## =====================================================================
diff_contrasts <- list(c("Teff","Tn"), c("Tmem","Tn"), c("Tex","Tn"),
                       c("Tex","Tmem"), c("Tex","Teff"))
res_all <- rbindlist(lapply(diff_contrasts, function(cc) {
  r <- as.data.table(results(dds, contrast = c("state", cc[1], cc[2])), keep.rownames = "miRNA")
  r[, contrast := paste0(cc[1], "_vs_", cc[2])][]
}))
res_all[, sig := !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1]
csv(res_all[order(contrast, padj)], "differential_loading_all_contrasts.csv")
lab4 <- res_all[sig == TRUE | miRNA %in% CURATED]
m4 <- ggplot(res_all, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(colour = sig), size = .8, alpha = .6) +
  scale_colour_manual(values = c(`TRUE` = "#C0413B", `FALSE` = "grey75"), guide = "none") +
  geom_vline(xintercept = c(-1, 1), lty = 2, colour = "grey60") +
  geom_hline(yintercept = -log10(0.05), lty = 2, colour = "grey60") +
  geom_text_repel(data = lab4, aes(label = sub("mmu-", "", miRNA)), size = 2.4, max.overlaps = 14) +
  facet_wrap(~ contrast, nrow = 2, scales = "free") +
  labs(title = "Differential miRNA loading across states (DESeq2)",
       subtitle = "log2FC > 0 = higher in first state. n=2/state: treat as exploratory ranking",
       x = "log2 fold change", y = "-log10 adjusted p") + theme_pub(11)
save_png(m4, "M4_differential_loading_volcano", 11, 7.5)

## =====================================================================
## 7.  M5 — curated immunology miRNAs across the differentiation axis
## =====================================================================
cur <- abundance[miRNA %in% CURATED]
cur[, miRNA := factor(miRNA, levels = CURATED[CURATED %in% cur$miRNA])]
csv(dcast(cur, miRNA ~ state, value.var = "mean_logN"), "curated_immune_miRNAs.csv")
m5 <- ggplot(cur, aes(state, mean_logN, group = miRNA)) +
  geom_line(colour = "grey70", linewidth = .4) +
  geom_point(aes(colour = state), size = 2.6) +
  scale_colour_manual(values = STATE_COLS, guide = "none") +
  scale_x_discrete(labels = STATE_LABEL[STATES]) +
  facet_wrap(~ miRNA, scales = "free_y") +
  labs(title = "Curated immunology miRNAs across CD8 differentiation",
       subtitle = "log2 DESeq2-normalized chimeric coverage (engaged loading)",
       x = NULL, y = "log2 norm coverage") +
  theme_pub(10) + theme(strip.text = element_text(size = 8))
save_png(m5, "M5_curated_immune_miRNAs", 12, 8)

## =====================================================================
## 8.  Save objects
## =====================================================================
save_obj(dds,        "loading_dds")
save_obj(abundance,  "miR_abundance_by_state")
save_obj(res_all,    "miR_differential_loading")
save_obj(spec,       "miR_state_specificity")

## =====================================================================
## 9.  README explainer
## =====================================================================
n_eng <- length(engaged); n_test <- nrow(M)
topfmt <- function(st) {
  d <- top_state[state == st][order(-rel_pct)][1:5]
  paste0("    ", st, ": ", paste(sprintf("%s (%.1f%%)", sub("mmu-","",d$miRNA), d$rel_pct), collapse = ", "))
}
sig_n <- res_all[sig == TRUE, .N, by = contrast]
L <- c(
  "=====================================================================",
  " STEP 2 / AIM 1 — per-miRNA engaged abundance per CD8 state",
  "=====================================================================",
  "",
  "WHAT THIS ANSWERS: which miRNAs are abundantly engaged (loaded + actively",
  "forming chimeras with targets) at each state, and how loading shifts across",
  "naive -> effector -> memory and into exhaustion.",
  "",
  "QUANTITY & NORMALIZATION:",
  "  - Unit = chimeric-cluster COVERAGE per miRNA (per step-1 QC, the BAM has no",
  "    miRNA identity, so coverage is the definitive per-miRNA read support).",
  sprintf("  - Universe = %d miRNAs reproducible in >=1 state; %d pass the count", n_eng, n_test),
  "    prefilter and are tested.",
  "  - DESeq2 size factors normalize per-library coverage -> removes the naive",
  "    IP7 depth spike so states are comparable.",
  "  - n=2/state: DIFFERENTIAL LOADING IS EXPLORATORY (ranking aid, not proof).",
  "",
  "--- TOP ENGAGED miRNAs PER STATE (within-state %) ------------------",
  topfmt("Tn"), topfmt("Teff"), topfmt("Tmem"), topfmt("Tex"),
  "",
  "--- SIGNIFICANT DIFFERENTIAL LOADING (padj<0.05 & |log2FC|>=1) ------",
  if (nrow(sig_n)) apply(sig_n, 1, function(r) sprintf("    %-12s %s miRNAs", r[["contrast"]], r[["N"]]))
  else "    (none cleared the threshold — expected with n=2; use ranking)",
  "",
  "--- HOW TO READ THE FIGURES ----------------------------------------",
  "  M1 top engaged miRNAs per state (what dominates each state's signal)",
  "  M2 z-scored heatmap of the most state-variable miRNAs (loading shifts)",
  "  M3 tau state-specificity (which miRNAs are restricted to one state)",
  "  M4 differential loading volcanoes across 5 contrasts (exploratory)",
  "  M5 curated immunology miRNAs across the axis (incl. miR-29 family for the",
  "     CAR-T link; miR-155/-31/-210 exhaustion; miR-150/let-7 naive/quiescence)",
  "",
  "--- BIOLOGY TO CHECK -----------------------------------------------",
  "  * miR-150 / let-7 high in naive (quiescence) and falling on activation?",
  "  * miR-155 / miR-21 / miR-31 / miR-210 rising into exhaustion?",
  "  * miR-17~92 (miR-17/19b/92a) up in effector (proliferation)?",
  "  * miR-29a/b/c levels — feeds directly into the CAR-T miR-29a overexpression",
  "    work (its targets are validated in step 6).",
  "  NOTE: naive abundance is dominated by few miRNAs (e.g. miR-466i); interpret",
  "  with the IP7/IP8 imbalance from step-1 QC in mind. Tex has the thinnest",
  "  engaged layer, so low Tex abundance is partly capture, not only biology.",
  "",
  "--- CSVs TO UPLOAD (eCLIP/02_miR_loading/tables) -------------------",
  "  miRNA_abundance_by_state_long.csv, top_miRNAs_per_state.csv,",
  "  miRNA_state_specificity.csv, differential_loading_all_contrasts.csv,",
  "  curated_immune_miRNAs.csv",
  "",
  "NEXT (step 3, aim 2 prep): annotate reproducible clusters to genes (GENCODE",
  "  vM25 / mm10) so we know WHAT each miRNA targets per state.",
  "=====================================================================")
writeLines(L, file.path(OUT, "README_02_miR_loading.txt"))
message("02_miR_loading COMPLETE -> ", OUT)