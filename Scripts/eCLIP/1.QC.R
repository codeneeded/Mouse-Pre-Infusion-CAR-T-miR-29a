## =====================================================================
## 01_load_qc.R  —  STEP 1: Load all layers + QC (self-contained)
##
## Ground truth = the DATA FILES. The Eclipsebio HTML report is used only for
## read-level metrics not present in the BEDs (depth, trim, mapping, PCR-dup,
## % chimeras); its cluster counts are shown beside file-derived counts as a
## cross-check (they do NOT match — see data_vs_report_check/).
##
## Output layout (flat):
##   eCLIP/01_QC/QC1..QC8_*.png   the eight QC figures (PNG only)
##   eCLIP/01_QC/tables/*.csv     all CSV tables (upload these back)
##   eCLIP/01_QC/README_01_QC.txt plain-text explainer
## Objects (qs2) -> saved_R_data/.
## =====================================================================

## ---- libraries (plain loads) -------------------------------------------
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(rvest)
library(ggrepel)
library(qs2)

## =====================================================================
## CONFIG (self-contained — edit these paths)
## =====================================================================
DATA_DIR   <- "/media/akshay-iyer/Elements/data_from_hpc/QN-0000916_miR-eCLIP_No-Gel/files"
REPORT_DIR <- "/media/akshay-iyer/Elements/data_from_hpc/QN-0000916_miR-eCLIP_No-Gel/deliverables"
SUMMARY_REPORT <- file.path(REPORT_DIR, "miR_eCLIP_summary_report.html")
ROOT_OUT <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/eCLIP"
OBJ_DIR  <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/saved_R_data"
## Genome/annotation (from report): UCSC GRCm38/mm10, GENCODE vM25. Used in step 3.

STATES      <- c("Tn","Teff","Tmem","Tex")
STATE_LABEL <- c(Tn="Naive", Teff="Effector", Tmem="Memory", Tex="Exhausted")
state_factor <- function(x) factor(x, levels = STATES)
STATE_COLS  <- c(Tn="#3B6FB6", Teff="#E08A3C", Tmem="#4FA168", Tex="#C0413B")
FEATURE_COLS <- c("3' UTR"="#1F77B4","5' UTR"="#17BECF","CDS"="#FF7F0E",
                  "Intron"="#2CA02C","miRNA"="#9467BD","ncRNA"="#8C564B","Other"="#7F7F7F")

SAMPLES <- data.table(
  ip    = c("IP1_Teff_1_RR2","IP3_Teff_2_RR2","IP2_Tex_1_RR2","IP4_Tex_2_RR2",
            "IP5_Tmem_1_RR2","IP6_Tmem_2_RR2","IP7_Tn_1_2_RR2","IP8_Tn_3_RR2"),
  state = c("Teff","Teff","Tex","Tex","Tmem","Tmem","Tn","Tn"),
  rep   = c(1L,2L,1L,2L,1L,2L,1L,2L))
SAMPLES[, state := state_factor(state)]

f_persample_clusters <- function(ip) file.path(DATA_DIR, paste0(ip, ".mir_targets_clusters.bed"))
f_repro_clusters     <- function(st) file.path(DATA_DIR, paste0(st, ".mir_targets_reproducible_clusters.bed"))
f_ago2_peaks         <- function(st) file.path(DATA_DIR, paste0(st, ".ago2_reproducible_peaks.bed"))
f_mt_compare    <- function(a,b) file.path(DATA_DIR, sprintf("%s_vs_%s.miRNA_target_sample_comparisons.final.tsv", a, b))
f_eclip_compare <- function(a,b) file.path(DATA_DIR, sprintf("%s_vs_%s.eclip_sample_comparisons.final.tsv", a, b))
CONTRASTS <- list(c("Teff","Tex"),c("Teff","Tmem"),c("Teff","Tn"),
                  c("Tex","Tmem"),c("Tex","Tn"),c("Tmem","Tn"))

num_clean <- function(x) as.numeric(gsub(",", "", gsub("%", "", trimws(as.character(x)))))
read_cluster_bed <- function(path, type = c("mir","ago2")) {
  type <- match.arg(type); dt <- fread(path, header = FALSE, sep = "\t")
  if (type == "mir") setnames(dt, 1:6, c("chr","start","end","miRNA","coverage","strand"))
  else               setnames(dt, 1:6, c("chr","start","end","repl_freq","log2fc_enrich","strand"))
  dt[, start := start + 1L]; dt[]
}
theme_pub <- function(base_size = 13) theme_minimal(base_size) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
        axis.line  = element_line(linewidth = 0.4, colour = "grey20"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey40"),
        plot.title = element_text(face = "bold", size = rel(1.05)),
        plot.subtitle = element_text(colour = "grey35", size = rel(0.85)),
        strip.text = element_text(face = "bold"), legend.key.size = unit(0.8,"lines"),
        plot.margin = margin(10,14,10,10))

## ---- output scaffolding: flat 01_QC/ + single tables/ folder ----------
QC_ROOT <- file.path(ROOT_OUT, "01_QC")          # all QC PNGs land here
TAB     <- file.path(QC_ROOT, "tables")          # all CSVs land here
dir.create(TAB,     recursive = TRUE, showWarnings = FALSE)
dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)
save_png <- function(p, name, w = 8, h = 5.5)
  ggsave(file.path(QC_ROOT, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
csv      <- function(x, f) fwrite(x, file.path(TAB, f))
save_obj <- function(x, n) { qs2::qs_save(x, file.path(OBJ_DIR, paste0(n, ".qs2"))); message("saved ", n) }

## =====================================================================
## 2.  DATA-DERIVED metrics (ground truth — straight from the files)
## =====================================================================
clusters_persample <- rbindlist(lapply(SAMPLES$ip, function(ip) {
  dt <- read_cluster_bed(f_persample_clusters(ip), "mir"); dt[, ip := ip]; dt }))
clusters_persample <- merge(clusters_persample, SAMPLES[, .(ip, state, rep)], by = "ip")

clusters_repro <- rbindlist(lapply(STATES, function(st) {
  dt <- read_cluster_bed(f_repro_clusters(st), "mir"); dt[, state := st]; dt }))
clusters_repro[, state := state_factor(state)]

ago2_peaks <- rbindlist(lapply(STATES, function(st) {
  dt <- read_cluster_bed(f_ago2_peaks(st), "ago2"); dt[, state := st]; dt }))
ago2_peaks[, state := state_factor(state)]

obs_persample <- clusters_persample[, .(observed_chimeric_clusters = .N,
                                        observed_miRNAs = uniqueN(miRNA),
                                        median_coverage = as.numeric(median(coverage))),
                                    by = .(ip, state, rep)]
repro_yield <- clusters_repro[, .(reproducible_clusters = .N,
                                  unique_miRNAs = uniqueN(miRNA)), by = state]
csv(repro_yield, "reproducible_yield_per_state.csv")

feat_layer <- function(files, lay) {
  dt <- rbindlist(lapply(files, function(f) fread(f, sep = "\t", select = "Feature")))
  dt <- dt[Feature != "Feature"]; dt[, .N, by = Feature][, layer := lay] }
feature_dist <- rbind(
  feat_layer(sapply(CONTRASTS, function(c) f_mt_compare(c[1], c[2])),    "chimeric (miR target)"),
  feat_layer(sapply(CONTRASTS, function(c) f_eclip_compare(c[1], c[2])), "AGO2 occupancy"))
setnames(feature_dist, "Feature", "feature"); feature_dist[, frac := N/sum(N), by = layer]
csv(feature_dist, "feature_distribution.csv")

## =====================================================================
## 3.  REPORT metrics (read-level only) + reconcile against the data
## =====================================================================
if (!file.exists(SUMMARY_REPORT)) {
  hit <- list.files(REPORT_DIR, pattern = "summary_report.*\\.html$", full.names = TRUE, recursive = TRUE)
  if (length(hit)) SUMMARY_REPORT <- hit[1]
}
have_report <- file.exists(SUMMARY_REPORT)
qc_report <- NULL
if (have_report) {
  tabs <- rvest::html_table(rvest::read_html(SUMMARY_REPORT), fill = TRUE)
  for (t in tabs) if ("Initial reads" %in% names(t)) { qc_report <- as.data.table(t); break }
  setnames(qc_report,
           c("Library","Initial reads","% Pass trim","% Repetitive elements","% Uniquely aligned",
             "% PCR duplicates","Final nonchimeric reads","AGO2 clusters","AGO2 peaks",
             "Final chimeric reads","% Chimeras","Chimeric clusters"),
           c("library","initial_reads","pct_pass_trim","pct_repetitive","pct_unique_aln",
             "pct_pcr_dup","final_nonchimeric","ago2_clusters","ago2_peaks",
             "final_chimeric","pct_chimeras","report_chimeric_clusters"), skip_absent = TRUE)
  qc_report <- qc_report[grepl("^IP|^Input", library)]
  nc <- setdiff(names(qc_report), "library"); qc_report[, (nc) := lapply(.SD, num_clean), .SDcols = nc]
  qc_report[, kind  := fifelse(grepl("^IP", library), "IP", "Input")]
  qc_report[, state := state_factor(str_extract(library, "Teff|Tex|Tmem|Tn"))]
  csv(qc_report[kind=="IP"], "qc_report_metrics_IP.csv")
  csv(qc_report[kind=="Input", 1:7, with=FALSE], "qc_report_metrics_input.csv")
} else {
  warning("Summary report not found under ", REPORT_DIR,
          " — read-level QC skipped; data-derived QC still runs.")
}

recon <- merge(obs_persample,
               if (have_report) qc_report[kind=="IP", .(ip=library, report_chimeric_clusters,
                                                        final_chimeric, pct_chimeras, ago2_peaks)]
               else obs_persample[, .(ip)][, report_chimeric_clusters := NA_real_],
               by = "ip", all.x = TRUE)
recon[, diff_pct := 100*(observed_chimeric_clusters - report_chimeric_clusters)/report_chimeric_clusters]
csv(recon, "data_vs_report_reconciliation.csv")

## =====================================================================
## 4.  PUBLICATION FIGURES (PNG only) — saved into topic folders
## =====================================================================
if (have_report) {
  qc_report[, library := factor(library, levels = library[order(state, kind)])]
  
  p1 <- ggplot(qc_report, aes(library, initial_reads/1e6, fill = state, alpha = kind)) +
    geom_col(width=.72, colour="grey25", linewidth=.2) +
    scale_fill_manual(values=STATE_COLS, labels=STATE_LABEL) +
    scale_alpha_manual(values=c(IP=1, Input=.45)) + coord_flip() +
    labs(title="Sequencing depth per library", x=NULL, y="initial reads (millions)",
         fill="state", alpha=NULL) + theme_pub()
  save_png(p1, "QC1_sequencing_depth", 8, 6)
  
  ratecols <- c("pct_pass_trim","pct_unique_aln","pct_pcr_dup","pct_repetitive")
  ratelab  <- c(pct_pass_trim="% pass trim", pct_unique_aln="% uniquely aligned",
                pct_pcr_dup="% PCR duplicates", pct_repetitive="% repetitive elements")
  ratel <- melt(qc_report, id.vars=c("library","state","kind"), measure.vars=ratecols,
                variable.name="metric", value.name="pct")
  ratel[, metric := factor(ratelab[as.character(metric)], levels=ratelab)]
  p2 <- ggplot(ratel, aes(library, pct, fill=state, alpha=kind)) +
    geom_col(width=.72, colour="grey25", linewidth=.2) + facet_wrap(~metric, nrow=1) +
    scale_fill_manual(values=STATE_COLS, guide="none") +
    scale_alpha_manual(values=c(IP=1, Input=.45)) + coord_flip() +
    labs(title="Library processing & mapping rates", x=NULL, y="percent", alpha=NULL) + theme_pub(11)
  save_png(p2, "QC2_processing_rates", 13, 6)
  
  qIP <- qc_report[kind=="IP"]
  p3 <- ggplot(qIP, aes(reorder(library, pct_chimeras), pct_chimeras, fill=state)) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=.5, ymax=1.5, fill="grey80", alpha=.35) +
    geom_col(width=.7, colour="grey25", linewidth=.2) +
    geom_text(aes(label=sprintf("%.2f%%", pct_chimeras)), hjust=-.15, size=3.1) +
    scale_fill_manual(values=STATE_COLS, labels=STATE_LABEL) + coord_flip() +
    expand_limits(y=max(qIP$pct_chimeras)*1.15) +
    labs(title="Chimeric capture per IP library",
         subtitle="shaded band ~ typical miR-eCLIP range (0.5-1.5% of mapped reads)",
         x=NULL, y="% chimeric reads", fill="state") + theme_pub()
  save_png(p3, "QC3_chimeric_rate", 8.5, 5.5)
  
  p5 <- ggplot(qIP, aes(ago2_clusters/1e3, ago2_peaks/1e3, colour=state)) +
    geom_point(size=3.4) + geom_text_repel(aes(label=library), size=3, show.legend=FALSE) +
    scale_colour_manual(values=STATE_COLS, labels=STATE_LABEL) +
    labs(title="AGO2 occupancy: clusters vs significant peaks",
         x="AGO2 clusters (thousands)", y="AGO2 peaks (thousands)", colour="state") + theme_pub()
  save_png(p5, "QC5_ago2_clusters_vs_peaks", 8, 6)
  
  rc <- melt(recon[, .(ip, state, observed=observed_chimeric_clusters,
                       report=report_chimeric_clusters)],
             id.vars=c("ip","state"), variable.name="source", value.name="clusters")
  p4 <- ggplot(rc, aes(reorder(ip, clusters), clusters, fill=source)) +
    geom_col(position=position_dodge(width=.75), width=.7, colour="grey25", linewidth=.2) +
    scale_fill_manual(values=c(observed="#2C3E50", report="#BDC3C7"),
                      labels=c(observed="BED files (analysis input)", report="Eclipsebio report")) +
    coord_flip() + scale_y_log10(labels=label_comma()) +
    labs(title="Chimeric clusters: data files vs report",
         subtitle="downstream analysis uses the BED rows; report summary runs ~30% lower",
         x=NULL, y="chimeric clusters (log10)", fill=NULL) + theme_pub()
  save_png(p4, "QC4_data_vs_report_clusters", 8.5, 5.5)
}

ry <- melt(repro_yield, id.vars="state")
p6 <- ggplot(ry, aes(state, value, fill=state)) +
  geom_col(width=.7, colour="grey25", linewidth=.2) +
  geom_text(aes(label=comma(value)), vjust=-.3, size=3.1) +
  facet_wrap(~variable, scales="free_y",
             labeller=as_labeller(c(reproducible_clusters="reproducible chimeric clusters",
                                    unique_miRNAs="distinct engaged miRNAs"))) +
  scale_fill_manual(values=STATE_COLS, guide="none") +
  scale_x_discrete(labels=STATE_LABEL[STATES]) +
  labs(title="Reproducible chimeric yield per CD8 state (from files)", x=NULL, y=NULL) +
  theme_pub() + expand_limits(y=0)
save_png(p6, "QC6_reproducible_yield", 9, 5)

mir_cov <- clusters_persample[, .(cov=sum(coverage)), by=.(state, rep, miRNA)]
conc <- dcast(mir_cov, state + miRNA ~ rep, value.var="cov", fill=0)
setnames(conc, c("1","2"), c("rep1","rep2"))
conc_stats <- conc[, .(spearman=cor(log1p(rep1), log1p(rep2), method="spearman"),
                       n_miRNA=.N), by=state]
csv(conc, "persample_miRNA_coverage_rep1_rep2.csv")
csv(conc_stats, "replicate_concordance_stats.csv")
p7 <- ggplot(conc, aes(rep1+1, rep2+1)) +
  geom_point(aes(colour=state), alpha=.5, size=1.6, show.legend=FALSE) +
  geom_abline(slope=1, intercept=0, lty=2, colour="grey50") +
  geom_text(data=conc_stats, aes(x=2, y=max(conc$rep2)+1,
                                 label=sprintf("rho = %.2f\nn = %d", spearman, n_miRNA)),
            hjust=0, vjust=1, size=3.2) +
  facet_wrap(~state, labeller=labeller(state=STATE_LABEL)) +
  scale_x_log10(labels=label_comma()) + scale_y_log10(labels=label_comma()) +
  scale_colour_manual(values=STATE_COLS) +
  labs(title="Replicate concordance of per-miRNA target coverage (from files)",
       x="replicate 1 coverage (+1, log10)", y="replicate 2 coverage (+1, log10)") + theme_pub()
save_png(p7, "QC7_replicate_concordance", 8.5, 7)

p8 <- ggplot(feature_dist, aes(reorder(feature, -frac), frac, fill=feature)) +
  geom_col(width=.72, colour="grey25", linewidth=.2) +
  geom_text(aes(label=sprintf("%.0f%%", 100*frac)), vjust=-.3, size=3) +
  facet_wrap(~layer) + scale_fill_manual(values=FEATURE_COLS, guide="none") +
  scale_y_continuous(labels=percent_format()) +
  labs(title="Where AGO2 binds vs where miRNAs make chimeras (from files)",
       subtitle="chimeric targets are not 3'UTR-dominated - CDS/intron carry most signal",
       x=NULL, y="fraction of sites") +
  theme_pub() + theme(axis.text.x=element_text(angle=25, hjust=1))
save_png(p8, "QC8_feature_distribution", 10, 5.5)

## =====================================================================
## 5.  Save objects
## =====================================================================
save_obj(clusters_persample, "clusters_persample")
save_obj(clusters_repro,     "clusters_repro")
save_obj(ago2_peaks,         "ago2_peaks")
save_obj(feature_dist,       "feature_dist")
save_obj(repro_yield,        "repro_yield")
save_obj(recon,              "qc_reconciliation")
if (have_report) save_obj(qc_report, "qc_report_metrics")

## =====================================================================
## 6.  Plain-text explainer
## =====================================================================
L <- c(
  "=====================================================================",
  " STEP 1 - QUALITY CONTROL : miR-eCLIP CD8 differentiation",
  "=====================================================================",
  "",
  "GROUND TRUTH = THE DATA FILES. The Eclipsebio HTML report supplies only",
  "read-level metrics absent from the BEDs (depth, trim, mapping, PCR-dup,",
  "% chimeras). Its cluster counts are a summary and DO NOT match the files.",
  "",
  "GENOME: UCSC GRCm38/mm10, GENCODE vM25 (differs from CAR-T's GRCm39 - step 3",
  "annotation must use a vM25/mm10 GTF).",
  "",
  "DESIGN: 4 CD8 states (Naive Teff Tmem Tex), 2 IP + 2 input each. AGO2 IP +",
  "ligation -> chimeric reads = DIRECT miR->target pairs (a small % of reads).",
  "",
  "OUTPUT (eCLIP/01_QC/): QC1..QC8 PNGs in the folder; all CSVs in 01_QC/tables/.",
  "",
  "--- DATA vs REPORT: chimeric clusters (per IP) ----------------------",
  "  The .mir_targets_clusters.bed files carry ~30% MORE rows than the report's",
  "  'Chimeric clusters' summary. Analysis runs on the BED rows; treat report",
  "  values as approximate. Per library (observed in file | report | diff):",
  apply(recon[order(state)], 1, function(r) sprintf("    %-16s %6s | %6s | %+5.0f%%",
                                                    r[["ip"]], r[["observed_chimeric_clusters"]],
                                                    ifelse(is.na(r[["report_chimeric_clusters"]]),"NA",r[["report_chimeric_clusters"]]),
                                                    as.numeric(r[["diff_pct"]]))),
  "",
  "--- REPRODUCIBLE YIELD PER STATE (from files) ----------------------",
  apply(repro_yield, 1, function(r) sprintf("    %-5s  %6s reproducible clusters | %4s distinct miRNAs",
                                            r[["state"]], r[["reproducible_clusters"]], r[["unique_miRNAs"]])),
  "",
  "--- REPLICATE CONCORDANCE (Spearman, per-miRNA coverage) -----------",
  apply(conc_stats, 1, function(r) sprintf("    %-5s  rho = %.2f (n=%s miRNAs)",
                                           r[["state"]], as.numeric(r[["spearman"]]), r[["n_miRNA"]])),
  "",
  "--- KEY FLAGS ------------------------------------------------------",
  "  * NAIVE IS UNBALANCED: IP7_Tn_1_2 >> IP8_Tn_3 in chimeric reads/clusters.",
  "    Raw per-library counts are not comparable. Defenses: reproducible set",
  "    (must be in BOTH reps) + within-state normalization (step 2).",
  "  * Tex has the lowest chimeric yield - absolute target counts into exhaustion",
  "    read low partly for capture reasons, not only biology.",
  "  * Chimeric targets span CDS/intron/3'UTR: do NOT restrict to 3'UTR downstream.",
  "",
  "--- BAM CHECK (IP1_Teff_1 .mir_targets.bam inspected) --------------",
  "  The chimeric BAM is STAR-aligned target reads (mm10) with a UMI in the read",
  "  name and NO miRNA tag - miRNA identity is not in the BAM. The cluster BED",
  "  'coverage' column IS the definitive per-miRNA read support; no BAM needed",
  "  for aim 1. (BAM read count = report 'Final chimeric reads', a useful check.)",
  "",
  "NEXT (step 2, aim 1): per-miRNA engaged abundance per state, normalized for the",
  "  naive imbalance, from the reproducible set + per-sample coverage.",
  "=====================================================================")
writeLines(L, file.path(QC_ROOT, "README_01_QC.txt"))
message("01_load_qc COMPLETE -> ", QC_ROOT)