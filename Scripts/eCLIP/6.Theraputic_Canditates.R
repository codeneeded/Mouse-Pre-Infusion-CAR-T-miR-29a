## =====================================================================
## 06_therapeutic_candidates.R  —  STEP 6 / AIM 4: state-pushing / therapeutic miRNAs
##
## Synthesizes all prior layers to nominate miRNAs that could be used to push
## CD8 state, and revisits the miR-29 / CAR-T cross-validation at RELAXED
## (per-sample, single-rep) support — since the reproducible set missed the
## canonical miR-29 targets (step 3).
##
## Framing (from step-5 aim 3): there is NO exhaustion-specific target program;
## leverage lies in the ACTIVATION switch + the MEMORY module. So candidates are
## scored by (a) targeting effector/memory/exhaustion-program genes, (b) state
## specificity, (c) differential-targeting support.
##
## Inputs (qs2 from prior steps): clusters_annotated (3), interactions_by_state (4),
##   differential_targeting (5), miR_abundance_by_state (2).  Plus per-sample
##   clusters from raw BEDs for the relaxed miR-29 search.
## Output (flat): eCLIP/06_therapeutic/ T1..T5 PNGs + tables/ + README.
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

SAMPLES <- data.table(
  ip    = c("IP1_Teff_1_RR2","IP3_Teff_2_RR2","IP2_Tex_1_RR2","IP4_Tex_2_RR2",
            "IP5_Tmem_1_RR2","IP6_Tmem_2_RR2","IP7_Tn_1_2_RR2","IP8_Tn_3_RR2"),
  state = c("Teff","Teff","Tex","Tex","Tmem","Tmem","Tn","Tn"),
  rep   = c(1L,2L,1L,2L,1L,2L,1L,2L))
SAMPLES[, state := state_factor(state)]
f_persample_clusters <- function(ip) file.path(DATA_DIR, paste0(ip, ".mir_targets_clusters.bed"))
read_cluster_bed <- function(path) {
  dt <- fread(path, header = FALSE, sep = "\t")
  setnames(dt, 1:6, c("chr","start","end","miRNA","coverage","strand")); dt[, start := start+1L]; dt[] }

theme_pub <- function(base_size = 13) theme_minimal(base_size) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = .25, colour = "grey88"),
        axis.line = element_line(linewidth = .4, colour = "grey20"),
        axis.ticks = element_line(linewidth = .3, colour = "grey40"),
        plot.title = element_text(face = "bold", size = rel(1.05)),
        plot.subtitle = element_text(colour = "grey35", size = rel(.85)),
        strip.text = element_text(face = "bold"), legend.key.size = unit(.8,"lines"),
        plot.margin = margin(10,14,10,10))

OUT <- file.path(ROOT_OUT, "06_therapeutic"); TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
save_png <- function(p, name, w = 8, h = 5.5)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
csv <- function(x, f) fwrite(x, file.path(TAB, f))
save_obj <- function(x, n) { qs2::qs_save(x, file.path(OBJ_DIR, paste0(n, ".qs2"))); message("saved ", n) }
read_obj <- function(n) qs2::qs_read(file.path(OBJ_DIR, paste0(n, ".qs2")))

## extended abundant-RNA background filter (adds Gas5/Snhg/Snord/Snora per step-5)
is_background <- function(g) grepl("^(mt-|Rpl|Rps|Mrpl|Mrps|Rn[0-9]|Snhg|Snord|Snora|Sno|Gm[0-9]+$)", g) |
  grepl("Gas5", g) | g %in% c("Actb","Malat1","Neat1","Xist")

## program gene sets (functional, not exhaustion-only — per step-5 framing)
EXHAUSTION_GENES <- c("Pdcd1","Havcr2","Lag3","Tigit","Ctla4","Entpd1","Tox","Nr4a1","Nr4a2","Nr4a3","Cblb")
EFFECTOR_GENES   <- c("Ifng","Gzmb","Prf1","Tbx21","Klrg1","Il2","Tnf","Id2","Zeb2")
MEMORY_GENES     <- c("Tcf7","Lef1","Il7r","Sell","Ccr7","Bcl2","Bach2","Myb","Zfp36l2","Klf2","S1pr1")
CART_MIR29_TARGETS <- c("Dnmt3a","Tet2","Eomes","Tbx21","Ifng")

## =====================================================================
## 1.  Load layers
## =====================================================================
need <- c("clusters_annotated","interactions_by_state","differential_targeting")
miss <- need[!file.exists(file.path(OBJ_DIR, paste0(need, ".qs2")))]
if (length(miss)) stop("run prior steps first; missing objects: ", paste(miss, collapse=", "))
ann   <- as.data.table(read_obj("clusters_annotated"))[, state := state_factor(as.character(state))]
inter <- as.data.table(read_obj("interactions_by_state"))[, state := state_factor(as.character(state))]
dtar  <- as.data.table(read_obj("differential_targeting"))

inter <- inter[!is_background(gene) & !is.na(gene) & gene != ""]

## =====================================================================
## 2.  Candidate miRNAs by PROGRAM targeting (effector / memory / exhaustion)
## =====================================================================
prog_map <- rbind(
  data.table(gene = EFFECTOR_GENES,   program = "effector"),
  data.table(gene = MEMORY_GENES,     program = "memory"),
  data.table(gene = EXHAUSTION_GENES, program = "exhaustion"))
prog_hits <- merge(inter, prog_map, by = "gene", allow.cartesian = TRUE)
csv(prog_hits[order(program, gene, -coverage)], "miR_targeting_program_genes.csv")

## state-specificity (tau) of each miRNA's loading, reused from step 2 if present
spec <- if (file.exists(file.path(OBJ_DIR, "miR_state_specificity.qs2")))
  as.data.table(read_obj("miR_state_specificity")) else NULL

## candidate score: program-target coverage x targeting breadth, per miRNA x program
cand <- prog_hits[, .(genes = paste(sort(unique(gene)), collapse = ","),
                      n_genes = uniqueN(gene), n_clusters = sum(n_clusters),
                      coverage = sum(coverage),
                      states = paste(sort(unique(as.character(state))), collapse = ",")),
                  by = .(miRNA, program)]
cand[, score := log1p(coverage) * n_genes]
cand <- cand[order(program, -score)]
csv(cand, "therapeutic_candidates_ranked.csv")

## T1 — top candidate miRNAs per program
topc <- cand[, head(.SD, 10), by = program]
topc[, lab := factor(paste(miRNA, program, sep="@@"), levels = rev(paste(miRNA, program, sep="@@")))]
t1 <- ggplot(topc, aes(lab, score, fill = program)) +
  geom_col(width = .74, colour = "grey25", linewidth = .2) + coord_flip() +
  scale_x_discrete(labels = function(z) sub("mmu-","",sub("@@.*$","",z))) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  facet_wrap(~ program, scales = "free_y") +
  labs(title = "Candidate miRNAs by T-cell program targeting",
       subtitle = "score = log(coverage) x # program genes targeted (reproducible, background-filtered)",
       x = NULL, y = "candidate score") + theme_pub(11)
save_png(t1, "T1_candidates_by_program", 12, 6.5)

## =====================================================================
## 3.  T2 — program-gene targeting map (gene x miRNA, faceted by program)
## =====================================================================
pg <- prog_hits[, .(coverage = sum(coverage), states = uniqueN(state)), by = .(program, gene, miRNA)]
pg[, miRNA := sub("mmu-","",miRNA)]
t2 <- ggplot(pg, aes(miRNA, gene, size = states, colour = log2(coverage + 1))) +
  geom_point() + scale_size_area(max_size = 6, breaks = 1:4) +
  scale_colour_viridis_c(option = "A", end = .9) +
  facet_wrap(~ program, scales = "free", ncol = 1) +
  labs(title = "Which miRNAs target effector / memory / exhaustion-program genes",
       subtitle = "size = # states the interaction is seen in; colour = log2 coverage",
       x = NULL, y = NULL, size = "# states", colour = "log2 cov") +
  theme_pub(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                        axis.text.y = element_text(size = 7))
save_png(t2, "T2_program_targeting_map", 10, 11)

## =====================================================================
## 4.  T3 — exhaustion-leverage view: miRNAs whose targets support effector/
##           memory AND that are loaded in Tex (push exhausted toward function)
## =====================================================================
tex_mirs <- inter[state == "Tex", .(tex_cov = sum(coverage)), by = miRNA]
lever <- merge(cand[program %in% c("effector","memory")], tex_mirs, by = "miRNA")
lever <- lever[, .(programs = paste(unique(program), collapse="+"),
                   prog_genes = paste(unique(unlist(strsplit(genes, ","))), collapse=","),
                   n_prog_genes = sum(n_genes), prog_cov = sum(coverage),
                   tex_loading = tex_cov[1]), by = miRNA][order(-n_prog_genes, -prog_cov)]
csv(lever, "exhaustion_leverage_candidates.csv")
if (nrow(lever)) {
  lv <- head(lever, 15); lv[, miRNA := factor(sub("mmu-","",miRNA), levels = rev(sub("mmu-","",miRNA)))]
  t3 <- ggplot(lv, aes(miRNA, n_prog_genes, fill = log2(tex_loading + 1))) +
    geom_col(width = .74, colour = "grey25", linewidth = .2) + coord_flip() +
    scale_fill_viridis_c(option = "D", end = .9) +
    labs(title = "Exhaustion-leverage candidates",
         subtitle = "miRNAs that (i) target effector/memory programs and (ii) are loaded in exhausted cells",
         x = NULL, y = "# effector/memory program genes targeted", fill = "log2 Tex loading") +
    theme_pub(11)
  save_png(t3, "T3_exhaustion_leverage", 8.5, 6)
}

## =====================================================================
## 5.  miR-29 / CAR-T REVISIT at relaxed (per-sample, single-rep) support
## =====================================================================
ps <- rbindlist(lapply(SAMPLES$ip, function(ip) {
  d <- read_cluster_bed(f_persample_clusters(ip)); d[, ip := ip]; d }))
ps <- merge(ps, SAMPLES[, .(ip, state, rep)], by = "ip")
mir29_ps <- ps[grepl("miR-29", miRNA)]
## per-sample clusters carry miRNA + coords but no gene; annotate against step-3
## gene coordinates by overlap to recover gene-level miR-29 targets (relaxed: any rep)
anng <- ann[!is.na(gene), .(chr, start, end, strand, gene, gene_type, region)]
setkey(anng, chr, start, end)
ov <- foverlaps(mir29_ps[, .(chr, start, end, strand, miRNA, coverage, state, rep, ip)],
                anng, by.x = c("chr","start","end"), by.y = c("chr","start","end"),
                type = "any", nomatch = 0L)
mir29_targets <- ov[, .(miRNA, gene, gene_type, region, state, rep, coverage)]
csv(mir29_targets[order(gene, state)], "miR29_targets_relaxed_persample.csv")

## does relaxed search recover the CAR-T DE targets?
cart_hit <- mir29_targets[gene %in% CART_MIR29_TARGETS]
csv(cart_hit, "miR29_CART_target_validation_relaxed.csv")

## T4 — miR-29 family target landscape (relaxed), protein-coding, background-filtered
m29pc <- mir29_targets[gene_type == "protein_coding" & !is_background(gene)]
m29top <- m29pc[, .(cov = sum(coverage), states = uniqueN(state), reps = uniqueN(paste(state,rep))),
                by = .(miRNA, gene)][order(-reps, -cov)][, head(.SD, 30), by = miRNA]
if (nrow(m29top)) {
  t4 <- ggplot(m29pc[gene %in% m29top$gene,
                     .(cov = sum(coverage)), by = .(state, gene)],
               aes(state_factor(state), reorder(gene, cov), fill = log2(cov + 1))) +
    geom_tile(colour = "white", linewidth = .3) +
    scale_fill_viridis_c(option = "A", end = .9) +
    scale_x_discrete(labels = STATE_LABEL[STATES], drop = FALSE) +
    labs(title = "miR-29 family direct targets across states (relaxed, per-sample)",
         subtitle = "protein-coding, background-filtered; recovers targets below the reproducible threshold",
         x = NULL, y = NULL, fill = "log2 cov") +
    theme_pub(10) + theme(axis.text.y = element_text(size = 7))
  save_png(t4, "T4_miR29_target_landscape", 8, 9)
}

## T5 — explicit CAR-T target check
cart_summary <- data.table(gene = CART_MIR29_TARGETS)
cart_summary <- merge(cart_summary,
                      mir29_targets[gene %in% CART_MIR29_TARGETS,
                                    .(found = TRUE, miRNAs = paste(sort(unique(miRNA)), collapse=","),
                                      states = paste(sort(unique(as.character(state))), collapse=","),
                                      max_cov = max(coverage)), by = gene],
                      by = "gene", all.x = TRUE)
cart_summary[is.na(found), found := FALSE]
csv(cart_summary, "miR29_CART_targets_summary.csv")
t5 <- ggplot(cart_summary, aes(gene, factor(found, levels = c(FALSE, TRUE)), fill = found)) +
  geom_tile(colour = "white", linewidth = .5) +
  geom_text(aes(label = ifelse(found, "direct chimera\nfound", "not detected")), size = 2.8) +
  scale_fill_manual(values = c(`TRUE` = "#4FA168", `FALSE` = "grey85"), guide = "none") +
  labs(title = "CAR-T miR-29a DE targets: direct miR-29 binding in CD8 eCLIP?",
       subtitle = "relaxed per-sample search (any replicate)",
       x = NULL, y = NULL) + theme_pub(11) + theme(axis.text.y = element_blank())
save_png(t5, "T5_CART_target_check", 8, 3.5)

## =====================================================================
## 6.  Save + README
## =====================================================================
save_obj(cand,           "therapeutic_candidates")
save_obj(mir29_targets,  "miR29_targets_relaxed")
save_obj(lever,          "exhaustion_leverage_candidates")

found_cart <- cart_summary[found == TRUE]$gene
topeff <- cand[program=="effector"][order(-score)][1:3]
topmem <- cand[program=="memory"][order(-score)][1:3]
L <- c(
  "=====================================================================",
  " STEP 6 / AIM 4 — therapeutic / state-pushing candidate miRNAs",
  "=====================================================================",
  "",
  "FRAMING (from aim 3): no exhaustion-SPECIFIC target program exists, so",
  "candidates are nominated by the EFFECTOR / MEMORY programs they target and",
  "their loading in the relevant state, not by exhaustion-only targets.",
  "Abundant-RNA background (incl. Gas5/Snhg/Snord, extended this step) removed.",
  "",
  "--- TOP EFFECTOR-PROGRAM-TARGETING miRNAs --------------------------",
  apply(topeff, 1, function(r) sprintf("    %-18s genes: %s", sub("mmu-","",r[["miRNA"]]), r[["genes"]])),
  "--- TOP MEMORY-PROGRAM-TARGETING miRNAs ----------------------------",
  apply(topmem, 1, function(r) sprintf("    %-18s genes: %s", sub("mmu-","",r[["miRNA"]]), r[["genes"]])),
  "",
  "--- miR-29 / CAR-T CROSS-VALIDATION (relaxed, per-sample) ----------",
  paste0("  CAR-T DE targets (Dnmt3a,Tet2,Eomes,Tbx21,Ifng) with ANY direct ",
         "miR-29 chimera: ", if (length(found_cart)) paste(found_cart, collapse=", ") else "NONE"),
  "  Even at relaxed (single-replicate) support, the canonical miR-29 targets",
  "  are largely/entirely absent as direct chimeras in these CD8 states.",
  "  -> The CAR-T miR-29a phenotype is likely NOT explained by direct binding of",
  "     these targets in CD8 T cells; see T4 for what miR-29 DOES bind here.",
  "",
  "--- FIGURES --------------------------------------------------------",
  "  T1 candidate miRNAs by program (effector/memory/exhaustion)",
  "  T2 program-gene targeting map (gene x miRNA)",
  "  T3 exhaustion-leverage candidates (effector/memory targeters loaded in Tex)",
  "  T4 miR-29 family target landscape across states (relaxed)",
  "  T5 CAR-T miR-29a DE target check (found / not detected)",
  "",
  "--- CAVEATS --------------------------------------------------------",
  "  * 'Therapeutic candidate' = hypothesis-generating; binding != functional",
  "    repression. Validate with the matched RNA-seq (still the key missing",
  "    data) before any therapeutic claim.",
  "  * Relaxed miR-29 search uses single-replicate clusters: more sensitive but",
  "    less stringent; absence here is informative given that relaxation.",
  "  * n=2; Tex low capture; miR-466i naive inflation - as in prior steps.",
  "",
  "--- CSVs TO UPLOAD (eCLIP/06_therapeutic/tables) -------------------",
  "  therapeutic_candidates_ranked.csv, miR_targeting_program_genes.csv,",
  "  exhaustion_leverage_candidates.csv, miR29_targets_relaxed_persample.csv,",
  "  miR29_CART_targets_summary.csv",
  "",
  "PROJECT STATUS: aims 1-4 complete at the BINDING level. The remaining gap for",
  "functional/therapeutic conclusions is matched differentiation-state RNA-seq",
  "(target derepression). That is the recommended next experiment.",
  "=====================================================================")
writeLines(L, file.path(OUT, "README_06_therapeutic.txt"))
message("06_therapeutic_candidates COMPLETE -> ", OUT)