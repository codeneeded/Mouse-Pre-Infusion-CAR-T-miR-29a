## =====================================================================
## 03_annotate_targets.R  —  STEP 3 / AIM 2 prep: cluster -> gene + region
##
## The reproducible chimeric clusters carry a miRNA + coordinates but NO gene.
## Here we assign each cluster a target GENE, REGION and BIOTYPE by genomic
## overlap against GENCODE vM25 (mm10, UCSC chr names). This also RESOLVES the
## ~30% "Other" feature bucket flagged in step-1 QC into real categories
## (sense ncRNA / antisense / intergenic).
##
## Region priority (sense, high->low): 3'UTR > 5'UTR > CDS > ncRNA_exon > intron
##   no sense gene -> antisense (opposite-strand gene) -> intergenic
##
## Output (flat): eCLIP/03_annotation/ A1..A4 PNGs + tables/ + README.
## Object (qs2): clusters_annotated -> saved_R_data/.
## =====================================================================

## ---- libraries (plain loads) -------------------------------------------
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(qs2)

## =====================================================================
## CONFIG (self-contained — edit paths)
## =====================================================================
DATA_DIR <- "/media/akshay-iyer/Elements/data_from_hpc/QN-0000916_miR-eCLIP_No-Gel/files"
ROOT_OUT <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/eCLIP"
OBJ_DIR  <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/saved_R_data"
## GENCODE vM25 (mm10) GTF. Download once:
##   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
##   gunzip gencode.vM25.annotation.gtf.gz
## On BioC >= 3.19 also: BiocManager::install("txdbmaker")
GTF_PATH <- "/home/akshay-iyer/Documents/refs/gencode.vM25.annotation.gtf"   # EDIT

STATES      <- c("Tn","Teff","Tmem","Tex")
STATE_LABEL <- c(Tn="Naive", Teff="Effector", Tmem="Memory", Tex="Exhausted")
state_factor <- function(x) factor(x, levels = STATES)
STATE_COLS  <- c(Tn="#3B6FB6", Teff="#E08A3C", Tmem="#4FA168", Tex="#C0413B")
REGION_LEVELS <- c("3' UTR","5' UTR","CDS","ncRNA_exon","intron","antisense","intergenic")
REGION_COLS <- c("3' UTR"="#1F77B4","5' UTR"="#17BECF","CDS"="#FF7F0E",
                 "ncRNA_exon"="#9467BD","intron"="#2CA02C","antisense"="#E377C2","intergenic"="#9AA0A6")
BIOTYPE_COLS <- c(protein_coding="#FF7F0E", lncRNA="#9467BD", miRNA="#D62728",
                  pseudogene="#8C564B", other_ncRNA="#17BECF", other="#BCBD22",
                  unannotated="#9AA0A6")

f_repro_clusters <- function(st) file.path(DATA_DIR, paste0(st, ".mir_targets_reproducible_clusters.bed"))
read_cluster_bed <- function(path) {
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

OUT <- file.path(ROOT_OUT, "03_annotation"); TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)
save_png <- function(p, name, w = 8, h = 5.5)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
csv <- function(x, f) fwrite(x, file.path(TAB, f))
save_obj <- function(x, n) { qs2::qs_save(x, file.path(OBJ_DIR, paste0(n, ".qs2"))); message("saved ", n) }

## focus sets (preview aims 2/4 + CAR-T link)
MASTER_TF <- c("Tcf7","Lef1","Myb","Bach2","Tbx21","Eomes","Prdm1","Tox",
               "Nr4a1","Nr4a2","Nr4a3","Irf4","Batf","Zeb2","Id2","Id3")
EXHAUSTION_GENES <- c("Pdcd1","Havcr2","Lag3","Tigit","Ctla4","Entpd1")
CART_MIR29_TARGETS <- c("Dnmt3a","Tet2","Eomes","Tbx21","Ifng")

## =====================================================================
## 1.  Load reproducible clusters + build GENCODE annotation
## =====================================================================
if (!file.exists(GTF_PATH)) stop(
  "GENCODE vM25 GTF not found at GTF_PATH.\n",
  "Download it once:\n",
  "  mkdir -p ", dirname(GTF_PATH), "\n",
  "  cd ", dirname(GTF_PATH), "\n",
  "  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz\n",
  "  gunzip gencode.vM25.annotation.gtf.gz")
clusters_repro <- rbindlist(lapply(STATES, function(st) {
  dt <- read_cluster_bed(f_repro_clusters(st)); dt[, state := st]; dt }))
clusters_repro[, state := state_factor(state)]

message("building TxDb from GENCODE vM25 (one-time, ~1-2 min)...")
txdb <- suppressWarnings(
  if (requireNamespace("txdbmaker", quietly = TRUE))
    txdbmaker::makeTxDbFromGFF(GTF_PATH, format = "gtf")
  else GenomicFeatures::makeTxDbFromGFF(GTF_PATH, format = "gtf")
)
gtf_genes <- rtracklayer::import(GTF_PATH, feature.type = "gene")
strip <- function(x) sub("\\..*$", "", x)
g2name <- setNames(gtf_genes$gene_name, strip(gtf_genes$gene_id))
g2type <- setNames(gtf_genes$gene_type, strip(gtf_genes$gene_id))

genes_gr <- genes(txdb)
utr3 <- unlist(threeUTRsByTranscript(txdb)); utr5 <- unlist(fiveUTRsByTranscript(txdb))
cds  <- unlist(cdsBy(txdb, by = "tx"));      exn  <- exons(txdb)

## seqlevel sanity (BEDs are UCSC chr; GENCODE is UCSC chr — should match)
common <- intersect(seqlevels(genes_gr), unique(clusters_repro$chr))
message(length(common), " shared seqlevels with annotation")

## =====================================================================
## 2.  Annotate (vectorized, strand-aware)
## =====================================================================
cl <- GRanges(clusters_repro$chr, IRanges(clusters_repro$start, clusters_repro$end),
              strand = clusters_repro$strand)
n <- length(cl)
region <- rep("intron", n); gene_id <- rep(NA_character_, n); strand_rel <- rep(NA_character_, n)

## (a) sense gene: pick best overlap (protein_coding first, then widest)
ov <- findOverlaps(cl, genes_gr, ignore.strand = FALSE)
if (length(ov)) {
  d <- data.table(q = queryHits(ov), s = subjectHits(ov),
                  w = width(pintersect(cl[queryHits(ov)], genes_gr[subjectHits(ov)])),
                  gid = strip(names(genes_gr))[subjectHits(ov)])
  d[, pc := as.integer(g2type[gid] == "protein_coding")]
  setorder(d, q, -pc, -w)
  best <- d[, .SD[1L], by = q]
  gene_id[best$q] <- best$gid; strand_rel[best$q] <- "sense"
}
has_gene <- !is.na(gene_id)
region[!has_gene] <- "intergenic"

## (b) sub-genic region, priority low->high (only where a sense gene exists)
setreg <- function(reg, feats, label) {
  idx <- unique(queryHits(findOverlaps(cl, feats, ignore.strand = FALSE)))
  idx <- idx[has_gene[idx]]; reg[idx] <- label; reg }
region <- setreg(region, exn,  "ncRNA_exon")   # generic exon; coding parts overwritten next
region <- setreg(region, cds,  "CDS")
region <- setreg(region, utr5, "5' UTR")
region <- setreg(region, utr3, "3' UTR")

## (c) antisense for the still-intergenic (opposite-strand gene overlap)
ai <- which(region == "intergenic")
if (length(ai)) {
  ova <- findOverlaps(cl[ai], genes_gr, ignore.strand = TRUE)
  if (length(ova)) {
    da <- data.table(qi = ai[queryHits(ova)], s = subjectHits(ova))
    da <- da[as.character(strand(cl))[qi] != as.character(strand(genes_gr))[s]]
    if (nrow(da)) {
      fa <- da[, .SD[1L], by = qi]
      region[fa$qi] <- "antisense"; gene_id[fa$qi] <- strip(names(genes_gr))[fa$s]
      strand_rel[fa$qi] <- "antisense"
    }
  }
}

clusters_annotated <- copy(clusters_repro)
clusters_annotated[, `:=`(gene_id = gene_id, gene = g2name[gene_id],
                          gene_type = g2type[gene_id],
                          region = factor(region, levels = REGION_LEVELS),
                          strand_rel = strand_rel)]
clusters_annotated[is.na(gene) & !is.na(gene_id), gene := gene_id]
collapse_biotype <- function(t) fcase(
  is.na(t), "unannotated", t == "protein_coding", "protein_coding", t == "miRNA", "miRNA",
  t %in% c("lncRNA","antisense","lincRNA","processed_transcript","sense_intronic",
           "sense_overlapping","bidirectional_promoter_lncRNA","macro_lncRNA",
           "3prime_overlapping_ncRNA","TEC"), "lncRNA",
  grepl("pseudogene", t), "pseudogene",
  t %in% c("snoRNA","snRNA","misc_RNA","scaRNA","rRNA","Mt_rRNA","Mt_tRNA","ribozyme",
           "sRNA","scRNA","vault_RNA","miRNA"), "other_ncRNA", default = "other")
clusters_annotated[, biotype := factor(collapse_biotype(gene_type),
                                       levels = names(BIOTYPE_COLS))]
csv(clusters_annotated, "annotated_reproducible_clusters.csv")
save_obj(clusters_annotated, "clusters_annotated")

## =====================================================================
## 3.  Summaries
## =====================================================================
region_comp <- clusters_annotated[, .N, by = .(state, region)][, frac := N/sum(N), by = state]
csv(region_comp, "region_composition_per_state.csv")
biotype_comp <- clusters_annotated[, .N, by = .(state, biotype)][, frac := N/sum(N), by = state]
csv(biotype_comp, "biotype_composition_per_state.csv")

genes_per_state <- clusters_annotated[!is.na(gene) & biotype == "protein_coding",
                                      .(unique_target_genes = uniqueN(gene)), by = state]
csv(genes_per_state, "unique_proteincoding_targets_per_state.csv")

## key-gene targeting teaser (previews aims 2/4 + CAR-T link)
key <- unique(c(MASTER_TF, EXHAUSTION_GENES, CART_MIR29_TARGETS))
key_hits <- clusters_annotated[gene %in% key,
                               .(state, miRNA, gene, region, coverage)][order(gene, state, -coverage)]
csv(key_hits, "key_gene_targeting_preview.csv")

## annotation-outcome sanity (resolves the old "Other" bucket)
outcome <- clusters_annotated[, .N, by = .(state,
                                           outcome = fcase(region %in% c("3' UTR","5' UTR","CDS"), "coding gene",
                                                           region == "ncRNA_exon", "ncRNA exon",
                                                           region == "intron", "intron",
                                                           region == "antisense", "antisense",
                                                           default = "intergenic"))][, frac := N/sum(N), by = state]
csv(outcome, "annotation_outcome_per_state.csv")

## =====================================================================
## 4.  FIGURES
## =====================================================================
## A1 region composition per state (resolves "Other")
a1 <- ggplot(region_comp, aes(state, frac, fill = region)) +
  geom_col(width = .72, colour = "grey30", linewidth = .15) +
  scale_fill_manual(values = REGION_COLS) +
  scale_x_discrete(labels = STATE_LABEL[STATES]) + scale_y_continuous(labels = percent_format()) +
  labs(title = "Target region composition per state",
       subtitle = "the step-1 'Other' bucket resolved: antisense + intergenic + ncRNA now explicit",
       x = NULL, y = "fraction of clusters", fill = "region") + theme_pub()
save_png(a1, "A1_region_composition", 8.5, 5.5)

## A2 target biotype per state
b2 <- ggplot(biotype_comp, aes(state, frac, fill = biotype)) +
  geom_col(width = .72, colour = "grey30", linewidth = .15) +
  scale_fill_manual(values = BIOTYPE_COLS) +
  scale_x_discrete(labels = STATE_LABEL[STATES]) + scale_y_continuous(labels = percent_format()) +
  labs(title = "Target gene biotype per state", x = NULL, y = "fraction of clusters",
       fill = "biotype") + theme_pub()
save_png(b2, "A2_target_biotype", 8.5, 5.5)

## A3 unique protein-coding target genes per state
a3 <- ggplot(genes_per_state, aes(state, unique_target_genes, fill = state)) +
  geom_col(width = .68, colour = "grey25", linewidth = .2) +
  geom_text(aes(label = comma(unique_target_genes)), vjust = -.3, size = 3.2) +
  scale_fill_manual(values = STATE_COLS, guide = "none") +
  scale_x_discrete(labels = STATE_LABEL[STATES]) +
  labs(title = "Distinct protein-coding target genes per state",
       x = NULL, y = "unique genes") + theme_pub() + expand_limits(y = 0)
save_png(a3, "A3_unique_target_genes", 7.5, 5)

## A4 key-gene targeting preview (immune/CAR-T genes x state, # miRNAs binding)
if (nrow(key_hits)) {
  kk <- key_hits[, .(n_miRNAs = uniqueN(miRNA), total_cov = sum(coverage)), by = .(gene, state)]
  kk[, state := state_factor(state)]
  a4 <- ggplot(kk, aes(state, gene, size = n_miRNAs, colour = log2(total_cov + 1))) +
    geom_point() + scale_size_area(max_size = 7) +
    scale_colour_viridis_c(option = "A", end = .9) +
    scale_x_discrete(labels = STATE_LABEL[STATES]) +
    labs(title = "Direct targeting of key immune / CAR-T genes (preview)",
         subtitle = "dot size = # miRNAs binding the gene in that state",
         x = NULL, y = NULL, size = "# miRNAs", colour = "log2 coverage") +
    theme_pub(11) + theme(axis.text.y = element_text(size = 8))
  save_png(a4, "A4_key_gene_targeting_preview", 8, max(5, 0.32 * uniqueN(kk$gene)))
}

## =====================================================================
## 5.  README explainer
## =====================================================================
asn <- clusters_annotated[, .(
  pct_sense = round(100*mean(strand_rel == "sense", na.rm = TRUE), 1),
  pct_antisense = round(100*mean(region == "antisense"), 1),
  pct_intergenic = round(100*mean(region == "intergenic"), 1)), by = state]
rc <- dcast(region_comp, region ~ state, value.var = "frac")
keygenes_found <- sort(unique(key_hits$gene))
mir29_found <- key_hits[gene %in% CART_MIR29_TARGETS & grepl("miR-29", miRNA)]
L <- c(
  "=====================================================================",
  " STEP 3 / AIM 2 prep — annotate reproducible clusters to genes",
  "=====================================================================",
  "",
  "WHAT THIS DOES: assigns every reproducible chimeric cluster a target GENE,",
  "REGION and BIOTYPE by overlap with GENCODE vM25 (mm10). Strand-aware; the",
  "~30% 'Other' bucket from step-1 QC is now resolved into sense ncRNA /",
  "antisense / intergenic.",
  "",
  "REGION PRIORITY (sense): 3'UTR > 5'UTR > CDS > ncRNA_exon > intron;",
  "  no sense gene -> antisense (opposite-strand gene) -> intergenic.",
  "",
  "--- ANNOTATION OUTCOME PER STATE (% sense / antisense / intergenic) -",
  apply(asn, 1, function(r) sprintf("    %-5s  sense %s%% | antisense %s%% | intergenic %s%%",
                                    r[["state"]], r[["pct_sense"]], r[["pct_antisense"]], r[["pct_intergenic"]])),
  "",
  "--- DISTINCT PROTEIN-CODING TARGET GENES PER STATE -----------------",
  apply(genes_per_state, 1, function(r) sprintf("    %-5s  %s genes", r[["state"]], r[["unique_target_genes"]])),
  "",
  "--- KEY GENES DIRECTLY TARGETED (preview of aims 2/4) --------------",
  paste0("  Found among reproducible targets: ",
         if (length(keygenes_found)) paste(keygenes_found, collapse = ", ") else "(none)"),
  paste0("  miR-29 family -> CAR-T DE targets bound here: ",
         if (nrow(mir29_found)) paste(sprintf("%s(%s)", mir29_found$gene, mir29_found$miRNA), collapse=", ")
         else "none directly (check step 4/6 with relaxed support)"),
  "",
  "--- FIGURES --------------------------------------------------------",
  "  A1 region composition per state (Other resolved)",
  "  A2 target gene biotype per state",
  "  A3 distinct protein-coding target genes per state",
  "  A4 key immune/CAR-T gene targeting preview (dot = #miRNAs binding)",
  "",
  "--- CAVEATS --------------------------------------------------------",
  "  * Clusters overlapping multiple genes keep the best (protein_coding, then",
  "    widest overlap); multi-gene loci are simplified.",
  "  * intergenic/antisense clusters are real eCLIP signal (AGO can bind ncRNA /",
  "    nascent transcripts); they are retained, not discarded.",
  "  * Annotation uses GENCODE vM25 (mm10) to MATCH the eCLIP processing genome,",
  "    NOT the CAR-T GRCm39 reference.",
  "",
  "--- CSVs TO UPLOAD (eCLIP/03_annotation/tables) --------------------",
  "  region_composition_per_state.csv, biotype_composition_per_state.csv,",
  "  unique_proteincoding_targets_per_state.csv, key_gene_targeting_preview.csv,",
  "  annotation_outcome_per_state.csv",
  "",
  "NEXT (step 4, aim 2): per-state miR->gene interaction networks, hub miRNAs,",
  "  most-targeted genes, and target switching across states.",
  "=====================================================================")
writeLines(L, file.path(OUT, "README_03_annotation.txt"))
message("03_annotate_targets COMPLETE -> ", OUT)
