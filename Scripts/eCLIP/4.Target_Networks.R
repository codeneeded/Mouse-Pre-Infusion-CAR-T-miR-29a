## =====================================================================
## 04_target_networks.R  —  STEP 4 / AIM 2: what each miRNA targets per state
##
## Builds per-state miR -> gene interaction tables from the annotated
## reproducible clusters (step 3), then: most-targeted genes, hub miRNAs
## (targeting breadth), bipartite networks (graphml for Cytoscape), and a
## first pass at TARGET SWITCHING (same miRNA, different gene set per state).
##
## Restricted to PROTEIN-CODING targets (per step-3 interpretation; ncRNA/
## structural-RNA hits set aside for the biological/therapeutic story).
##
## Input object: clusters_annotated (step 3).  Falls back to re-annotation
## note if absent.  Output (flat): eCLIP/04_networks/ N1..N5 PNGs + tables/.
## =====================================================================

## ---- libraries (plain loads) -------------------------------------------
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(ggrepel)
library(igraph)
library(qs2)

## =====================================================================
## CONFIG (self-contained)
## =====================================================================
ROOT_OUT <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/eCLIP"
OBJ_DIR  <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a/saved_R_data"

STATES      <- c("Tn","Teff","Tmem","Tex")
STATE_LABEL <- c(Tn="Naive", Teff="Effector", Tmem="Memory", Tex="Exhausted")
state_factor <- function(x) factor(x, levels = STATES)
STATE_COLS  <- c(Tn="#3B6FB6", Teff="#E08A3C", Tmem="#4FA168", Tex="#C0413B")

theme_pub <- function(base_size = 13) theme_minimal(base_size) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = .25, colour = "grey88"),
        axis.line = element_line(linewidth = .4, colour = "grey20"),
        axis.ticks = element_line(linewidth = .3, colour = "grey40"),
        plot.title = element_text(face = "bold", size = rel(1.05)),
        plot.subtitle = element_text(colour = "grey35", size = rel(.85)),
        strip.text = element_text(face = "bold"), legend.key.size = unit(.8,"lines"),
        plot.margin = margin(10,14,10,10))

OUT <- file.path(ROOT_OUT, "04_networks"); TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
save_png <- function(p, name, w = 8, h = 5.5)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
csv <- function(x, f) fwrite(x, file.path(TAB, f))
save_obj <- function(x, n) { qs2::qs_save(x, file.path(OBJ_DIR, paste0(n, ".qs2"))); message("saved ", n) }
read_obj <- function(n) qs2::qs_read(file.path(OBJ_DIR, paste0(n, ".qs2")))

MASTER_TF <- c("Tcf7","Lef1","Myb","Bach2","Tbx21","Eomes","Prdm1","Tox",
               "Nr4a1","Nr4a2","Nr4a3","Irf4","Batf","Zeb2","Id2","Id3")

## =====================================================================
## 1.  Load annotated clusters; keep protein-coding, sense, gene-named
## =====================================================================
stopifnot("run step 3 first (clusters_annotated.qs2 missing)" =
            file.exists(file.path(OBJ_DIR, "clusters_annotated.qs2")))
ann <- read_obj("clusters_annotated")
ann <- as.data.table(ann)
ann[, state := state_factor(as.character(state))]
pc <- ann[biotype == "protein_coding" & !is.na(gene) & gene != ""]
message(nrow(pc), " protein-coding cluster assignments across states")

## =====================================================================
## 2.  Per-state miR -> gene edges (collapse clusters of same miR on same gene)
## =====================================================================
edges <- pc[, .(n_clusters = .N, coverage = sum(coverage),
                regions = paste(sort(unique(as.character(region))), collapse = ",")),
            by = .(state, miRNA, gene)]
csv(edges, "interactions_by_state.csv")

## =====================================================================
## 3.  N1 — most-targeted genes per state (by # distinct miRNAs + coverage)
## =====================================================================
top_genes <- edges[, .(n_miRNAs = uniqueN(miRNA), total_clusters = sum(n_clusters),
                       total_cov = sum(coverage)), by = .(state, gene)][
                         order(state, -n_miRNAs, -total_cov)]
csv(top_genes, "top_targeted_genes_per_state.csv")
tg <- top_genes[, head(.SD, 15), by = state]
tg[, lab := factor(paste(gene, state, sep="@@"), levels = rev(paste(gene, state, sep="@@")))]
n1 <- ggplot(tg, aes(lab, n_miRNAs, fill = state)) +
  geom_col(width = .74, colour = "grey25", linewidth = .2) + coord_flip() +
  scale_x_discrete(labels = function(z) sub("@@.*$", "", z)) +
  scale_fill_manual(values = STATE_COLS, guide = "none") +
  facet_wrap(~ state, scales = "free_y", labeller = labeller(state = STATE_LABEL)) +
  labs(title = "Most-targeted genes per state (miRNA convergence)",
       subtitle = "# distinct miRNAs binding each gene's transcript",
       x = NULL, y = "# distinct miRNAs") + theme_pub(11)
save_png(n1, "N1_top_targeted_genes", 11, 8)

## =====================================================================
## 4.  N2 — hub miRNAs per state (targeting breadth)
## =====================================================================
hub <- edges[, .(n_targets = uniqueN(gene), n_clusters = sum(n_clusters),
                 coverage = sum(coverage)), by = .(state, miRNA)][
                   order(state, -n_targets)]
csv(hub, "hub_miRNAs_per_state.csv")
hb <- hub[, head(.SD, 15), by = state]
hb[, lab := factor(paste(miRNA, state, sep="@@"), levels = rev(paste(miRNA, state, sep="@@")))]
n2 <- ggplot(hb, aes(lab, n_targets, fill = state)) +
  geom_col(width = .74, colour = "grey25", linewidth = .2) + coord_flip() +
  scale_x_discrete(labels = function(z) sub("mmu-","",sub("@@.*$", "", z))) +
  scale_fill_manual(values = STATE_COLS, guide = "none") +
  facet_wrap(~ state, scales = "free_y", labeller = labeller(state = STATE_LABEL)) +
  labs(title = "Hub miRNAs per state (targeting breadth)",
       subtitle = "# distinct protein-coding genes targeted",
       x = NULL, y = "# target genes") + theme_pub(11)
save_png(n2, "N2_hub_miRNAs", 11, 8)

## =====================================================================
## 5.  N3 — per-state networks: curated PNG (all states) + interactive HTML + graphml
## =====================================================================
MIN_EDGE_WEIGHT <- 1   # raise (e.g. 2-3) if a state's interactive HTML is too dense
have_vis <- requireNamespace("visNetwork", quietly = TRUE) &&
  requireNamespace("htmlwidgets", quietly = TRUE)
if (!have_vis) message("NOTE: install.packages(c('visNetwork','htmlwidgets','pandoc')) for interactive HTML networks")

## self-contained HTML needs pandoc; provision a user-local pandoc if none is found
ensure_pandoc <- function() {
  if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) return(TRUE)
  if (requireNamespace("pandoc", quietly = TRUE)) {
    if (!isTRUE(try(pandoc::pandoc_is_installed(), silent = TRUE))) try(pandoc::pandoc_install(), silent = TRUE)
    try(pandoc::pandoc_activate(), silent = TRUE)
    return(isTRUE(try(rmarkdown::pandoc_available(), silent = TRUE)))
  }
  FALSE
}
pandoc_ok <- have_vis && ensure_pandoc()
if (have_vis && !pandoc_ok)
  message("NOTE: pandoc not found -> cannot write SELF-CONTAINED HTML. Install once with\n",
          "      install.packages('pandoc'); pandoc::pandoc_install()\n",
          "      Skipping HTML this run (PNG + graphml still written).")

## curated, legible static PNG: top hub miRNAs x (master TFs + top targets)
curated_network_png <- function(st) {
  ep <- edges[state == st & n_clusters >= MIN_EDGE_WEIGHT]
  keep_mir  <- hub[state == st][order(-n_targets)][1:8]$miRNA
  key_genes <- union(MASTER_TF, top_genes[state == st][1:40]$gene)
  ep <- ep[miRNA %in% keep_mir & gene %in% key_genes]
  if (nrow(ep) == 0) { message("no curated edges for ", st); return(invisible(NULL)) }
  gg <- graph_from_data_frame(ep[, .(from = miRNA, to = gene, weight = n_clusters)], directed = TRUE)
  V(gg)$kind <- ifelse(V(gg)$name %in% ep$miRNA, "miRNA", "gene")
  set.seed(1); lay <- layout_with_fr(gg)
  dfv <- data.table(name = V(gg)$name, kind = V(gg)$kind, x = lay[,1], y = lay[,2])
  de  <- as.data.table(igraph::as_data_frame(gg, "edges"))
  de  <- merge(de, dfv[, .(from = name, x1 = x, y1 = y)], by = "from")
  de  <- merge(de, dfv[, .(to   = name, x2 = x, y2 = y)], by = "to")
  p <- ggplot() +
    geom_segment(data = de, aes(x1, y1, xend = x2, yend = y2, linewidth = weight),
                 colour = "grey75", alpha = .6) +
    geom_point(data = dfv, aes(x, y, colour = kind, size = kind)) +
    geom_text_repel(data = dfv, aes(x, y, label = sub("mmu-", "", name)), size = 2.6, max.overlaps = 40) +
    scale_colour_manual(values = c(miRNA = "#C0413B", gene = "#3B6FB6")) +
    scale_size_manual(values = c(miRNA = 4, gene = 2.5), guide = "none") +
    scale_linewidth(range = c(.3, 2), guide = "none") +
    labs(title = paste0(STATE_LABEL[st], " miR->target network (hub miRNAs x key genes)"),
         colour = NULL) +
    theme_void() + theme(plot.title = element_text(face = "bold"))
  save_png(p, paste0("N3_network_", st), 9, 8)
}

## full interactive HTML network (draggable / zoomable / searchable)
interactive_network_html <- function(st) {
  if (!have_vis || !pandoc_ok) return(invisible(NULL))
  e <- edges[state == st & n_clusters >= MIN_EDGE_WEIGHT]
  if (nrow(e) == 0) return(invisible(NULL))
  mirs  <- unique(e$miRNA); genes <- setdiff(unique(e$gene), mirs)
  deg   <- rbind(e[, .(id = miRNA)], e[, .(id = gene)])[, .N, by = id]
  nodes <- data.table(id = c(mirs, genes))
  nodes[, group := ifelse(id %in% mirs, "miRNA", "gene")]
  nodes <- merge(nodes, deg, by = "id", all.x = TRUE)
  nodes[, `:=`(label = sub("mmu-", "", id), value = N,
               title = paste0("<b>", id, "</b><br>", group, "<br>connections: ", N))]
  edg <- e[, .(from = miRNA, to = gene, value = n_clusters,
               title = paste0(sub("mmu-","",miRNA), " -> ", gene,
                              "<br>clusters: ", n_clusters, "<br>coverage: ", round(coverage,1)))]
  vn <- visNetwork::visNetwork(nodes, edg,
                               main = paste0(STATE_LABEL[st], " miRNA->target network (full, reproducible)"),
                               width = "100%", height = "800px")
  vn <- visNetwork::visGroups(vn, groupname = "miRNA", color = "#C0413B", shape = "dot")
  vn <- visNetwork::visGroups(vn, groupname = "gene",  color = "#3B6FB6", shape = "square")
  vn <- visNetwork::visEdges(vn, arrows = "to", color = list(opacity = 0.4), smooth = FALSE)
  vn <- visNetwork::visOptions(vn, highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                               nodesIdSelection = TRUE)
  vn <- visNetwork::visPhysics(vn, stabilization = list(iterations = 150),
                               barnesHut = list(gravitationalConstant = -3000))
  vn <- visNetwork::visLegend(vn)
  f <- file.path(OUT, paste0("network_", st, ".html"))
  htmlwidgets::saveWidget(vn, f, selfcontained = TRUE)   # single self-contained file (no _files/ dir)
  message("interactive: ", f)
}

for (st in STATES) {
  ## graphml (Cytoscape) — full network
  e <- edges[state == st, .(from = miRNA, to = gene, weight = n_clusters, coverage)]
  g <- graph_from_data_frame(e, directed = TRUE)
  V(g)$type <- ifelse(V(g)$name %in% e$from, "miRNA", "gene")
  write_graph(g, file.path(OUT, paste0("network_", st, ".graphml")), format = "graphml")
  curated_network_png(st)        # legible static PNG, all 4 states
  interactive_network_html(st)   # full interactive HTML, all 4 states
}

## =====================================================================
## 6.  N4 — master-TF targeting map across states
## =====================================================================
tf <- edges[gene %in% MASTER_TF]
csv(tf[order(gene, state, -coverage)], "masterTF_targeting.csv")
if (nrow(tf)) {
  tfm <- tf[, .(n_miRNAs = uniqueN(miRNA), coverage = sum(coverage)), by = .(gene, state)]
  tfm[, state := state_factor(state)]
  n4 <- ggplot(tfm, aes(state, gene, size = n_miRNAs, colour = log2(coverage + 1))) +
    geom_point() + scale_size_area(max_size = 8) +
    scale_colour_viridis_c(option = "A", end = .9) +
    scale_x_discrete(labels = STATE_LABEL[STATES]) +
    labs(title = "Direct targeting of master TFs across states",
         subtitle = "dot size = # miRNAs binding; colour = log2 coverage",
         x = NULL, y = NULL, size = "# miRNAs", colour = "log2 cov") +
    theme_pub(11) + theme(axis.text.y = element_text(size = 8))
  save_png(n4, "N4_masterTF_targeting", 7.5, max(5, .35*uniqueN(tfm$gene)))
}

## =====================================================================
## 7.  N5 — target switching (same miRNA, different gene set across states)
## =====================================================================
gsets <- pc[, .(genes = list(sort(unique(gene)))), by = .(miRNA, state)]
mir_multi <- gsets[, .N, by = miRNA][N >= 2]$miRNA
switch_tab <- rbindlist(lapply(mir_multi, function(m) {
  sl <- gsets[miRNA == m]; cb <- combn(as.character(sl$state), 2)
  rbindlist(lapply(seq_len(ncol(cb)), function(k) {
    ga <- sl[state == cb[1,k]]$genes[[1]]; gb <- sl[state == cb[2,k]]$genes[[1]]
    inter <- length(intersect(ga, gb)); uni <- length(union(ga, gb))
    data.table(miRNA = m, stateA = cb[1,k], stateB = cb[2,k],
               nA = length(ga), nB = length(gb), shared = inter,
               jaccard = ifelse(uni > 0, inter/uni, NA_real_))
  }))
}))
switch_tab <- switch_tab[order(jaccard)]
csv(switch_tab, "target_switching_pairwise.csv")
## per-miRNA mean jaccard across state pairs (low = strong switcher), require breadth
sw_sum <- switch_tab[, .(mean_jaccard = mean(jaccard, na.rm = TRUE),
                         max_targets = max(c(nA, nB))), by = miRNA][
                           max_targets >= 5][order(mean_jaccard)]
csv(sw_sum, "target_switching_summary.csv")
sws <- head(sw_sum, 20); sws[, miRNA := factor(sub("mmu-","",miRNA), levels = rev(sub("mmu-","",miRNA)))]
n5 <- ggplot(sws, aes(miRNA, 1 - mean_jaccard, fill = max_targets)) +
  geom_col(width = .74, colour = "grey25", linewidth = .2) + coord_flip() +
  scale_fill_viridis_c(option = "D", end = .9) +
  labs(title = "Strongest target-switching miRNAs across states",
       subtitle = "1 - mean Jaccard of target sets between state pairs (higher = more switching)",
       x = NULL, y = "target-set divergence", fill = "max # targets") + theme_pub(11)
save_png(n5, "N5_target_switching", 8.5, 6.5)

## =====================================================================
## 8.  Save objects
## =====================================================================
save_obj(edges,       "interactions_by_state")
save_obj(top_genes,   "top_targeted_genes")
save_obj(hub,         "hub_miRNAs")
save_obj(switch_tab,  "target_switching")

## =====================================================================
## 9.  README explainer
## =====================================================================
ntargets <- edges[, .(genes = uniqueN(gene), miRNAs = uniqueN(miRNA), edges = .N), by = state]
tophub <- hub[, head(.SD,3), by=state]
topgene <- top_genes[, head(.SD,3), by=state]
L <- c(
  "=====================================================================",
  " STEP 4 / AIM 2 — per-state miR->gene interaction networks",
  "=====================================================================",
  "",
  "WHAT THIS ANSWERS: what each miRNA targets at each CD8 state. Built from the",
  "step-3 annotated reproducible clusters, restricted to PROTEIN-CODING targets.",
  "",
  "--- NETWORK SIZE PER STATE -----------------------------------------",
  apply(ntargets, 1, function(r) sprintf("    %-5s  %s edges | %s genes | %s miRNAs",
                                         r[["state"]], r[["edges"]], r[["genes"]], r[["miRNAs"]])),
  "",
  "--- HUB miRNAs (top 3 by targeting breadth) ------------------------",
  apply(tophub, 1, function(r) sprintf("    %-5s  %-18s %s targets",
                                       r[["state"]], sub("mmu-","",r[["miRNA"]]), r[["n_targets"]])),
  "",
  "--- MOST-TARGETED GENES (top 3 by miRNA convergence) ---------------",
  apply(topgene, 1, function(r) sprintf("    %-5s  %-12s %s miRNAs",
                                        r[["state"]], r[["gene"]], r[["n_miRNAs"]])),
  "",
  "--- FIGURES --------------------------------------------------------",
  "  N1 most-targeted genes per state (miRNA convergence)",
  "  N2 hub miRNAs per state (targeting breadth)",
  "  N3 curated network PNG per state: N3_network_<state>.png (all 4 states)",
  "     + interactive full network per state: network_<state>.html (open in browser)",
  "     + Cytoscape file per state: network_<state>.graphml",
  "  N4 master-TF targeting map across states (Tcf7/Lef1/Bach2/Id2/Tox/...)",
  "  N5 strongest target-switching miRNAs (low target-set overlap across states)",
  "",
  "--- LEADS CARRIED FROM STEP 3 --------------------------------------",
  "  * miR-466i -> Tcf7/Lef1/Bach2 (naive stemness TFs) - check N4/N1 (IP7 caveat).",
  "  * miR-142a-5p -> Id2 3'UTR (memory/exhausted) - check N4.",
  "  * Target switching (N5) is the core of aim 3: same miRNA redeployed onto a",
  "    different gene set per state = differentiation-stage-specific regulation.",
  "",
  "--- CAVEATS --------------------------------------------------------",
  "  * Protein-coding only; ncRNA/structural-RNA targets excluded here.",
  "  * Exhausted network is smallest (low capture, step-1 QC) - compare breadth",
  "    cautiously; normalize interpretation against capture depth.",
  "  * Edge weight = # reproducible clusters (shallow, median ~5 reads each).",
  "",
  "--- CSVs TO UPLOAD (eCLIP/04_networks/tables) ----------------------",
  "  top_targeted_genes_per_state.csv, hub_miRNAs_per_state.csv,",
  "  masterTF_targeting.csv, target_switching_summary.csv,",
  "  target_switching_pairwise.csv, interactions_by_state.csv",
  "",
  "NEXT (step 5, aim 3): differential targeting across states using the",
  "  Eclipsebio DESeq2 comparison TSVs + the switching results, focused on",
  "  differentiation master regulators.",
  "=====================================================================")
writeLines(L, file.path(OUT, "README_04_networks.txt"))
message("04_target_networks COMPLETE -> ", OUT)

