###############################################################################
## 03_integration_mouse_CARTmiR29a.R
## Murine pre-infusion CAR-T (miR-29a) — single-cell, GEX only.
## Runs on the cell-cycle-scored, doublet-clean object from script 02.
##
## (EV/Scr/miR29a) is CONFOUNDED with capture (one library per condition per
## replicate), the choice of batch variable can either preserve or erase the
## miR-29a effect. So we run BOTH batch strategies under identical settings and
## compare:
##   - batch_var = "replicate"  -> removes run-date effect, KEEPS condition
##   - batch_var = "orig.ident" -> treats each capture as a batch; risks
##                                 collapsing condition (the confound)
##
## Integration methods compared: HARMONY (gentle, embedding-space correction)
## and FastMNN (mutual-nearest-neighbour correction). Everything except
## batch_var is held constant (same HVGs, dims, resolution, cc_regression), so any
## difference in the outputs is attributable to the batch choice. Feed both
## saved objects into the pseudobulk/composition referee (script 04).
##
###############################################################################

# ---- Libraries ----
library(Seurat)
library(SeuratWrappers)  # FastMNNIntegration
library(SeuratExtend)    # DimPlot2
library(tidyverse)
library(patchwork)
library(harmony)         # HarmonyIntegration
library(clustree)
library(qs2)

set.seed(123)

# ============================ Parameters (held constant) =====================
batch_vars        <- c("replicate", "orig.ident")  # run both, compare
cc_regression     <- "difference" # "none" | "full" | "difference"
#  none       = keep all cell-cycle signal
#  full       = regress S.Score + G2M.Score (removes ALL proliferation)
#  difference = regress S.Score - G2M.Score (keeps cycling-vs-resting,
#               removes only S-vs-G2M phase fragmentation)  [recommended here]
dims              <- 1:30
nfeatures         <- 2000
final_res         <- 0.4          # chosen from clustree (stable; 0.6 over-split)
primary_reduction <- "harmony"    # "harmony" or "integrated.mnn"

# ============================ Paths ==========================================
project_root <- "/home/akshay-iyer/Documents/Mouse-Pre-Infusion-CAR-T-miR-29a"
saved_dir    <- file.path(project_root, "saved_R_data")
integ_dir    <- file.path(project_root, "Integration")
dir.create(integ_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load doublet-clean object ONCE; each run starts from this clean base ----
base_obj <- qs_read(file.path(saved_dir,
                              "Mouse_CARTmiR29a_CellCycle_DoubletClean.qs2"))
DefaultAssay(base_obj) <- "RNA"
base_obj[["RNA"]] <- JoinLayers(base_obj[["RNA"]])   # guarantee a clean start

# cell-cycle regression target (computed once on the base object)
base_obj$CC.Difference <- base_obj$S.Score - base_obj$G2M.Score
vars_regress <- switch(cc_regression,
                       none       = NULL,
                       full       = c("S.Score", "G2M.Score"),
                       difference = "CC.Difference",
                       stop("cc_regression must be 'none', 'full', or 'difference'"))
message("Cell-cycle regression mode: ", cc_regression,
        " (vars.to.regress = ",
        if (is.null(vars_regress)) "none" else paste(vars_regress, collapse = ", "), ")")

# ============================ Integration routine ============================
run_integration <- function(obj, batch_var) {
  
  tag       <- gsub("[^A-Za-z0-9]", "", batch_var)    # "replicate" / "origident"
  out_dir   <- file.path(integ_dir, paste0("by_", tag))
  ctree_dir <- file.path(out_dir, "Clustree")
  dir.create(ctree_dir, recursive = TRUE, showWarnings = FALSE)
  message("\n==================  INTEGRATING over: ", batch_var, "  ==================")
  
  # --- split by batch + identical preprocessing ---
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[[batch_var]])
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  obj <- ScaleData(obj, vars.to.regress = vars_regress)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  
  # unintegrated baseline
  obj <- FindNeighbors(obj, reduction = "pca", dims = dims)
  obj <- FindClusters(obj, resolution = final_res, cluster.name = "unintegrated_clusters")
  obj <- RunUMAP(obj, reduction = "pca", dims = dims, reduction.name = "umap.unintegrated")
  
  # --- integrate: Harmony + FastMNN ---
  obj <- IntegrateLayers(obj, method = HarmonyIntegration,
                         orig.reduction = "pca",
                         new.reduction = "harmony", verbose = FALSE)
  obj <- IntegrateLayers(obj, method = FastMNNIntegration,
                         new.reduction = "integrated.mnn", verbose = FALSE)
  
  obj <- FindNeighbors(obj, reduction = "harmony", dims = dims)
  obj <- FindClusters(obj, resolution = final_res, cluster.name = "harmony_clusters")
  obj <- RunUMAP(obj, reduction = "harmony", dims = dims, reduction.name = "umap.harmony")
  
  obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = dims)
  obj <- FindClusters(obj, resolution = final_res, cluster.name = "mnn_clusters")
  obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = dims, reduction.name = "umap.mnn")
  
  # --- method comparison grid: left col = batch (should MIX),
  #     right col = condition (should STAY structured if real) ---
  mk <- function(red, grp, ttl)
    DimPlot2(obj, reduction = red, group.by = grp, label = FALSE, pt.size = 0.4) +
    ggtitle(ttl)
  
  comparison <- wrap_plots(
    mk("umap.unintegrated", batch_var,   paste0("Unintegrated - ", batch_var)),
    mk("umap.unintegrated", "condition", "Unintegrated - condition"),
    mk("umap.harmony",      batch_var,   paste0("Harmony - ", batch_var)),
    mk("umap.harmony",      "condition", "Harmony - condition"),
    mk("umap.mnn",          batch_var,   paste0("FastMNN - ", batch_var)),
    mk("umap.mnn",          "condition", "FastMNN - condition"),
    ncol = 2
  )
  ggsave(file.path(out_dir, "Integration_method_comparison.png"),
         comparison, width = 14, height = 18, dpi = 300)
  
  # Phase on primary reduction (cell cycle is regressed per cc_regression above)
  ggsave(file.path(out_dir, paste0("UMAP_", primary_reduction, "_Phase.png")),
         DimPlot2(obj, reduction = paste0("umap.",
                                          ifelse(primary_reduction == "harmony", "harmony", "mnn")),
                  group.by = "Phase", pt.size = 0.4),
         width = 11, height = 8, dpi = 300)
  
  # --- clustree resolution sweep on primary reduction ---
  obj <- FindNeighbors(obj, reduction = primary_reduction, dims = dims)
  for (res in seq(0.2, 1.2, by = 0.2))
    obj <- FindClusters(obj, resolution = res)        # writes RNA_snn_res.<res>
  
  ggsave(file.path(ctree_dir, "clustree.png"),
         clustree(obj, prefix = "RNA_snn_res."),
         width = 15, height = 9, dpi = 300, bg = "white")
  
  # lock working clusters at chosen resolution
  umap_primary <- ifelse(primary_reduction == "harmony", "umap.harmony", "umap.mnn")
  obj$clusters <- obj[[paste0("RNA_snn_res.", final_res)]][, 1]
  # order cluster levels NUMERICALLY (0,1,2,...,10 — not 0,1,10,...,2)
  lv <- levels(factor(obj$clusters))
  obj$clusters <- factor(obj$clusters, levels = lv[order(as.numeric(lv))])
  Idents(obj)  <- "clusters"
  
  ggsave(file.path(out_dir, paste0("UMAP_clusters_res", final_res, ".png")),
         DimPlot2(obj, reduction = umap_primary,
                  group.by = "clusters", label = TRUE, box = TRUE, repel = TRUE,
                  label.color = "black", pt.size = 0.6),
         width = 11, height = 8, dpi = 300)
  
  # --- composition: each cluster's makeup by condition (the flattening check) ---
  comp <- obj@meta.data %>%
    dplyr::count(clusters, condition) %>%
    dplyr::group_by(clusters) %>%
    dplyr::mutate(frac = n / sum(n)) %>%
    dplyr::ungroup()
  write.csv(comp, file.path(out_dir, "Cluster_composition_by_condition.csv"),
            row.names = FALSE)
  
  ggsave(file.path(out_dir, "Cluster_composition_by_condition.png"),
         ggplot(comp, aes(x = clusters, y = frac, fill = condition)) +
           geom_col() + theme_classic() +
           labs(y = "fraction of cluster", x = "cluster",
                title = paste0("Cluster composition - integrated by ", batch_var)),
         width = 12, height = 6, dpi = 300)
  
  # --- rejoin RNA layers (single matrix for downstream DE) + save ---
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  out_file <- file.path(saved_dir,
                        paste0("Mouse_CARTmiR29a_integrated_by_", tag, ".qs2"))
  qs_save(obj, file = out_file)
  message("Saved -> ", basename(out_file),
          " | cells: ", ncol(obj),
          " | clusters @res ", final_res, ": ", length(unique(obj$clusters)))
  
  invisible(NULL)
}

# ============================ Run both strategies ============================
for (bv in batch_vars) {
  run_integration(base_obj, bv)   # fresh copy each call -> identical start
  gc()
}

message("\nDone. Compare:\n",
        "  Integration/by_replicate/   vs   Integration/by_origident/\n",
        "Watch the composition bars: if 'origident' flattens clusters toward\n",
        "even EV/Scr/miR29a fractions while 'replicate' keeps real skew, that\n",
        "skew IS the condition signal orig.ident integration removes.\n",
        "Feed both saved objects to the pseudobulk/composition referee (script 04).")
###############################################################################