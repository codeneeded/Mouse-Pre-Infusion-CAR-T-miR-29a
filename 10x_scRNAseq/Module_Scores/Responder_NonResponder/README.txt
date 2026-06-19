07b -- Responder vs Non-Responder module scores (CD4/CD8 subsets)

Inputs:
  saved_R_data/Mouse_CARTmiR29a_WithModuleScores.qs2  (script 07 output)
  Resources/CAR-T_Responder_v_Non-Responder_Meta.xlsx

Human->mouse conversion: babelgene orthologs + validated title-case fallback.
Scores computed on the FULL object, then plotted on CD4 / CD8 subsets
(subset only -- existing umap.harmony coords reused, no re-clustering).

Plots/UMAP_CD4, Plots/UMAP_CD8      per-set score FeaturePlots + net_* scores
Plots/VlnByCluster_CD4 / _CD8       per-cluster violins within each lineage
Plots/UMAP_headline_RvNR_CD4/CD8    Fraietta + Combined composite grids
Plots/Hypoxia_ROS_Stress            hypoxia/ROS/stress FeaturePlots + violins

Tables/RvNR_human_to_mouse_audit.csv         per-gene mapping + method + in_object
Tables/RvNR_conversion_summary.csv           per-set mapping counts
Tables/RvNR_module_genes_used.csv            mouse genes actually scored
Tables/RvNR_module_condition_wilcoxon_CD4_CD8.csv  miR29a vs EV/Scr per subset
Tables/hypoxia_ROS_stress_subset_cluster_summary.csv
Tables/hypoxia_ROS_stress_CD8_vs_CD4.csv
