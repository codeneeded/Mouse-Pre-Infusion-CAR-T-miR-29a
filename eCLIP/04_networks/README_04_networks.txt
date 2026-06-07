=====================================================================
 STEP 4 / AIM 2 — per-state miR->gene interaction networks
=====================================================================

WHAT THIS ANSWERS: what each miRNA targets at each CD8 state. Built from the
step-3 annotated reproducible clusters, restricted to PROTEIN-CODING targets.

--- NETWORK SIZE PER STATE -----------------------------------------
    Tn     543 edges | 427 genes | 49 miRNAs
    Teff   273 edges | 208 genes | 59 miRNAs
    Tmem   771 edges | 536 genes | 82 miRNAs
    Tex    174 edges | 137 genes | 48 miRNAs

--- HUB miRNAs (top 3 by targeting breadth) ------------------------
    Tn     miR-466i-5p        317 targets
    Tn     miR-142a-3p         50 targets
    Tn     miR-466f-3p         20 targets
    Teff   miR-466i-5p         75 targets
    Teff   miR-142a-3p         26 targets
    Teff   miR-466f-3p         24 targets
    Tmem   miR-466i-5p         96 targets
    Tmem   miR-142a-3p         79 targets
    Tmem   miR-92a-3p          63 targets
    Tex    miR-466i-5p         36 targets
    Tex    miR-142a-3p         16 targets
    Tex    miR-466f-3p         13 targets

--- MOST-TARGETED GENES (top 3 by miRNA convergence) ---------------
    Tn     Dock2        11 miRNAs
    Tn     Lars2         6 miRNAs
    Tn     Hexb          5 miRNAs
    Teff   Actb          5 miRNAs
    Teff   Rpl13a        4 miRNAs
    Teff   Gimap4        4 miRNAs
    Tmem   mt-Nd1       13 miRNAs
    Tmem   Actb         11 miRNAs
    Tmem   Pim1          9 miRNAs
    Tex    mt-Nd1        7 miRNAs
    Tex    Hexb          5 miRNAs
    Tex    Akna          3 miRNAs

--- FIGURES --------------------------------------------------------
  N1 most-targeted genes per state (miRNA convergence)
  N2 hub miRNAs per state (targeting breadth)
  N3 curated network PNG per state: N3_network_<state>.png (all 4 states)
     + interactive full network per state: network_<state>.html (open in browser)
     + Cytoscape file per state: network_<state>.graphml
  N4 master-TF targeting map across states (Tcf7/Lef1/Bach2/Id2/Tox/...)
  N5 strongest target-switching miRNAs (low target-set overlap across states)

--- LEADS CARRIED FROM STEP 3 --------------------------------------
  * miR-466i -> Tcf7/Lef1/Bach2 (naive stemness TFs) - check N4/N1 (IP7 caveat).
  * miR-142a-5p -> Id2 3'UTR (memory/exhausted) - check N4.
  * Target switching (N5) is the core of aim 3: same miRNA redeployed onto a
    different gene set per state = differentiation-stage-specific regulation.

--- CAVEATS --------------------------------------------------------
  * Protein-coding only; ncRNA/structural-RNA targets excluded here.
  * Exhausted network is smallest (low capture, step-1 QC) - compare breadth
    cautiously; normalize interpretation against capture depth.
  * Edge weight = # reproducible clusters (shallow, median ~5 reads each).

--- CSVs TO UPLOAD (eCLIP/04_networks/tables) ----------------------
  top_targeted_genes_per_state.csv, hub_miRNAs_per_state.csv,
  masterTF_targeting.csv, target_switching_summary.csv,
  target_switching_pairwise.csv, interactions_by_state.csv

NEXT (step 5, aim 3): differential targeting across states using the
  Eclipsebio DESeq2 comparison TSVs + the switching results, focused on
  differentiation master regulators.
=====================================================================
