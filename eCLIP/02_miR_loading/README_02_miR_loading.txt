=====================================================================
 STEP 2 / AIM 1 — per-miRNA engaged abundance per CD8 state
=====================================================================

WHAT THIS ANSWERS: which miRNAs are abundantly engaged (loaded + actively
forming chimeras with targets) at each state, and how loading shifts across
naive -> effector -> memory and into exhaustion.

QUANTITY & NORMALIZATION:
  - Unit = chimeric-cluster COVERAGE per miRNA (per step-1 QC, the BAM has no
    miRNA identity, so coverage is the definitive per-miRNA read support).
  - Universe = 120 miRNAs reproducible in >=1 state; 120 pass the count
    prefilter and are tested.
  - DESeq2 size factors normalize per-library coverage -> removes the naive
    IP7 depth spike so states are comparable.
  - n=2/state: DIFFERENTIAL LOADING IS EXPLORATORY (ranking aid, not proof).

--- TOP ENGAGED miRNAs PER STATE (within-state %) ------------------
    Tn: miR-466i-5p (35.7%), miR-24-3p (12.4%), miR-142a-3p (11.3%), miR-26a-5p (3.2%), miR-92a-3p (3.2%)
    Teff: miR-24-3p (25.1%), miR-466i-5p (12.0%), miR-21a-5p (11.1%), miR-142a-3p (7.0%), miR-25-3p (4.9%)
    Tmem: miR-24-3p (31.0%), miR-155-5p (8.7%), miR-466i-5p (6.7%), miR-92a-3p (6.4%), miR-142a-3p (6.0%)
    Tex: miR-155-5p (21.6%), miR-24-3p (18.8%), miR-466i-5p (9.3%), miR-142a-3p (5.6%), miR-25-3p (5.3%)

--- SIGNIFICANT DIFFERENTIAL LOADING (padj<0.05 & |log2FC|>=1) ------
    Teff_vs_Tn   14 miRNAs
    Tmem_vs_Tn   17 miRNAs
    Tex_vs_Tn    26 miRNAs
    Tex_vs_Teff   2 miRNAs

--- HOW TO READ THE FIGURES ----------------------------------------
  M1 top engaged miRNAs per state (what dominates each state's signal)
  M2 z-scored heatmap of the most state-variable miRNAs (loading shifts)
  M3 tau state-specificity (which miRNAs are restricted to one state)
  M4 differential loading volcanoes across 5 contrasts (exploratory)
  M5 curated immunology miRNAs across the axis (incl. miR-29 family for the
     CAR-T link; miR-155/-31/-210 exhaustion; miR-150/let-7 naive/quiescence)

--- BIOLOGY TO CHECK -----------------------------------------------
  * miR-150 / let-7 high in naive (quiescence) and falling on activation?
  * miR-155 / miR-21 / miR-31 / miR-210 rising into exhaustion?
  * miR-17~92 (miR-17/19b/92a) up in effector (proliferation)?
  * miR-29a/b/c levels — feeds directly into the CAR-T miR-29a overexpression
    work (its targets are validated in step 6).
  NOTE: naive abundance is dominated by few miRNAs (e.g. miR-466i); interpret
  with the IP7/IP8 imbalance from step-1 QC in mind. Tex has the thinnest
  engaged layer, so low Tex abundance is partly capture, not only biology.

--- CSVs TO UPLOAD (eCLIP/02_miR_loading/tables) -------------------
  miRNA_abundance_by_state_long.csv, top_miRNAs_per_state.csv,
  miRNA_state_specificity.csv, differential_loading_all_contrasts.csv,
  curated_immune_miRNAs.csv

NEXT (step 3, aim 2 prep): annotate reproducible clusters to genes (GENCODE
  vM25 / mm10) so we know WHAT each miRNA targets per state.
=====================================================================
