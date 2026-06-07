=====================================================================
 STEP 5 / AIM 3 — differential targeting across CD8 states
=====================================================================

LENS: Eclipsebio DESeq2 differential TSVs (IP-vs-input normalized), which
control for transcript abundance — the clean complement to step-4 raw
convergence. Abundant-RNA background (mt-/Rpl/Rps/Mrp/Actb/Malat1/Neat1)
is removed from gene-level views. log2FC sign = first state vs second.

Significant differential interactions (padj<0.05,|log2FC|>=1, non-bg): 154

--- GAINED / LOST PER CONTRAST -------------------------------------
    Teff_vs_Tn   n=27   gained=19   lost= 8
    Tmem_vs_Tn   n=85   gained=54   lost=31
    Tex_vs_Tn    n=29   gained=18   lost=11
    Tex_vs_Tmem  n= 5   gained= 4   lost= 1
    Tex_vs_Teff  n= 2   gained= 1   lost= 1
    Tmem_vs_Teff n= 6   gained= 3   lost= 3

--- MASTER-REGULATOR DIFFERENTIAL TARGETING ------------------------
    Bach2/miR-142a-3p      Teff_vs_Tn   log2FC -3.68
    Bach2/miR-466i-5p      Teff_vs_Tn   log2FC -3.02
    Bach2/miR-466i-5p      Teff_vs_Tn   log2FC -2.89
    Bach2/miR-466i-5p      Teff_vs_Tn   log2FC -2.5
    Bach2/miR-466i-5p      Tmem_vs_Tn   log2FC -4.24
    Id2/miR-142a-5p        Tmem_vs_Tn   log2FC 3.45
    Bach2/miR-466i-5p      Tmem_vs_Tn   log2FC -3.73
    Bach2/miR-466i-5p      Tmem_vs_Tn   log2FC -3.2
    Bach2/miR-466i-5p      Tmem_vs_Tn   log2FC -3.08
    Bach2/miR-142a-3p      Tmem_vs_Tn   log2FC -3.12
    Bach2/miR-466i-5p      Tex_vs_Tn    log2FC -3.61
    Bach2/miR-142a-3p      Tex_vs_Tn    log2FC -3.53
    Bach2/miR-466i-5p      Tex_vs_Tn    log2FC -3.37
    Bach2/miR-466i-5p      Tex_vs_Tn    log2FC -2.96
    Id2/miR-142a-5p        Tex_vs_Tn    log2FC 2.72

--- FIGURES --------------------------------------------------------
  D1 volcanoes per contrast (master TFs labelled)
  D2 gained vs lost targeting per contrast (directionality)
  D3 master-regulator differential targeting heatmap
  D4 target switching: top switchers' genes per state (miR-210/-155/-23a/-150)
  D5 targeting gained INTO exhaustion (Tex vs others)

--- CAVEATS --------------------------------------------------------
  * n=2/state (Eclipsebio DESeq2): treat p-values as exploratory; weight by
    effect size + biological coherence.
  * The differential TSVs are the abundance-controlled view; cross-check
    candidate targets against step-4 reproducible interactions.
  * Tex has lowest capture (step-1 QC): 'lost into exhaustion' can be partly
    technical — interpret losses more cautiously than gains.

--- CSVs TO UPLOAD (eCLIP/05_differential/tables) ------------------
  significant_differential_targeting.csv, gained_lost_targeting_per_contrast.csv,
  masterTF_differential_targeting.csv, top_switcher_targets_by_state.csv,
  exhaustion_gained_targeting_genes.csv

NEXT (step 6, aim 4): therapeutic candidates — miRNAs targeting exhaustion
  drivers / effector programs, ranked; miR-29 family revisited at relaxed
  support against the CAR-T DE targets.
=====================================================================
