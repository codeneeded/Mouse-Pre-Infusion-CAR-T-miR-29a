=====================================================================
 STEP 6 / AIM 4 — therapeutic / state-pushing candidate miRNAs
=====================================================================

FRAMING (from aim 3): no exhaustion-SPECIFIC target program exists, so
candidates are nominated by the EFFECTOR / MEMORY programs they target and
their loading in the relevant state, not by exhaustion-only targets.
Abundant-RNA background (incl. Gas5/Snhg/Snord, extended this step) removed.

--- TOP EFFECTOR-PROGRAM-TARGETING miRNAs --------------------------
    miR-142a-5p        genes: Id2
    NA                 genes: NA
    NA                 genes: NA
--- TOP MEMORY-PROGRAM-TARGETING miRNAs ----------------------------
    miR-466i-5p        genes: Bach2,Bcl2,Il7r,Lef1,Tcf7
    miR-15b-5p         genes: Bcl2,Ccr7,Il7r
    miR-142a-3p        genes: Bach2,Bcl2

--- miR-29 / CAR-T CROSS-VALIDATION (relaxed, per-sample) ----------
  CAR-T DE targets (Dnmt3a,Tet2,Eomes,Tbx21,Ifng) with ANY direct miR-29 chimera: NONE
  Even at relaxed (single-replicate) support, the canonical miR-29 targets
  are largely/entirely absent as direct chimeras in these CD8 states.
  -> The CAR-T miR-29a phenotype is likely NOT explained by direct binding of
     these targets in CD8 T cells; see T4 for what miR-29 DOES bind here.

--- FIGURES --------------------------------------------------------
  T1 candidate miRNAs by program (effector/memory/exhaustion)
  T2 program-gene targeting map (gene x miRNA)
  T3 exhaustion-leverage candidates (effector/memory targeters loaded in Tex)
  T4 miR-29 family target landscape across states (relaxed)
  T5 CAR-T miR-29a DE target check (found / not detected)

--- CAVEATS --------------------------------------------------------
  * 'Therapeutic candidate' = hypothesis-generating; binding != functional
    repression. Validate with the matched RNA-seq (still the key missing
    data) before any therapeutic claim.
  * Relaxed miR-29 search uses single-replicate clusters: more sensitive but
    less stringent; absence here is informative given that relaxation.
  * n=2; Tex low capture; miR-466i naive inflation - as in prior steps.

--- CSVs TO UPLOAD (eCLIP/06_therapeutic/tables) -------------------
  therapeutic_candidates_ranked.csv, miR_targeting_program_genes.csv,
  exhaustion_leverage_candidates.csv, miR29_targets_relaxed_persample.csv,
  miR29_CART_targets_summary.csv

PROJECT STATUS: aims 1-4 complete at the BINDING level. The remaining gap for
functional/therapeutic conclusions is matched differentiation-state RNA-seq
(target derepression). That is the recommended next experiment.
=====================================================================
