=====================================================================
 STEP 1 - QUALITY CONTROL : miR-eCLIP CD8 differentiation
=====================================================================

GROUND TRUTH = THE DATA FILES. The Eclipsebio HTML report supplies only
read-level metrics absent from the BEDs (depth, trim, mapping, PCR-dup,
% chimeras). Its cluster counts are a summary and DO NOT match the files.

GENOME: UCSC GRCm38/mm10, GENCODE vM25 (differs from CAR-T's GRCm39 - step 3
annotation must use a vM25/mm10 GTF).

DESIGN: 4 CD8 states (Naive Teff Tmem Tex), 2 IP + 2 input each. AGO2 IP +
ligation -> chimeric reads = DIRECT miR->target pairs (a small % of reads).

OUTPUT (eCLIP/01_QC/): QC1..QC8 PNGs in the folder; all CSVs in 01_QC/tables/.

--- DATA vs REPORT: chimeric clusters (per IP) ----------------------
  The .mir_targets_clusters.bed files carry ~30% MORE rows than the report's
  'Chimeric clusters' summary. Analysis runs on the BED rows; treat report
  values as approximate. Per library (observed in file | report | diff):
    IP7_Tn_1_2_RR2     8004 |   6311 |   +27%
    IP8_Tn_3_RR2       1869 |   1514 |   +23%
    IP1_Teff_1_RR2     1039 |    780 |   +33%
    IP3_Teff_2_RR2      886 |    586 |   +51%
    IP5_Tmem_1_RR2     1932 |   1399 |   +38%
    IP6_Tmem_2_RR2     1985 |   1351 |   +47%
    IP2_Tex_1_RR2       463 |    336 |   +38%
    IP4_Tex_2_RR2       939 |    691 |   +36%

--- REPRODUCIBLE YIELD PER STATE (from files) ----------------------
    Tn       1373 reproducible clusters |   68 distinct miRNAs
    Teff      549 reproducible clusters |   81 distinct miRNAs
    Tmem     1152 reproducible clusters |  102 distinct miRNAs
    Tex       322 reproducible clusters |   65 distinct miRNAs

--- REPLICATE CONCORDANCE (Spearman, per-miRNA coverage) -----------
    Tn     rho = 0.71 (n=151 miRNAs)
    Teff   rho = 0.87 (n=118 miRNAs)
    Tmem   rho = 0.86 (n=150 miRNAs)
    Tex    rho = 0.79 (n=115 miRNAs)

--- KEY FLAGS ------------------------------------------------------
  * NAIVE IS UNBALANCED: IP7_Tn_1_2 >> IP8_Tn_3 in chimeric reads/clusters.
    Raw per-library counts are not comparable. Defenses: reproducible set
    (must be in BOTH reps) + within-state normalization (step 2).
  * Tex has the lowest chimeric yield - absolute target counts into exhaustion
    read low partly for capture reasons, not only biology.
  * Chimeric targets span CDS/intron/3'UTR: do NOT restrict to 3'UTR downstream.

--- BAM CHECK (IP1_Teff_1 .mir_targets.bam inspected) --------------
  The chimeric BAM is STAR-aligned target reads (mm10) with a UMI in the read
  name and NO miRNA tag - miRNA identity is not in the BAM. The cluster BED
  'coverage' column IS the definitive per-miRNA read support; no BAM needed
  for aim 1. (BAM read count = report 'Final chimeric reads', a useful check.)

NEXT (step 2, aim 1): per-miRNA engaged abundance per state, normalized for the
  naive imbalance, from the reproducible set + per-sample coverage.
=====================================================================
