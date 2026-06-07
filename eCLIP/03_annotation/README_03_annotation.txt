=====================================================================
 STEP 3 / AIM 2 prep — annotate reproducible clusters to genes
=====================================================================

WHAT THIS DOES: assigns every reproducible chimeric cluster a target GENE,
REGION and BIOTYPE by overlap with GENCODE vM25 (mm10). Strand-aware; the
~30% 'Other' bucket from step-1 QC is now resolved into sense ncRNA /
antisense / intergenic.

REGION PRIORITY (sense): 3'UTR > 5'UTR > CDS > ncRNA_exon > intron;
  no sense gene -> antisense (opposite-strand gene) -> intergenic.

--- ANNOTATION OUTCOME PER STATE (% sense / antisense / intergenic) -
    Tn     sense 100% | antisense 0% | intergenic 0%
    Teff   sense 100% | antisense 0% | intergenic 0%
    Tmem   sense 100% | antisense 0% | intergenic 0%
    Tex    sense 100% | antisense 0% | intergenic 0%

--- DISTINCT PROTEIN-CODING TARGET GENES PER STATE -----------------
    Tn     427 genes
    Teff   208 genes
    Tmem   536 genes
    Tex    137 genes

--- KEY GENES DIRECTLY TARGETED (preview of aims 2/4) --------------
  Found among reproducible targets: Bach2, Id2, Lef1, Tcf7
  miR-29 family -> CAR-T DE targets bound here: none directly (check step 4/6 with relaxed support)

--- FIGURES --------------------------------------------------------
  A1 region composition per state (Other resolved)
  A2 target gene biotype per state
  A3 distinct protein-coding target genes per state
  A4 key immune/CAR-T gene targeting preview (dot = #miRNAs binding)

--- CAVEATS --------------------------------------------------------
  * Clusters overlapping multiple genes keep the best (protein_coding, then
    widest overlap); multi-gene loci are simplified.
  * intergenic/antisense clusters are real eCLIP signal (AGO can bind ncRNA /
    nascent transcripts); they are retained, not discarded.
  * Annotation uses GENCODE vM25 (mm10) to MATCH the eCLIP processing genome,
    NOT the CAR-T GRCm39 reference.

--- CSVs TO UPLOAD (eCLIP/03_annotation/tables) --------------------
  region_composition_per_state.csv, biotype_composition_per_state.csv,
  unique_proteincoding_targets_per_state.csv, key_gene_targeting_preview.csv,
  annotation_outcome_per_state.csv

NEXT (step 4, aim 2): per-state miR->gene interaction networks, hub miRNAs,
  most-targeted genes, and target switching across states.
=====================================================================
