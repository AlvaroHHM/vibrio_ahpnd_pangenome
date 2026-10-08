#!/bin/bash
# ============================================================
# generar_README.sh
# Crea README.md dentro de Supplementary_Files/
# ============================================================

DEST=~/bacterial-genomics-tutorial/Supplementary_Files
mkdir -p "$DEST"

cat > "$DEST/README.md" << 'EOF'
# Supplementary Files

This folder contains all supplementary tables and figures referenced in the main manuscript.

Below is a detailed description of each item, the analysis from which it was generated, and a suggested caption.

## Supplementary Tables

### Table S1. Strains used in this study
- **File:** `Table_S1_strains.tsv`
- **Content:** strain name, clade/group, AHPND status, host, origin, year, total genes, accessions.
- **Suggested legend:** *"Vibrio parahaemolyticus strains used in this study..."*

### Table S2. Candidate accessory-genome families
- **File:** `Table_S2_candidate_genes.tsv`
- **Content:** gene family, counts, odds ratio, 95% CI, p-value, FDR.
- **Suggested legend:** *"Candidate accessory-genome families with nominal association..."*

### Table S3. Sensitivity analyses excluding L2181 and R13
- **Files:** `Table_S3_sensitivity_excluding_L2181.tsv`, `Table_S3_sensitivity_excluding_R13_annotated.tsv`
- **Content:** Fisher results after excluding L2181 or R13.
- **Suggested legend:** *"Sensitivity analyses of the gene-phenotype association..."*

### Table S4. Validation of selected gene absences
- **File:** `Table_S4_validation_absences.tsv`
- **Content:** gene, strain, BLASTp/TBLASTN hits, synteny, mapping, PCR, conclusion.
- **Suggested legend:** *"Validation of selected gene absences in 15_CESAIBC..."*

### Table S5. Replicate-level mortality data
- **File:** `Table_S5_mortality_by_replicate.tsv`
- **Content:** strain, replicate, initial n, deaths, mortality.
- **Suggested legend:** *"Replicate-level mortality data for the pathogenicity bioassay..."*

### Table S6. COG category enrichment analysis
- **Files:** `Table_S6_COG_enrichment.tsv`, `Table_S6_COG_enrichment_Roary.tsv`
- **Content:** COG category, counts, odds ratio, p-value, FDR.
- **Suggested legend:** *"Enrichment analysis of COG categories..."*

### Table S7. Annotated genes significant excluding R13
- **File:** `Table_S7_annotated_15genes_excluding_R13.tsv`
- **Content:** 15 genes with FDR < 0.05 after excluding R13, with annotations.
- **Suggested legend:** *"Annotated genes with FDR < 0.05 when R13 was excluded..."*

## Supplementary Figures

### Figure S1. Randomized accumulation curves
- **File:** `Figure_S1_randomized_accumulation_curves.png`
- **Suggested legend:** *"Randomized accumulation curves for the pangenome and core genome..."*

### Figure S2. Ancestral state reconstruction of four marker genes
- **Files:** `Figure_S2_ancestral_reconstruction_4genes.png`, `.svg`
- **Suggested legend:** *"Ancestral state reconstruction of the four genes perfectly associated with AHPND+..."*

### Figure S3. Heatmap of top 50 filtered genes
- **Files:** `Figure_S3_heatmap_top50_filtered.png`, `.svg`
- **Suggested legend:** *"Heatmap of the top 50 genes with the highest discriminatory score..."*

### Figure S4. Heatmap of COG V (defense) genes
- **Files:** `Figure_S4_heatmap_COG_V.png`, `.svg`
- **Suggested legend:** *"Heatmap of defense-related (COG category V) genes..."*

### Figure S5. Clustermap of top 50 genes
- **Files:** `Figure_S5_clustermap_top50.png`, `.svg`
- **Suggested legend:** *"Clustermap showing hierarchical clustering..."*

### Figure S6. Co-localization of pirA, pirB and associated genes
- **Files:** `Figure_S6_colocalization.png`, `.svg`
- **Suggested legend:** *"Co-localization of pirA, pirB, g07720, g06662, and g07221..."*

## Notes
- All analyses were performed in R v4.3.3 unless otherwise stated.
- External software included PIRATE, Prokka, Gubbins, IQ-TREE, MAFFT, BLAST, and eggNOG-mapper.
EOF

echo "README.md generado en $DEST"
