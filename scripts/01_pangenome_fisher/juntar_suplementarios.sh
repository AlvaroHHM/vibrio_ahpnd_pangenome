#!/bin/bash
# ============================================================
# juntar_suplementarios.sh
# Copia tablas y figuras suplementarias a Supplementary_Files/
# ============================================================

BASE=~/bacterial-genomics-tutorial
DEST="$BASE/Supplementary_Files"
mkdir -p "$DEST"

copy_if_exists() {
    local src="$1"
    local dest_name="$2"
    if [ -f "$src" ]; then
        cp "$src" "$DEST/$dest_name"
        echo "✔ Copiado: $dest_name"
    else
        echo "✘ No encontrado: $src"
    fi
}

echo "=== Tablas suplementarias ==="

copy_if_exists "$BASE/pangenome/Supplementary_Table_exclusivos_ausentes.tsv" "Table_S1_strains.tsv"
copy_if_exists "$BASE/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3/R_analysis/Table_S2_candidatos.tsv" "Table_S2_candidate_genes.tsv"
copy_if_exists "$BASE/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3/sensitivity_excluding_L2181.tsv" "Table_S3_sensitivity_excluding_L2181.tsv"
copy_if_exists "$BASE/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3/sensitivity_excluding_R13_FDR005_annotated.tsv" "Table_S3_sensitivity_excluding_R13_annotated.tsv"
copy_if_exists "$BASE/reference_genomes_pipeline/pirate_input_pirAB/Table_S4_validation_absences.tsv" "Table_S4_validation_absences.tsv"
copy_if_exists "$BASE/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3/Table_S5_mortality_template.tsv" "Table_S5_mortality_by_replicate.tsv"
copy_if_exists "$BASE/pangenome/enrichment_cog.tsv" "Table_S6_COG_enrichment.tsv"
copy_if_exists "$BASE/pangenome/enriquecimiento_COG_roary.tsv" "Table_S6_COG_enrichment_Roary.tsv"
copy_if_exists "$BASE/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3/sensitivity_excluding_R13_FDR005_annotated.tsv" "Table_S7_annotated_15genes_excluding_R13.tsv"

echo ""
echo "=== Figuras suplementarias ==="

copy_if_exists "$BASE/pangenome/curvas_acumulacion_aleatorizadas.png" "Figure_S1_randomized_accumulation_curves.png"
copy_if_exists "$BASE/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3/ancestral_svg/panel_ancestral_4genes.png" "Figure_S2_ancestral_reconstruction_4genes.png"
copy_if_exists "$BASE/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3/ancestral_svg/panel_ancestral_4genes.svg" "Figure_S2_ancestral_reconstruction_4genes.svg"
copy_if_exists "$BASE/pangenome/R_analysis_huella/heatmap_top50_pirate_R.png" "Figure_S3_heatmap_top50_filtered.png"
copy_if_exists "$BASE/pangenome/R_analysis_huella/heatmap_top50_pirate_R.svg" "Figure_S3_heatmap_top50_filtered.svg"
copy_if_exists "$BASE/pangenome/heatmap_COG_V.png" "Figure_S4_heatmap_COG_V.png"
copy_if_exists "$BASE/pangenome/heatmap_COG_V.svg" "Figure_S4_heatmap_COG_V.svg"
copy_if_exists "$BASE/pangenome/R_analysis_huella/clustermap_top50_pirate_R.png" "Figure_S5_clustermap_top50.png"
copy_if_exists "$BASE/pangenome/R_analysis_huella/clustermap_top50_pirate_R.svg" "Figure_S5_clustermap_top50.svg"
copy_if_exists "$BASE/reference_genomes_pipeline/pirate_input_pirAB/Figure_S6_colocalization.png" "Figure_S6_colocalization.png"
copy_if_exists "$BASE/reference_genomes_pipeline/pirate_input_pirAB/Figure_S6_colocalization.svg" "Figure_S6_colocalization.svg"

echo ""
echo "Archivos copiados a: $DEST"
