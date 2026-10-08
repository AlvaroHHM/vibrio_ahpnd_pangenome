#!/bin/bash
# process_fisher.sh – Filters significant genes from Fisher results and generates the supplementary table.

# Input files (edit names if needed)
FISHER_FILE="${1:-fisher_results_corrected.tsv}"
ANNOTATION="PIRATE.gene_families.ordered.tsv"
SIG_IDS="sig_ids_temp.txt"
ANNOTATED="sig_annotated_temp.tsv"
OUTPUT="supplementary_table_final.tsv"

echo "==> Using Fisher file: $FISHER_FILE"

# 1. Extract significant genes (p < 0.05)
awk -F'\t' 'NR>1 && $7 < 0.05' "$FISHER_FILE" > significant_temp.tsv
wc -l < significant_temp.tsv | xargs echo "Significant genes found:"

# 2. Extract gene IDs (first column)
tail -n +2 significant_temp.tsv | cut -f1 > "$SIG_IDS"

# 3. Annotate with product from PIRATE (column 2 = gene_family, column 4 = consensus_product)
awk -F'\t' 'NR==FNR {ids[$1]; next} $2 in ids {print $2 "\t" $4}' "$SIG_IDS" "$ANNOTATION" > "$ANNOTATED"

# 4. Merge annotation with significant results and add header
echo -e "gene\tpresent_in_AHPND+\tabsent_in_AHPND+\tpresent_in_AHPND-\tabsent_in_AHPND-\todds_ratio\tp_value\tproduct" > "$OUTPUT"
awk -F'\t' 'NR==FNR {annot[$1]=$2; next} {print $0 "\t" annot[$1]}' "$ANNOTATED" significant_temp.tsv >> "$OUTPUT"

# 5. Classify genes into gained/lost (based on counts)
awk -F'\t' 'NR==1 {print $0 "\tpattern"; next}
    {
        if      ($2>=5 && $4<=1) pat="Gained in AHPND+"
        else if ($2<=1 && $4>=3) pat="Lost in AHPND+"
        else                     pat="Other"
        print $0 "\t" pat
    }' "$OUTPUT" > tmp_classified.tsv
mv tmp_classified.tsv "$OUTPUT"

echo "==> Final supplementary table: $OUTPUT"
