#!/bin/bash
# filter_extreme_annotate.sh
# Filtra genes con odds ratio extremo (inf o 0/0.0/0.000) y anota con producto génico.

FISHER_FILE="significant_fisher_results.tsv"
PIRATE_FILE="PIRATE.gene_families.ordered.tsv"
OUTPUT="significant_extreme_annotated.tsv"

# 1. Crear mapa familyID -> product
awk -F'\t' 'NR>1 {print $2 "\t" $4}' "$PIRATE_FILE" | sort -u > family_product_map.tsv

# 2. Filtrar y anotar
awk -F'\t' '
BEGIN {
    while ((getline line < "family_product_map.tsv") > 0) {
        split(line, arr, "\t");
        product[arr[1]] = arr[2];
    }
    print "gene\tpresent_in_AHPND+\tabsent_in_AHPND+\tpresent_in_AHPND-\tabsent_in_AHPND-\todds_ratio\tp_value\tproduct"
}
NR>1 {
    # Captura inf y cualquier forma de 0: 0, 0.0, 0.000, etc.
    if ($6 == "inf" || $6 ~ /^0(\.0+)?$/) {
        prod = ($1 in product) ? product[$1] : "NA";
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, prod;
    }
}' "$FISHER_FILE" > "$OUTPUT"

rm -f family_product_map.tsv
echo "Listo. Archivo generado: $OUTPUT"
