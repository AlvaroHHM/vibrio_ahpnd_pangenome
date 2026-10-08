import csv
from scipy.stats import fisher_exact

# ----- Leer fenotipos -----
pheno = {}
with open('phenotype.csv') as pf:
    reader = csv.reader(pf)
    header = next(reader)
    for row in reader:
        if not row: continue
        strain = row[0].strip()
        ahpnd = row[1].strip()
        pheno[strain] = 1 if ahpnd == '1' else 0

# ----- Leer matriz de presencia/ausencia -----
gene_data = {}      # gene -> lista de 0/1 en el orden de aparición de las cepas
strains_order = []  # mismo orden que las filas del CSV
with open('gene_presence_absence.csv') as gf:
    reader = csv.reader(gf)
    header = next(reader)
    genes = header[1:]  # todos los nombres de genes
    # Inicializar listas vacías
    for gene in genes:
        gene_data[gene] = []
    for row in reader:
        strain = row[0].strip()
        strains_order.append(strain)
        for i, gene in enumerate(genes, 1):
            gene_data[gene].append(int(row[i]))

# ----- Calcular test de Fisher para cada gen -----
all_results = []   # lista de tuplas (p_value, gene, a, b, c, d, odds)
for gene in genes:
    a = b = c = d = 0
    for strain, pres in zip(strains_order, gene_data[gene]):
        if strain not in pheno:
            continue
        if pheno[strain] == 1:   # AHPND+
            if pres == 1: a += 1
            else: b += 1
        else:                    # AHPND-
            if pres == 1: c += 1
            else: d += 1
    # Solo calcular si hay variabilidad (al menos un conteo en cada fila/columna)
    if (a + b) > 0 and (c + d) > 0 and (a + c) > 0 and (b + d) > 0:
        odds, p = fisher_exact([[a, b], [c, d]])
    else:
        odds, p = float('inf'), 1.0  # o NA, pero mantenemos 1 para evitar errores
    all_results.append((p, gene, a, b, c, d, odds))

# ----- Escribir tabla completa -----
with open('fisher_results.tsv', 'w') as out:
    out.write('gene\tpresent_in_AHPND+\tabsent_in_AHPND+\tpresent_in_AHPND-\tabsent_in_AHPND-\todds_ratio\tp_value\n')
    for (p, gene, a, b, c, d, odds) in all_results:
        out.write(f'{gene}\t{a}\t{b}\t{c}\t{d}\t{odds}\t{p}\n')

# ----- Filtrar significativos (p < 0.05) y ordenar por p_value -----
significant = [(p, gene, a, b, c, d, odds) for (p, gene, a, b, c, d, odds) in all_results if p < 0.05]
significant.sort(key=lambda x: x[0])   # orden creciente de p

with open('significant_fisher_results.tsv', 'w') as out:
    out.write('gene\tpresent_in_AHPND+\tabsent_in_AHPND+\tpresent_in_AHPND-\tabsent_in_AHPND-\todds_ratio\tp_value\n')
    for (p, gene, a, b, c, d, odds) in significant:
        out.write(f'{gene}\t{a}\t{b}\t{c}\t{d}\t{odds:.3f}\t{p:.6e}\n')

print(f"Análisis completado. {len(all_results)} genes evaluados.")
print(f"Genes significativos (p < 0.05): {len(significant)}")
print("Tablas guardadas:")
print("  - fisher_results.tsv (todos)")
print("  - significant_fisher_results.tsv (solo p<0.05)")
