#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extrae genes candidatos asociados a AHPND+ a partir de los resultados de Fisher.
Usa feature_sequences/<familia>.aa.fasta como fuente de secuencias aminoacídicas.
Genera:
  - candidate_fdr4.fa           -> secuencias de los 4 genes con FDR < 0.05
  - candidate_extreme113.fa     -> secuencias de los 113 genes con p<0.05 y odds ratio extremo (0 o inf)
  - candidate_fdr4.tsv          -> tabla resumen de los 4 genes
  - candidate_extreme113.tsv    -> tabla resumen de los 113 genes

Los archivos .fa contienen una secuencia representativa (la primera del archivo)
por familia génica, lista para usar en eggNOG-mapper (Galaxy).
"""

import csv
import os

# ---------- Parámetros ----------
FISHER_ADJ = "fisher_results_direct_adjusted.tsv"   # resultados con FDR y Bonferroni
PIRATE_TSV = "PIRATE.gene_families.ordered.tsv"     # tabla de PIRATE con productos
FEATURE_DIR = "feature_sequences"                   # carpeta con secuencias por familia

# ---------- Leer anotación de productos desde PIRATE ----------
print(f"Leyendo anotaciones desde {PIRATE_TSV}...")
gene_product = {}
with open(PIRATE_TSV) as f:
    header = f.readline().rstrip('\n').split('\t')
    for line in f:
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 5:
            continue
        fam_id = cols[1]
        product = cols[4] if cols[4] else "hypothetical protein"
        gene_product[fam_id] = product

print(f"  Productos cargados: {len(gene_product)}")

# ---------- Leer resultados ajustados ----------
print(f"Leyendo {FISHER_ADJ}...")
rows = []
with open(FISHER_ADJ) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        rows.append(row)

# ---------- Seleccionar genes candidatos ----------
fdr_set = set()
extreme_set = set()

for row in rows:
    gene = row['gene']
    try:
        p = float(row['p_value'])
    except:
        p = 1.0
    try:
        q = float(row['q_value_FDR'])
    except:
        q = 1.0
    odds = row['odds_ratio'].strip()

    if q < 0.05:
        fdr_set.add(gene)

    is_inf = odds == 'inf'
    is_zero = False
    try:
        # Si es numérico, comprobar si es 0.0
        odds_float = float(odds)
        if odds_float == 0.0:
            is_zero = True
    except:
        pass

    if p < 0.05 and (is_inf or is_zero):
        extreme_set.add(gene)

candidate_fdr = sorted(fdr_set)
candidate_extreme = sorted(extreme_set)

print(f"Genes FDR < 0.05: {len(candidate_fdr)}")
print(f"Genes odds ratio extremo (p<0.05): {len(candidate_extreme)}")

# ---------- Función para determinar patrón ----------
def get_pattern(row):
    try:
        a = int(row['present_in_AHPND+'])
        b = int(row['absent_in_AHPND+'])
        c = int(row['present_in_AHPND-'])
        d = int(row['absent_in_AHPND-'])
    except:
        return "Other"
    if a >= 5 and c <= 1:
        return "Gained in AHPND+"
    elif a <= 1 and c >= 3:
        return "Lost in AHPND+"
    else:
        return "Other"

# ---------- Extraer secuencia desde feature_sequences ----------
def get_sequence(gene):
    """Lee la primera secuencia del archivo gXXXXX.aa.fasta en FEATURE_DIR."""
    fasta_path = os.path.join(FEATURE_DIR, f"{gene}.aa.fasta")
    if not os.path.exists(fasta_path):
        return None
    with open(fasta_path) as f:
        seq = []
        header = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is None:
                    header = line[1:]
                elif seq:
                    # Si ya tenemos una secuencia, devolvemos la primera
                    break
            elif header is not None:
                seq.append(line)
        return ''.join(seq) if seq else None

# ---------- Guardar tablas y FASTA ----------
def write_candidates(gene_list, fasta_name, tsv_name):
    fasta_out = open(fasta_name, 'w')
    tsv_out = open(tsv_name, 'w', newline='')
    writer = csv.writer(tsv_out, delimiter='\t')
    writer.writerow(['gene', 'product', 'present_AHPND+', 'absent_AHPND+',
                     'present_AHPND-', 'absent_AHPND-', 'odds_ratio', 'p_value',
                     'q_value_FDR', 'pattern'])

    for gene in gene_list:
        row = next((r for r in rows if r['gene'] == gene), None)
        if row is None:
            print(f"  Advertencia: {gene} no está en los resultados")
            continue

        product = gene_product.get(gene, "hypothetical protein")
        seq = get_sequence(gene)
        if seq:
            # Escribir en FASTA con cabecera informativa
            fasta_out.write(f">{gene} {product}\n{seq}\n")
        else:
            print(f"  Advertencia: {gene} no tiene secuencia en {FEATURE_DIR}/{gene}.aa.fasta")

        writer.writerow([
            gene,
            product,
            row['present_in_AHPND+'],
            row['absent_in_AHPND+'],
            row['present_in_AHPND-'],
            row['absent_in_AHPND-'],
            row['odds_ratio'],
            row['p_value'],
            row['q_value_FDR'],
            get_pattern(row)
        ])

    fasta_out.close()
    tsv_out.close()
    print(f"  Guardado {fasta_name} y {tsv_name}")

# ---------- Generar archivos ----------
write_candidates(candidate_fdr, "candidate_fdr4.fa", "candidate_fdr4.tsv")
write_candidates(candidate_extreme, "candidate_extreme113.fa", "candidate_extreme113.tsv")

print("\nListo. Sube los archivos .fa a Galaxy y usa eggNOG-mapper para anotarlos.")
