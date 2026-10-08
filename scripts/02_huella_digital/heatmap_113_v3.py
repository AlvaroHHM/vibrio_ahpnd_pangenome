#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Heatmap de los 113 genes candidatos usando PIRATE.gene_families.ordered.tsv como fuente.
- Extrae la matriz de presencia/ausencia para los genes seleccionados.
- Ordena los genes por categoría funcional y p-value.
- Ordena las cepas por AHPND+ primero, luego AHPND−; dentro de cada grupo por clado.
- Muestra solo la barra de fenotipo AHPND en la parte superior.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

# ---------- Parámetros ----------
CANDIDATES_FILE = "candidate_extreme113_annotated.tsv"
PIRATE_TSV = "PIRATE.gene_families.ordered.tsv"
CLADES_FILE = "clades_grupos.tsv"
PHENO_FILE = "phenotype.csv"

# Categorías funcionales
CATEGORY_RULES = [
    ("Conjugation/Secretion", ["conjug", "trb", "virb", "virD", "traG", "traD", "type II", "type IV", "secretion", "pilin", "N_methyl", "T2SS", "T4SS"]),
    ("Toxins/Virulence", ["toxin", "PirA", "PirB", "endotoxin", "pyocin", "zonular", "hok", "gef", "virulence"]),
    ("Regulation", ["transcriptional regulator", "regulator", "Fis", "LysR", "HTH", "cold shock"]),
    ("Mobile elements/DNA modification", ["transposase", "integrase", "resolvase", "methylase", "phage", "recombinase", "excisionase"]),
    ("Metabolism/Other", ["chitinase", "arylsulfatase", "metabolism", "ATPase", "hydrolase", "synthetase", "oxidoreductase", "transferase"]),
    ("Hypothetical/Unknown", ["hypothetical", "unknown", "DUF", "protein of unknown"])
]

def classify_gene(row):
    text = " ".join(str(x) for x in [row.get("eggNOG_Description", ""), row.get("Preferred_name", ""), row.get("product", "")]).lower()
    for category, keywords in CATEGORY_RULES:
        for kw in keywords:
            if kw.lower() in text:
                return category
    return "Hypothetical/Unknown"

# ---------- Leer candidatos ----------
print("Leyendo candidatos...")
cand = pd.read_csv(CANDIDATES_FILE, sep='\t', dtype=str)
cand['category'] = cand.apply(classify_gene, axis=1)
cand['p_value'] = pd.to_numeric(cand['p_value'], errors='coerce')
genes_ids = cand['gene'].tolist()
print(f"  {len(genes_ids)} genes candidatos")

# ---------- Leer fenotipos y clados ----------
print("Leyendo fenotipos...")
pheno = pd.read_csv(PHENO_FILE, index_col=0, dtype=str)
pheno.columns = ['AHPND']

print("Leyendo clados...")
clades = pd.read_csv(CLADES_FILE, sep='\t', index_col=0)

# ---------- Extraer matriz desde PIRATE TSV ----------
print(f"Extrayendo matriz desde {PIRATE_TSV}...")
with open(PIRATE_TSV) as f:
    header = f.readline().rstrip('\n').split('\t')
    # Detectar índices de cepas (las que están en pheno)
    strain_cols = []
    strain_names = []
    for i, col in enumerate(header):
        if col in pheno.index:
            strain_cols.append(i)
            strain_names.append(col)
    print(f"  Cepas detectadas: {len(strain_names)}")

    gene_presence = {gene: {strain: 0 for strain in strain_names} for gene in genes_ids}

    for line in f:
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 2:
            continue
        fam_id = cols[1]
        if fam_id not in gene_presence:
            continue
        for idx, strain in zip(strain_cols, strain_names):
            if idx < len(cols) and cols[idx].strip() != '':
                gene_presence[fam_id][strain] = 1

# Construir DataFrame genes x cepas
matrix = pd.DataFrame.from_dict(gene_presence, orient='index')
matrix = matrix.loc[genes_ids]

# Verificar que no esté vacío
print(f"  Tamaño de la matriz: {matrix.shape}")
print(f"  Total de 1s en la matriz: {matrix.values.sum()}")
print(f"  Primeras filas:\n{matrix.head().to_string()}")

# ---------- Ordenar genes por categoría y p-value ----------
category_order = ["Conjugation/Secretion", "Toxins/Virulence", "Regulation", "Mobile elements/DNA modification", "Metabolism/Other", "Hypothetical/Unknown"]
cand_sorted = cand.sort_values(['category', 'p_value'], key=lambda col: col if col.name != 'category' else col.map(lambda x: category_order.index(x) if x in category_order else 99))
gene_order = [g for g in cand_sorted['gene'] if g in matrix.index]
matrix_subset = matrix.loc[gene_order]

# ---------- Ordenar cepas ----------
strains = matrix.columns.tolist()
strain_info = pd.DataFrame(index=strains)
strain_info['AHPND'] = strain_info.index.map(lambda s: pheno.loc[s, 'AHPND'] if s in pheno.index else '0')
strain_info['AHPND'] = strain_info['AHPND'].map({'1': 1, '0': 0})
strain_info['clade'] = strain_info.index.map(lambda s: clades.loc[s, 'clade'] if s in clades.index else 'NA')
strain_info_sorted = strain_info.sort_values(['AHPND', 'clade'], ascending=[False, True])
ordered_strains = strain_info_sorted.index.tolist()
matrix_subset = matrix_subset[ordered_strains]

# ---------- Anotaciones ----------
category_colors = {
    "Conjugation/Secretion": "#3B4992",
    "Toxins/Virulence": "#EE0000",
    "Regulation": "#008B45",
    "Mobile elements/DNA modification": "#631879",
    "Metabolism/Other": "#008280",
    "Hypothetical/Unknown": "#BB0021"
}
row_colors = cand_sorted.set_index('gene')['category'].map(category_colors)
row_colors = row_colors.reindex(gene_order)

ahpnd_values = strain_info_sorted['AHPND'].map({1: 1, 0: 0}).values
ahpnd_cmap = ListedColormap(['#BBBBBB', '#EE0000'])  # 0 = gris, 1 = rojo

# ---------- Heatmap ----------
print("Generando heatmap...")
fig, ax = plt.subplots(figsize=(14, len(gene_order)*0.15 + 2))
sns.heatmap(matrix_subset, cmap="Blues", cbar=False, linewidths=0.2, linecolor='gray', ax=ax)

ax.set_xticks(np.arange(len(ordered_strains)))
ax.set_yticks(np.arange(len(gene_order)))
ax.set_xticklabels(ordered_strains, rotation=45, ha='right', fontsize=8)
ax.set_yticklabels(gene_order, fontsize=7)
ax.set_ylabel("")

for label in ax.get_yticklabels():
    gene = label.get_text()
    if gene in row_colors.index:
        label.set_color(row_colors.loc[gene])
        label.set_fontweight('bold')

# Solo barra AHPND superior
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
ax_ahpnd = divider.append_axes("top", size="2%", pad=0.1)

ax_ahpnd.imshow(np.array([ahpnd_values]), aspect='auto', cmap=ahpnd_cmap, interpolation='nearest')
ax_ahpnd.set_xticks(range(len(ordered_strains)))
ax_ahpnd.set_xticklabels([])
ax_ahpnd.set_yticks([])
ax_ahpnd.set_title("AHPND", fontsize=8)

plt.tight_layout()
plt.savefig("heatmap_113.png", dpi=300, bbox_inches='tight')
plt.savefig("heatmap_113.svg", dpi=300, bbox_inches='tight')
print("Heatmap guardado como heatmap_113.png y heatmap_113.svg")
