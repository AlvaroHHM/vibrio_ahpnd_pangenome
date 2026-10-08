#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure S1: COG V (defense) gene heatmap across the three Mexican strains.

Uses Roary outputs (3 Mexican strains only).
Inputs:
  - input_data/pangenome/roary_3mexicanas/Supplementary_Table_exclusivos_ausentes.tsv
  - input_data/pangenome/roary_3mexicanas/gene_presence_absence.csv

Outputs:
  - outputs/supplementary_figures/Figure_S1_COG_V_heatmap.pdf
  - outputs/supplementary_figures/Figure_S1_COG_V_heatmap.png
  - outputs/huella/genes_COG_V.tsv (data table)
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# --- Rutas relativas al repo ---
INPUT_DIR = "input_data/pangenome/roary_3mexicanas"
OUT_FIG_DIR = "outputs/supplementary_figures"
OUT_TABLE_DIR = "outputs/huella"

os.makedirs(OUT_FIG_DIR, exist_ok=True)
os.makedirs(OUT_TABLE_DIR, exist_ok=True)

# --- Cargar tabla fusionada y filtrar COG V ---
df = pd.read_csv(os.path.join(INPUT_DIR, "Supplementary_Table_exclusivos_ausentes.tsv"),
                 sep='\t')

mask = df["COG_category"].fillna("").str.contains("V", case=False, regex=False)
genes_v = df[mask].copy()
print(f"Genes con COG V: {len(genes_v)}")

# --- Cargar Roary original ---
roary = pd.read_csv(os.path.join(INPUT_DIR, "gene_presence_absence.csv"))

def presence(value):
    if pd.isna(value) or str(value).strip() == "":
        return 0
    return 1

presence_rows = []
for _, row in genes_v.iterrows():
    gene = row["gene"]
    set_name = row["Set"]
    annotation = row.get("annotation_roary", row.get("annotation", ""))
    match = roary[roary["Gene"] == gene]
    if match.empty:
        print(f"Advertencia: {gene} no encontrado en Roary")
        continue
    r = match.iloc[0]
    presence_rows.append({
        "gene": gene,
        "annotation": annotation,
        "Set": set_name,
        "11_VM": presence(r["11_VM"]),
        "15_CESAIBC": presence(r["15_CESAIBC"]),
        "6_VM": presence(r["6_VM"]),
    })

presence_df = pd.DataFrame(presence_rows)
presence_df.to_csv(os.path.join(OUT_TABLE_DIR, "genes_COG_V.tsv"),
                   sep="\t", index=False)
print("Tabla guardada: outputs/huella/genes_COG_V.tsv")

# --- Preparar heatmap ---
order_df = presence_df.sort_values(["Set", "annotation"])
order = order_df["gene"].tolist()

x_labels = order_df.set_index("gene").apply(
    lambda x: f"{x.name}\n{x['annotation']}", axis=1
).reindex(order).tolist()

heatmap_data = order_df.set_index("gene")[["11_VM", "15_CESAIBC", "6_VM"]].T[order]

set_colors = {"Exclusive_15": "#3B4992", "Absent_15": "#EE0000"}
set_category = order_df.set_index("gene")["Set"].reindex(order)

# --- Graficar (sin seaborn) ---
fig, ax = plt.subplots(figsize=(max(12, len(order) * 0.45), 7))
im = ax.imshow(heatmap_data.values, cmap="Blues", vmin=0, vmax=1, aspect='auto')

# Anotaciones numéricas en cada celda
for i in range(heatmap_data.shape[0]):
    for j in range(heatmap_data.shape[1]):
        val = heatmap_data.values[i, j]
        color = 'white' if val > 0.5 else 'black'
        ax.text(j, i, str(int(val)), ha='center', va='center',
                fontsize=9, color=color)

# Etiquetas X
ax.set_xticks(np.arange(len(order)))
ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=7)
for xtick, gene in zip(ax.get_xticklabels(), order):
    xtick.set_color(set_colors[set_category[gene]])
    xtick.set_fontweight('bold')

# Eje Y
ax.set_yticks(np.arange(heatmap_data.shape[0]))
ax.set_yticklabels(["11_VM", "15_CESAIBC", "6_VM"], rotation=0, fontsize=12)
ax.set_ylabel("Strain", fontsize=12)
ax.set_title("COG V (defense) genes in the three Mexican strains", fontsize=14)

# Grid sutil
ax.set_xticks(np.arange(-0.5, len(order), 1), minor=True)
ax.set_yticks(np.arange(-0.5, heatmap_data.shape[0], 1), minor=True)
ax.grid(which='minor', color='gray', linewidth=0.5)
ax.tick_params(which='minor', length=0)

# Leyenda
legend_elements = [
    Line2D([0], [0], color=set_colors["Exclusive_15"], lw=4,
           label='Exclusive to 15_CESAIBC'),
    Line2D([0], [0], color=set_colors["Absent_15"], lw=4,
           label='Absent from 15_CESAIBC'),
]
ax.legend(handles=legend_elements, bbox_to_anchor=(1.02, 1),
          loc='upper left', borderaxespad=0.)

plt.tight_layout()
out_pdf = os.path.join(OUT_FIG_DIR, "Figure_S1_COG_V_heatmap.pdf")
out_png = os.path.join(OUT_FIG_DIR, "Figure_S1_COG_V_heatmap.png")
plt.savefig(out_png, dpi=300, format='png', bbox_inches='tight', facecolor='white')
plt.savefig(out_pdf, dpi=300, format='pdf', bbox_inches='tight', facecolor='white')
plt.close()

print(f"Heatmap guardado: {out_pdf}")
print(f"Heatmap guardado: {out_png}")
