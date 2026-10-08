#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Añade correcciones por comparaciones múltiples a fisher_results_direct.tsv.
Calcula:
  - p_bonferroni
  - q_value (FDR Benjamini-Hochberg)
Guarda fisher_results_direct_adjusted.tsv
"""

import csv

# ---------- Leer resultados ----------
rows = []
with open('fisher_results_direct.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['p_value'] in ('', 'NA', 'nan'):
            p = 1.0
        else:
            p = float(row['p_value'])
        row['p_value'] = p
        rows.append(row)

# ---------- Ordenar por p-valor ----------
rows_sorted = sorted(rows, key=lambda x: x['p_value'])
n_tests = len(rows_sorted)

# ---------- Calcular FDR (BH) ----------
q_values = [0.0] * n_tests
prev_q = 0.0
for i in range(n_tests - 1, -1, -1):   # desde el último (p mayor) hacia el primero
    p = rows_sorted[i]['p_value']
    rank = i + 1
    q = min(p * n_tests / rank, prev_q) if i < n_tests - 1 else p * n_tests / rank
    q_values[i] = q
    prev_q = q

for i, row in enumerate(rows_sorted):
    row['p_bonferroni'] = min(row['p_value'] * n_tests, 1.0)
    row['q_value_FDR'] = q_values[i]

# ---------- Guardar con todas las columnas originales + nuevas ----------
fieldnames = list(rows_sorted[0].keys())

with open('fisher_results_direct_adjusted.tsv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(rows_sorted)

print(f"Se guardó fisher_results_direct_adjusted.tsv con {n_tests} pruebas.")
print("Recuerda que el umbral de Bonferroni es:", 0.05 / n_tests)
