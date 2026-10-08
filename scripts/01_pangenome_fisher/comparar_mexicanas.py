#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compara genes entre las cepas mexicanas:
11_VM, 15_CESAIBC, 6_VM, T13P
"""

import csv

# Cepas de interés
mex = ["11_VM", "15_CESAIBC", "6_VM", "T13P"]

# Leer TSV y localizar columnas
with open('PIRATE.gene_families.ordered.tsv') as f:
    header = f.readline().rstrip('\n').split('\t')
    col_idx = {strain: header.index(strain) for strain in mex if strain in header}
    missing = [s for s in mex if s not in col_idx]
    if missing:
        print(f"Faltan cepas en el TSV: {missing}")
        exit(1)

    # Almacenar resultados
    gene_info = []   # (family_id, product, dict_presence)
    for line in f:
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 2:
            continue
        fam_id = cols[1]
        product = cols[4] if len(cols) > 4 else ''
        pres = {}
        for strain, idx in col_idx.items():
            pres[strain] = 1 if idx < len(cols) and cols[idx].strip() != '' else 0
        gene_info.append((fam_id, product, pres))

# Conjuntos
shared_15_T13P = []       # presentes en 15 y T13P, ausentes en 11 y 6
exclusive_15 = []         # presente solo en 15
exclusive_T13P = []       # presente solo en T13P
lost_15_present_11_6 = [] # ausente en 15, presente en 11 y 6
lost_15_present_all3 = [] # ausente en 15, presente en 11, 6 y T13P

for fam, prod, pres in gene_info:
    p11 = pres["11_VM"]; p15 = pres["15_CESAIBC"]; p6 = pres["6_VM"]; pT = pres["T13P"]

    if p15 == 1 and pT == 1 and p11 == 0 and p6 == 0:
        shared_15_T13P.append((fam, prod))
    if p15 == 1 and p11 == 0 and p6 == 0 and pT == 0:
        exclusive_15.append((fam, prod))
    if pT == 1 and p11 == 0 and p6 == 0 and p15 == 0:
        exclusive_T13P.append((fam, prod))
    if p15 == 0 and p11 == 1 and p6 == 1:
        lost_15_present_11_6.append((fam, prod))
    if p15 == 0 and p11 == 1 and p6 == 1 and pT == 1:
        lost_15_present_all3.append((fam, prod))

# Escribir resultados
def write_list(filename, data):
    with open(filename, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["gene_family", "product"])
        writer.writerows(data)
    print(f"{filename}: {len(data)} genes")

write_list("shared_15_T13P.tsv", shared_15_T13P)
write_list("exclusive_15.tsv", exclusive_15)
write_list("exclusive_T13P.tsv", exclusive_T13P)
write_list("lost_15_present_11_6.tsv", lost_15_present_11_6)
write_list("lost_15_present_all3.tsv", lost_15_present_all3)

print("\nResumen:")
print(f"Genes compartidos 15+T13P y ausentes en 11/6: {len(shared_15_T13P)}")
print(f"Genes exclusivos de 15: {len(exclusive_15)}")
print(f"Genes exclusivos de T13P: {len(exclusive_T13P)}")
print(f"Genes ausentes en 15 y presentes en 11/6: {len(lost_15_present_11_6)}")
print(f"Genes ausentes en 15 y presentes en 11/6/T13P: {len(lost_15_present_all3)}")
