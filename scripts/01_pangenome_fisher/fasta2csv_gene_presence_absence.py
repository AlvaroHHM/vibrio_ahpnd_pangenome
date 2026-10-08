#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convierte binary_presence_absence.fasta de PIRATE a CSV con 0/1.
Mapeo confirmado: A = ausente (0), C = presente (1).
"""

from Bio import SeqIO
import os

gene_names_file = 'gene_names.txt'

# Si el archivo de nombres no existe, se genera desde PIRATE.gene_families.ordered.tsv
if not os.path.isfile(gene_names_file):
    print("Generando gene_names.txt desde PIRATE.gene_families.ordered.tsv...")
    with open('PIRATE.gene_families.ordered.tsv') as f:
        header = f.readline()  # saltar cabecera
        gene_names = []
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) > 1 and cols[1]:
                gene_names.append(cols[1])
    with open(gene_names_file, 'w') as out:
        out.write('\n'.join(gene_names) + '\n')
    print(f"Se escribieron {len(gene_names)} genes en {gene_names_file}")

fasta_file = "binary_presence_absence.fasta"
gene_names_file = "gene_names.txt"
output_csv = "gene_presence_absence.csv"

# Leer nombres de genes
with open(gene_names_file) as f:
    gene_names = [line.strip() for line in f if line.strip()]

# Leer secuencias y convertir a 1/0
strains = []
presence_matrix = []
for record in SeqIO.parse(fasta_file, "fasta"):
    strains.append(record.id)
    seq = str(record.seq)
    if len(seq) != len(gene_names):
        print(f"WARNING: length mismatch for {record.id}: seq len {len(seq)} vs genes {len(gene_names)}")
    # Convertir caracteres a 1/0
    binary_seq = []
    for ch in seq:
        if ch in ('1', 'A', 'a'):
            binary_seq.append('0')   # A = ausente
        elif ch in ('0', 'C', 'c'):
            binary_seq.append('1')   # C = presente
        else:
            # Si aparece otro carácter, asumimos presente (1) por defecto
            print(f"  [aviso] carácter inesperado '{ch}' en {record.id}; se interpreta como 1")
            binary_seq.append('1')
    presence_matrix.append(''.join(binary_seq))

# Escribir CSV con 0/1
with open(output_csv, "w") as out:
    out.write("Strain," + ",".join(gene_names) + "\n")
    for strain, seq in zip(strains, presence_matrix):
        out.write(strain + "," + ",".join(seq) + "\n")

print(f"CSV file created: {output_csv} (con valores 0/1)")
