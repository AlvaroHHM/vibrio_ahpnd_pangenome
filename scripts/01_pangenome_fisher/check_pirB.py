#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verifica qué cepas poseen el gen pirB (g05919) usando PIRATE.gene_families.ordered.tsv.
Imprime la lista de cepas con presencia/ausencia y guarda un archivo TSV.
"""

import csv

# Parámetros
GENE_ID = "g05919"          # ID de familia génica de pirB
TSV_FILE = "PIRATE.gene_families.ordered.tsv"
OUTPUT_FILE = "presencia_pirB.tsv"

# Leer cabecera y localizar columnas de cepas
with open(TSV_FILE) as f:
    header = f.readline().rstrip('\n').split('\t')
    # Las columnas de cepas empiezan después de las primeras 21 columnas de metadatos
    # Pero mejor detectar automáticamente: buscar columnas cuyo nombre no es vacío y no son metadatos conocidos.
    # Sabemos que las cepas están al final y contienen un locus_tag en alguna fila.
    # Usaremos el criterio de que después de la columna 21 (índice 20) están las cepas.
    # Ajuste: en versiones anteriores vimos que las cepas empiezan en la columna 22 (índice 21).
    # Verificaremos con el nombre "11_VM" o "T13P".
    strain_cols = []
    strain_names = []
    for i, col in enumerate(header):
        if col in ("11_VM", "15_CESAIBC", "6_VM", "T13P", "13028_A3", "AG1"):
            # Encontramos una columna de cepa, asumimos que desde aquí todas son cepas
            strain_cols = list(range(i, len(header)))
            strain_names = header[i:]
            break
    
    if not strain_cols:
        # Fallback: usar índices fijos 22 en adelante
        strain_cols = list(range(22, len(header)))
        strain_names = header[22:]
        print("No se detectaron columnas de cepas automáticamente, usando desde columna 22")
    
    print(f"Columnas de cepas detectadas: {len(strain_names)}")

    # Buscar la fila de g05919
    gene_row = None
    for line in f:
        cols = line.rstrip('\n').split('\t')
        if len(cols) > 1 and cols[1] == GENE_ID:
            gene_row = cols
            break

if gene_row is None:
    print(f"No se encontró la familia génica {GENE_ID}")
else:
    results = []
    print(f"\nPresencia de {GENE_ID} (pirB):")
    for i, strain in zip(strain_cols, strain_names):
        present = False
        if i < len(gene_row) and gene_row[i].strip() != '':
            present = True
        results.append((strain, present))
        status = "PRESENTE" if present else "ausente"
        print(f"  {strain:20s} {status}")
    
    # Guardar archivo
    with open(OUTPUT_FILE, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["cepa", "pirB_presente"])
        for strain, present in results:
            writer.writerow([strain, 1 if present else 0])
    
    print(f"\nResultados guardados en {OUTPUT_FILE}")
