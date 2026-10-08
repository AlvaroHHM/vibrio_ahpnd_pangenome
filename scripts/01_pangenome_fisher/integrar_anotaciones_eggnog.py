#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integra anotaciones de eggNOG (annotation y seed_orthologs) a las tablas de candidatos.
Uso:
    python3 integrar_anotaciones_eggnog.py \
        <candidatos.tsv> \
        <eggNOG_annotations.tabular> \
        <eggNOG_seed_orthologs.tabular> \
        <output.tsv>
"""

import csv
import sys

def read_annotations(annot_file):
    """Lee el archivo de anotaciones de eggNOG y devuelve dict gene -> fila con columnas relevantes."""
    annot = {}
    with open(annot_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row.get('#query', row.get('query', '')).strip()
            if not gene:
                continue
            # Seleccionar columnas de interés
            entry = {
                'eggNOG_Description': row.get('Description', ''),
                'Preferred_name': row.get('Preferred_name', ''),
                'COG_category': row.get('COG_category', ''),
                'GOs': row.get('GOs', ''),
                'EC': row.get('EC', ''),
                'KEGG_ko': row.get('KEGG_ko', ''),
                'PFAMs': row.get('PFAMs', ''),
            }
            annot[gene] = entry
    return annot

def read_orthologs(ortho_file):
    """Lee el archivo de ortólogos (seed_orthologs) y devuelve dict gene -> mejor hit."""
    ortho = {}
    with open(ortho_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row.get('#qseqid', row.get('qseqid', '')).strip()
            if not gene:
                continue
            entry = {
                'ortholog_sseqid': row.get('sseqid', ''),
                'ortholog_pident': row.get('pident', ''),
                'ortholog_qcov': row.get('qcov', ''),
                'ortholog_evalue': row.get('evalue', ''),
                'ortholog_bitscore': row.get('bitscore', ''),
            }
            # Si hay varios hits, conservar el de menor e-value (asumimos que ya vienen ordenados)
            if gene not in ortho:
                ortho[gene] = entry
    return ortho

def main():
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)
    
    candidatos_file = sys.argv[1]
    annot_file = sys.argv[2]
    ortho_file = sys.argv[3]
    output_file = sys.argv[4]

    print(f"Leyendo candidatos desde {candidatos_file}...")
    # Leer tabla de candidatos
    with open(candidatos_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        candidates = list(reader)

    print(f"  {len(candidates)} genes candidatos cargados.")

    print(f"Leyendo anotaciones de eggNOG desde {annot_file}...")
    annot = read_annotations(annot_file)
    print(f"  {len(annot)} genes anotados.")

    print(f"Leyendo ortólogos desde {ortho_file}...")
    ortho = read_orthologs(ortho_file)
    print(f"  {len(ortho)} genes con ortólogo.")

    # Definir las nuevas columnas
    new_cols = [
        'eggNOG_Description', 'Preferred_name', 'COG_category', 'GOs', 'EC', 'KEGG_ko', 'PFAMs',
        'ortholog_sseqid', 'ortholog_pident', 'ortholog_qcov', 'ortholog_evalue', 'ortholog_bitscore'
    ]

    # Cabecera original + nuevas columnas
    if candidates:
        fieldnames = list(candidates[0].keys()) + new_cols
    else:
        fieldnames = new_cols

    with open(output_file, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in candidates:
            gene = row.get('gene', '').strip()
            # Añadir anotación y ortólogo
            row.update(annot.get(gene, {k: '' for k in new_cols[:7]}))
            row.update(ortho.get(gene, {k: '' for k in new_cols[7:]}))
            writer.writerow(row)

    print(f"Tabla final guardada en {output_file}")

if __name__ == "__main__":
    main()
