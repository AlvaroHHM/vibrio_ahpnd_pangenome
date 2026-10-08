#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Busca secuencias de genes de interés en las cepas mexicanas y las compara con T13P.
Los genes de interés incluyen:
  - vpadF (surface adhesin)
  - vscF / T3SS needle (variante de 15_CESAIBC y copia de 11/6)
  - hcp-2 (T6SS substrate)
  - pilin, MshA, MshP, etc. (opcional, si se encuentran)

El script:
  1. Localiza automáticamente los archivos .faa de Prokka en ~/bacterial-genomics-tutorial/
  2. Busca genes por locus_tag conocido o por nombre de producto en los .gff
  3. Extrae la secuencia proteica y la BLASTa contra T13P
  4. Genera tabla comparativa
"""

import os
import glob
import subprocess
import sys
from pathlib import Path

# ---------- Configuración ----------
T13P_DIR = None
CEPA_DIRS = {}
GFF_DIRS = []

# Intentar localizar las carpetas Prokka
for root, dirs, files in os.walk(os.path.expanduser("~/bacterial-genomics-tutorial")):
    if "T13P_prokka" in dirs:
        T13P_DIR = os.path.join(root, "T13P_prokka")
    for cepa in ["11_VM", "6_VM", "15_CESAIBC"]:
        if cepa in dirs and cepa not in CEPA_DIRS:
            CEPA_DIRS[cepa] = os.path.join(root, cepa)

print("Carpetas encontradas:")
for k,v in CEPA_DIRS.items():
    print(f"  {k}: {v}")
if T13P_DIR:
    print(f"  T13P: {T13P_DIR}")

# ---------- Genes de interés ----------
# Lista de locus_tags conocidos (puedes añadir más)
GENE_TARGETS = {
    "vpadF_11VM": ("11_VM", "11_VM_03384"),
    "vpadF_6VM":  ("6_VM",  "6_VM_03257"),
    "vscF_15":    ("15_CESAIBC", "JNOILPDC_02661"),   # Ajusta si no es correcto
    "vscF_11VM":  ("11_VM", "DKELKKMJ_03372"),        # Ajusta si no es correcto
    "vscF_6VM":   ("6_VM",  "CLIBEICH_03498"),        # Ajusta si no es correcto
    "hcp2_11VM":  ("11_VM", None),  # se buscará por producto
    "hcp2_6VM":   ("6_VM",  None),
}

# También buscar por nombre de producto en .gff
PRODUCT_SEARCH = {
    "hcp-2": "hcp",
    "surface adhesin": "surface adhesin",
    "T3SS needle": "needle",
}

def extract_sequence_from_faa(faa_path, locus_tag):
    """Extrae la secuencia FASTA de un locus_tag dado."""
    seq = ""
    with open(faa_path) as f:
        found = False
        for line in f:
            if line.startswith(">") and locus_tag in line:
                found = True
                seq += line
            elif line.startswith(">") and found:
                break
            elif found:
                seq += line
    return seq if seq else None

def run_blastp(query_file, db_path):
    """Ejecuta blastp y devuelve el mejor hit."""
    cmd = ["blastp", "-query", query_file, "-db", db_path,
           "-outfmt", "6 qseqid sseqid pident length qcovs evalue bitscore",
           "-evalue", "1e-5", "-max_target_seqs", "1"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout.strip():
        return result.stdout.strip().split("\n")[0]
    return None

def find_gene_by_product(gff_path, product_pattern):
    """Busca un gen por producto y devuelve su locus_tag."""
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            attr = cols[8]
            if product_pattern.lower() in attr.lower():
                # Extraer locus_tag
                for part in attr.split(";"):
                    if part.startswith("locus_tag="):
                        return part.split("=")[1].strip()
    return None

# ---------- Localizar archivos ----------
def find_faa(dirpath):
    return os.path.join(dirpath, "faa") if False else None  # placeholder

# Realmente buscaremos .faa dentro de cada carpeta
faa_paths = {}
for cepa, dirpath in CEPA_DIRS.items():
    faa = os.path.join(dirpath, f"{cepa}.faa")
    if os.path.exists(faa):
        faa_paths[cepa] = faa
    else:
        # buscar cualquier .faa
        for f in glob.glob(os.path.join(dirpath, "*.faa")):
            faa_paths[cepa] = f
            break

gff_paths = {}
for cepa, dirpath in CEPA_DIRS.items():
    gff = os.path.join(dirpath, f"{cepa}.gff")
    if os.path.exists(gff):
        gff_paths[cepa] = gff
    else:
        for f in glob.glob(os.path.join(dirpath, "*.gff")):
            gff_paths[cepa] = f
            break

# Base de datos de T13P
t13p_faa = None
if T13P_DIR:
    for f in glob.glob(os.path.join(T13P_DIR, "*.faa")):
        t13p_faa = f
        break
if not t13p_faa:
    print("No se encontró T13P.faa")
    sys.exit(1)

# Crear base de datos BLAST de T13P
db_name = "T13P_blastdb"
subprocess.run(["makeblastdb", "-in", t13p_faa, "-dbtype", "prot", "-out", db_name], capture_output=True)

# ---------- Procesar genes ----------
results = []
for gene_name, (cepa, locus_tag) in GENE_TARGETS.items():
    # Si locus_tag es None, buscar por producto
    if locus_tag is None:
        if cepa in gff_paths:
            # Probar con algunos términos
            for prod in ["hcp", "type VI secretion"]:
                locus_tag = find_gene_by_product(gff_paths[cepa], prod)
                if locus_tag:
                    break
    if not locus_tag:
        print(f"No se encontró locus_tag para {gene_name} en {cepa}, saltando.")
        continue
    # Extraer secuencia
    seq = extract_sequence_from_faa(faa_paths[cepa], locus_tag) if cepa in faa_paths else None
    if not seq:
        print(f"No se pudo extraer secuencia para {locus_tag} en {cepa}.")
        continue
    # Guardar temporal
    query_file = f"temp_{gene_name}.faa"
    with open(query_file, "w") as out:
        out.write(seq)
    # BLAST contra T13P
    blast_result = run_blastp(query_file, db_name)
    if blast_result:
        results.append((gene_name, cepa, locus_tag, blast_result))
    else:
        results.append((gene_name, cepa, locus_tag, "NO HIT"))
    os.remove(query_file)

# ---------- Mostrar resultados ----------
print("\nResultados BLAST contra T13P:")
print("="*100)
print(f"{'Gen':<15} {'Cepa':<12} {'Locus':<15} {'Mejor hit T13P':<20} {'%id':<6} {'Cobertura':<10} {'E-value':<10}")
for gene_name, cepa, locus_tag, blast in results:
    if blast != "NO HIT":
        parts = blast.split("\t")
        sseqid, pident, length, qcovs, evalue = parts[1], parts[2], parts[3], parts[4], parts[5]
        print(f"{gene_name:<15} {cepa:<12} {locus_tag:<15} {sseqid:<20} {pident:<6} {qcovs:<10} {evalue:<10}")
    else:
        print(f"{gene_name:<15} {cepa:<12} {locus_tag:<15} {'NO HIT':<20} {'-':<6} {'-':<10} {'-':<10}")

# Guardar tabla
with open("resultados_T13P.tsv", "w") as out:
    out.write("gen\tcepa\tlocus_tag\tblast_result\n")
    for gene_name, cepa, locus_tag, blast in results:
        out.write(f"{gene_name}\t{cepa}\t{locus_tag}\t{blast}\n")
print("\nResultados guardados en resultados_T13P.tsv")
