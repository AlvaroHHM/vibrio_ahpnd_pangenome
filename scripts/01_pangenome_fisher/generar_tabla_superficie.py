#!/usr/bin/env python3
import csv, os

# ============================================================
# Datos fijos: los 7 genes de superficie (locus tag real y producto)
# ============================================================
genes = [
    ("pilin",                                  "JNOILPDC_01230"),
    ("T3SS needle (EscF/YscF)",                "JNOILPDC_02661"),
    ("maltoporin",                             "JNOILPDC_03102"),
    ("BamA/TamA outer membrane protein",       "JNOILPDC_03135"),
    ("MshA",                                   "JNOILPDC_03761"),
    ("MshP",                                   "JNOILPDC_03765"),
    ("HlyD secretion protein",                 "JNOILPDC_03924"),
]

# ============================================================
# 1. Extraer ortólogos e identidades desde los archivos BLAST
# ============================================================
blast_files = {
    "6VM": "blast_6VM.txt",
    "11VM": "blast_11VM.txt",
}
ortologos = {}
identidades = {}

for cepa, fname in blast_files.items():
    ortologos[cepa] = {}
    identidades[cepa] = {}
    with open(fname) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                qseqid = parts[0]
                sseqid = parts[1]
                pident = parts[2]
                # Guardamos el primer hit (mejor) para cada query
                if qseqid not in ortologos[cepa]:
                    ortologos[cepa][qseqid] = sseqid
                    identidades[cepa][qseqid] = pident

# ============================================================
# 2. Obtener contig desde el archivo GFF de 15_CESAIBC
# ============================================================
gff_path = os.path.expanduser(
    "~/bacterial-genomics-tutorial/reference_genomes_pipeline/"
    "prokka_annotations_pirAB/15_CESAIBC/15_CESAIBC.gff"
)
contig_map = {}
with open(gff_path) as gff:
    for line in gff:
        if line.startswith('#') or not line.strip():
            continue
        cols = line.split('\t')
        if len(cols) < 9:
            continue
        contig = cols[0]
        attr = cols[8]
        # Buscar locus_tag
        for part in attr.split(';'):
            if part.startswith('locus_tag='):
                locus = part.split('=')[1].strip()
                contig_map[locus] = contig
                break

# ============================================================
# 3. Generar la tabla CSV
# ============================================================
output = "tabla_genes_superficie.csv"
with open(output, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "producto", "locus_15CESAIBC", "contig_15",
        "ortologo_6VM", "ortologo_11VM",
        "identidad_6VM(%)", "identidad_11VM(%)"
    ])
    for producto, locus15 in genes:
        contig = contig_map.get(locus15, "no_encontrado")
        orto6 = ortologos["6VM"].get(locus15, "")
        orto11 = ortologos["11VM"].get(locus15, "")
        id6 = identidades["6VM"].get(locus15, "")
        id11 = identidades["11VM"].get(locus15, "")
        writer.writerow([producto, locus15, contig, orto6, orto11, id6, id11])

print(f"Tabla generada: {output}")
