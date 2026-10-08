#!/bin/bash
# download_genomes_v2.sh - Descarga genomas de referencia usando efetch y wget
# Ejecutar: bash download_genomes_v2.sh

set -e

mkdir -p raw_genomes
cd raw_genomes

echo "=== Descargando genomas de referencia ==="

# 1. L2171 (Sudamérica, AHPND+)
# Usamos los números de acceso de los cromosomas (CP176031-CP176034)
echo "Descargando L2171 (4 cromosomas/plásmidos)..."
efetch -db nuccore -id CP176031,CP176032,CP176033,CP176034 -format fasta > L2171.fna

# 2. L2181 (Sudamérica, AHPND+)
efetch -db nuccore -id CP176035,CP176036,CP176037 -format fasta > L2181.fna

# 3. CGVP3 (Corea, AHPND+)
# Número de acceso del ensamblaje (obtenido de la literatura)
efetch -db assembly -id GCF_009859285.1 -format fasta > CGVP3.fna

# 4. CGVP8 (Corea, AHPND+)
efetch -db assembly -id GCA_016913245.1 -format fasta > CGVP8.fna

# 5. CGVP22 (Corea, AHPND+)
efetch -db assembly -id GCA_016913245.1 -format fasta > CGVP22.fna  # Ajustar si es otro

# 6. M1-1 (Vietnam, AHPND+)
# Obtener el número de acceso del genoma desde el artículo
efetch -db nucleotide -id CP022157 -format fasta > M1-1.fna

# 7. Ba94C2 (Sudamérica, AHPND+)
# Usar el BioProject para obtener el genoma
efetch -db assembly -id GCA_001717385.1 -format fasta > Ba94C2.fna

# 8. 19-021-D1 (Corea, AHPND+)
efetch -db assembly -id GCF_031498635.1 -format fasta > 19-021-D1.fna

# 9. 13028/A3 (Vietnam, AHPND+)
efetch -db assembly -id GCF_000579255.1 -format fasta > 13028_A3.fna

# 10. FORC_022 (Corea, AHPND-)
efetch -db assembly -id GCF_001886195.1 -format fasta > FORC_022.fna

# 11. FORC_023 (Corea, AHPND-)
efetch -db nucleotide -id CP012950 -format fasta > FORC_023.fna

# 12. RIMD 2210633 (Japón, AHPND-)
efetch -db assembly -id GCF_000196705.1 -format fasta > RIMD2210633.fna

# 13. VP291 (México, AHPND?) - No encontrado
echo "Advertencia: No se encontró acceso para VP291. Se omite."

echo "=== Descarga completada. Genomas disponibles en raw_genomes/ ==="
ls -lh *.fna
