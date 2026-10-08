#!/bin/bash
set -e

mkdir -p raw_genomes
cd raw_genomes

echo "=== Descargando genomas con comandos corregidos ==="

# 1. L2171: usar acceso de ensamblaje completo (RefSeq)
wget -O L2171.fna.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/435/455/GCF_030435455.1_ASM30435455v1/GCF_030435455.1_ASM30435455v1_genomic.fna.gz"
gunzip L2171.fna.gz

# 2. L2181: usar acceso de ensamblaje completo
wget -O L2181.fna.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/435/535/GCF_030435535.1_ASM30435535v1/GCF_030435535.1_ASM30435535v1_genomic.fna.gz"
gunzip L2181.fna.gz

# 3. FORC_023 (ya se descargó correctamente, pero lo dejamos como está)
#    En el script original, ya tienes un comando que funciona: efetch -db nucleotide -id CP012950 -format fasta > FORC_023.fna

# 4. FORC_022: usar acceso de ensamblaje completo (RefSeq)
wget -O FORC_022.fna.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/886/195/GCF_001886195.1_ASM188619v1/GCF_001886195.1_ASM188619v1_genomic.fna.gz"
gunzip FORC_022.fna.gz

# 5. RIMD 2210633: usar acceso de ensamblaje (RefSeq)
wget -O RIMD2210633.fna.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/196/705/GCF_000196705.1_ASM19670v1/GCF_000196705.1_ASM19670v1_genomic.fna.gz"
gunzip RIMD2210633.fna.gz

# 6. Ba94C2: usar acceso del BioProject (GenBank)
wget -O Ba94C2.fna "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/717/385/GCA_001717385.1_ASM171738v1/GCA_001717385.1_ASM171738v1_genomic.fna.gz"
gunzip Ba94C2.fna.gz

# 7. M1-1: acceso de nucleótido (CP022157) pero es un solo cromosoma; se descarga bien con efetch
efetch -db nucleotide -id CP022157 -format fasta > M1-1.fna

# 8. CGVP3, CGVP8, CGVP22: los números de acceso que teníamos no son correctos. 
#    Si no encuentras los accesos, omítelos o búscalos manualmente en NCBI.

echo "=== Descarga completada. Archivos disponibles ==="
ls -lh *.fna
