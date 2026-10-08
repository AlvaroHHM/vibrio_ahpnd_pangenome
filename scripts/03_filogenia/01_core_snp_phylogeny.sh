#!/bin/bash
# ============================================================
# 05_core_snp_phylogeny.sh
# Filtrado de recombinación (Gubbins 2.4.1) e inferencia filogenética (IQ-TREE)
# para el core-genome alineado de 19 cepas de Vibrio parahaemolyticus.
#
# Uso:
#   bash 05_core_snp_phylogeny.sh
#
# Entradas:
#   ../core_alignment.fasta  (desde PIRATE_results_v3)
#
# Salidas:
#   gubbins_out.*  (archivos generados con prefijo)
#   iqtree_out/    (filogenia IQ-TREE)
# ============================================================

set -e

# ---------- Configuración de rutas ----------
PIRATE_DIR=~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3
WORKDIR="$PIRATE_DIR/core_snp_phylogeny"
ALN_IN="$PIRATE_DIR/core_alignment.fasta"

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# ---------- 1. Localizar ejecutable de Gubbins ----------
GUBBINS_BIN=""
if command -v run_gubbins.py &> /dev/null; then
    GUBBINS_BIN=$(command -v run_gubbins.py)
elif command -v gubbins &> /dev/null; then
    GUBBINS_BIN=$(command -v gubbins)
else
    CONDA_ENV_BIN="$HOME/anaconda3/envs/bacterial-genomics-tutorial/bin"
    if [ -x "$CONDA_ENV_BIN/run_gubbins.py" ]; then
        GUBBINS_BIN="$CONDA_ENV_BIN/run_gubbins.py"
    elif [ -x "$CONDA_ENV_BIN/gubbins" ]; then
        GUBBINS_BIN="$CONDA_ENV_BIN/gubbins"
    fi
fi

if [ -z "$GUBBINS_BIN" ]; then
    echo "ERROR: No se encontró el ejecutable de Gubbins."
    exit 1
fi

echo "Usando Gubbins: $GUBBINS_BIN"
"$GUBBINS_BIN" --version 2>&1 | head -n 2 || true

# ---------- 2. Filtrar recombinación con Gubbins ----------
echo "Ejecutando Gubbins ..."
if [ -f "gubbins_out.filtered_polymorphic_sites.fasta" ]; then
    echo "El alineamiento filtrado ya existe."
else
    "$GUBBINS_BIN" \
        --tree_builder raxml \
        --prefix gubbins_out \
        --threads 4 \
        "$ALN_IN"
fi

# ---------- 3. Buscar el alineamiento filtrado ----------
FILTERED_ALN="gubbins_out.filtered_polymorphic_sites.fasta"
if [ ! -f "$FILTERED_ALN" ]; then
    # Buscar variante .fa
    FILTERED_ALN=$(find . -maxdepth 1 -name "gubbins_out.filtered_polymorphic_sites.fa" | head -n 1)
fi

if [ -z "$FILTERED_ALN" ]; then
    echo "ERROR: No se encontró el alineamiento filtrado."
    exit 1
fi

echo "Alineamiento filtrado: $FILTERED_ALN"

# ---------- 4. Inferir filogenia con IQ-TREE ----------
echo "Ejecutando IQ-TREE ..."
mkdir -p iqtree_out
iqtree \
    -s "$FILTERED_ALN" \
    -m GTR+F+I+G4 \
    -bb 1000 \
    -alrt 1000 \
    -nt AUTO \
    -pre iqtree_out/core_snp

echo "Resultados en:"
echo "  gubbins_out.*"
echo "  iqtree_out/"
