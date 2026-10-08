#!/bin/bash
# process_genomes.sh - Anota y ejecuta PIRATE sobre genomas en raw_genomes/
# Ejecutar: bash process_genomes.sh

set -e

# Configuración
THREADS=4
GENUS="Vibrio"
SPECIES="parahaemolyticus"

# Crear directorios de trabajo
mkdir -p prokka_annotations pirate_input final_results

# ----------------------------------------------------------------------
# PASO 1: Anotar todos los genomas de raw_genomes/ con Prokka
# ----------------------------------------------------------------------
echo "=== Anotando genomas con Prokka ==="
for fna in raw_genomes/*.fna; do
    base=$(basename "$fna" .fna)
    echo "Procesando $base ..."
    prokka --outdir prokka_annotations/${base} \
           --prefix ${base} \
           --kingdom Bacteria \
           --genus ${GENUS} \
           --species ${SPECIES} \
           --usegenus \
           --compliant \
           --cpus ${THREADS} \
           "$fna"
    # Copiar el archivo .gff a la carpeta de entrada de PIRATE
    cp prokka_annotations/${base}/${base}.gff pirate_input/
done

# ----------------------------------------------------------------------
# PASO 2: Preparar archivos para PIRATE (asegurar extensión .gff)
# ----------------------------------------------------------------------
cd pirate_input
# Si hay archivos con extensión .gff3 o .gff, normalizamos a .gff
for f in *.gff; do
    mv "$f" "${f%.gff}.gff"
done
# Asegurar que todos los .gff tengan un nombre limpio (sin espacios ni caracteres extraños)
rename 's/ /_/g' *.gff 2>/dev/null || true

# ----------------------------------------------------------------------
# PASO 3: Ejecutar PIRATE
# ----------------------------------------------------------------------
echo "=== Ejecutando PIRATE (puede tardar varios minutos) ==="
PIRATE -i ./ -t 2 -s "90,95,98" -a --para-off -o PIRATE_results

# ----------------------------------------------------------------------
# PASO 4: Generar figuras (requiere pirate_plots.py en el PATH)
# ----------------------------------------------------------------------
# Asumiendo que pirate_plots.py está en ~/bacterial-genomics-tutorial/gffs/PIRATE/
if [ -f ~/bacterial-genomics-tutorial/gffs/PIRATE/pirate_plots.py ]; then
    cp ~/bacterial-genomics-tutorial/gffs/PIRATE/pirate_plots.py .
    python pirate_plots.py PIRATE_results/binary_presence_absence.nwk \
                         PIRATE_results/binary_presence_absence.fasta \
                         --labels --format svg
    mv pangenome_*.svg pangenome_report.txt final_results/ 2>/dev/null || true
else
    echo "Advertencia: pirate_plots.py no encontrado. Las figuras no se generaron automáticamente."
    echo "Puedes generar las figuras manualmente después con el script pirate_plots.py."
fi

echo "=== Proceso completado ==="
echo "Resultados en:"
echo "  - Anotaciones: prokka_annotations/"
echo "  - PIRATE: pirate_input/PIRATE_results/"
echo "  - Figuras y reporte: final_results/"
