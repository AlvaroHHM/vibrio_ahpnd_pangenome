#!/bin/bash
# Ejecutar PIRATE con todos los archivos GFF de la carpeta gffs/

# Directorio de trabajo
cd ~/bacterial-genomics-tutorial/gffs

# Crear enlaces simbólicos .gff para todos los .gff3 (si no existen ya)
echo "Creando enlaces .gff para archivos .gff3..."
for f in *.gff3; do
    base=$(basename "$f" .gff3)
    if [ ! -f "${base}.gff" ]; then
        ln -s "$f" "${base}.gff"
    fi
done

# Verificar que hay archivos .gff
num_gff=$(ls -1 *.gff 2>/dev/null | wc -l)
if [ "$num_gff" -eq 0 ]; then
    echo "Error: No se encontraron archivos .gff en el directorio."
    exit 1
fi

echo "Se encontraron $num_gff archivos .gff. Ejecutando PIRATE..."

# Ejecutar PIRATE con parámetros optimizados para memoria
PIRATE -i ./ \
       -t 2 \
       -s "90,95,98" \
       -a \
       --para-off

echo "PIRATE finalizado. Los resultados están en ./PIRATE/"
