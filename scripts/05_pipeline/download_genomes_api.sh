#!/bin/bash
# download_genomes_api.sh
# Descarga genomas de referencia usando la API de NCBI

# Crear directorio de trabajo
mkdir -p raw_genomes
cd raw_genomes

# --- Función para descargar ensamblaje dado el nombre de la cepa ---
download_assembly() {
    local strain="$1"
    echo ">> Buscando ensamblaje para: $strain"
    
    # Buscar el assembly en NCBI y obtener el FTP path
    local ftp_path=$(esearch -db assembly -query "$strain" | \
                      esummary | \
                      xtract -pattern DocumentSummary -element FtpPath_GenBank)
    
    if [ -z "$ftp_path" ]; then
        echo "Error: No se encontró ensamblaje para $strain"
        return 1
    fi
    
    local fname=$(basename "$ftp_path")_genomic.fna.gz
    local url="${ftp_path}/${fname}"
    
    echo "Descargando: $url"
    wget -q --show-progress -O "${strain}.fna.gz" "$url"
    gunzip "${strain}.fna.gz"
}

# --- Lista de cepas a descargar ---
download_assembly "L2171"
download_assembly "L2181"
download_assembly "FORC_022"
download_assembly "FORC_023"
download_assembly "RIMD 2210633"
download_assembly "Ba94C2"
download_assembly "M1-1"
download_assembly "CGVP3"
download_assembly "CGVP8"
download_assembly "CGVP22"

echo "Descarga completada. Archivos .fna disponibles:"
ls -lh
