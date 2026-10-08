import csv
from scipy.stats import fisher_exact

# Leer fenotipos
pheno = {}
with open('phenotype.csv') as pf:
    reader = csv.reader(pf)
    next(reader)  # saltar encabezado
    for row in reader:
        if not row:
            continue
        strain = row[0].strip()
        pheno[strain] = 1 if row[1].strip() == '1' else 0

# Abrir PIRATE.gene_families.ordered.tsv
with open('PIRATE.gene_families.ordered.tsv') as f:
    header = f.readline().strip().split('\t')
    
    # Encontrar índices de las columnas de cepas (aquellas cuyos nombres están en pheno)
    strain_cols = []
    strain_names = []
    for idx, col_name in enumerate(header):
        if col_name in pheno:
            strain_cols.append(idx)
            strain_names.append(col_name)
    
    # Verificar que tenemos todas las cepas del fenotipo
    missing = [s for s in pheno if s not in strain_names]
    if missing:
        print(f"Warning: las siguientes cepas no se encontraron en la cabecera: {missing}")
    print(f"Se detectaron {len(strain_names)} cepas en las columnas {strain_cols[0]} a {strain_cols[-1]}")
    
    # Preparar archivo de salida
    out = open('fisher_results_corrected.tsv', 'w')
    out.write('gene\tpresent_in_AHPND+\tabsent_in_AHPND+\tpresent_in_AHPND-\tabsent_in_AHPND-\todds_ratio\tp_value\n')
    
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) < 2:
            continue
        gene = cols[1]  # gene_family (columna 2)
        
        # Extraer presencia/ausencia de las columnas detectadas
        presence = []
        for idx in strain_cols:
            if idx < len(cols) and cols[idx].strip() != '':
                presence.append(1)
            else:
                presence.append(0)
        
        # Contar a, b, c, d
        a = b = c = d = 0
        for strain, pres in zip(strain_names, presence):
            if pheno[strain] == 1:
                if pres == 1:
                    a += 1
                else:
                    b += 1
            else:
                if pres == 1:
                    c += 1
                else:
                    d += 1
        
        # Test exacto de Fisher
        if a+b > 0 and c+d > 0 and a+c > 0 and b+d > 0:
            odds, p = fisher_exact([[a, b], [c, d]])
        else:
            odds, p = float('inf'), 1.0
        
        out.write(f'{gene}\t{a}\t{b}\t{c}\t{d}\t{odds}\t{p}\n')
    
    out.close()
    print("Done. Results in fisher_results_corrected.tsv")
