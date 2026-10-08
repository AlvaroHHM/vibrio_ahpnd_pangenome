#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script para visualización de resultados de PIRATE (versión corregida).
Genera:
  - pangenome_frequency (histograma)
  - pangenome_matrix (matriz de presencia/ausencia con árbol + nombres + colores por clado)
  - pangenome_pie (gráfico circular de categorías)
  - pangenome_accumulation (curvas de acumulación del pangenoma y core)
  - pangenome_stacked_bars (barras apiladas por genoma) absolutas y relativas
  - pangenome_panel (collage de las figuras principales, sin la matriz)
  - pangenome_report.txt (interpretación en texto)

Todas las figuras se guardan en formato PNG y SVG.
Opcional: --clades <archivo.tsv> para colorear ramas según grupos.
"""

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from Bio import Phylo
import sys

sns.set_style('white')

# Paleta AAAS de ggsci
AAAS_COLORS = [
    '#3B4992',  # azul
    '#EE0000',  # rojo
    '#008B45',  # verde
    '#631879',  # púrpura
    '#008280',  # verde azulado
    '#BB0021',  # rojo oscuro
    '#5F559B',  # azul violáceo
    '#A20056',  # magenta
    '#808180',  # gris
    '#1B1919',  # negro
]

def read_clade_mapping(filepath):
    """Lee archivo TSV con columnas 'strain' y 'clade'."""
    mapping = {}
    with open(filepath) as f:
        header = f.readline().strip().split('\t')
        # Asumir que la primera columna es strain y la segunda es clade
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping

def assign_colors_to_tree(tree, strain_to_clade, clade_to_color):
    """Asigna colores a los clados del árbol según el grupo de sus hojas."""
    # Primero, asignar color y grupo a las hojas
    for leaf in tree.get_terminals():
        clade = strain_to_clade.get(leaf.name, None)
        if clade:
            leaf.group = clade
            leaf.color = clade_to_color[clade]
        else:
            leaf.group = None
            leaf.color = 'black'

    # Propagar colores hacia arriba si todas las hojas de un subárbol comparten grupo
    def propagate(clade):
        if clade.is_terminal():
            return clade.group
        child_groups = set()
        for child in clade.clades:
            group = propagate(child)
            if group:
                child_groups.add(group)
        if len(child_groups) == 1:
            clade.group = next(iter(child_groups))
            clade.color = clade_to_color[clade.group]
        else:
            clade.group = None
            clade.color = 'black'
        return clade.group

    propagate(tree.root)
    # También colorear la rama desde la raíz (opcional)
    if tree.root.group:
        tree.root.color = clade_to_color[tree.root.group]

def save_fig(base_name, dpi=300):
    """Guarda la figura actual en PNG y SVG."""
    plt.savefig(f'{base_name}.png', dpi=dpi, format='png', bbox_inches='tight')
    plt.savefig(f'{base_name}.svg', dpi=dpi, format='svg', bbox_inches='tight')
    plt.close()

def read_binary_fasta(fasta_file):
    with open(fasta_file) as f:
        lines = [l.strip() for l in f if l.strip()]
    genomes = []
    seqs = []
    current_genome = None
    current_seq = []
    for line in lines:
        if line.startswith('>'):
            if current_genome is not None:
                seqs.append(''.join(current_seq))
            current_genome = line[1:]
            genomes.append(current_genome)
            current_seq = []
        else:
            current_seq.append(line)
    if current_genome is not None:
        seqs.append(''.join(current_seq))
    lengths = [len(s) for s in seqs]
    if len(set(lengths)) != 1:
        raise ValueError("Las secuencias binarias no tienen la misma longitud.")
    num_genes = lengths[0]
    chars = set(seqs[0])
    if chars.issubset({'0','1'}):
        mapping = {'0':0, '1':1}
    elif chars.issubset({'A','C'}):
        mapping = {'A':1, 'C':0}
    elif chars.issubset({'a','c'}):
        mapping = {'a':1, 'c':0}
    else:
        print(f"Advertencia: caracteres {chars}. Se asume presencia si no es '0'.")
        mapping = {c:0 if c=='0' else 1 for c in chars}
    data = []
    for seq in seqs:
        row = [mapping.get(ch, 1) for ch in seq]
        data.append(row)
    df = pd.DataFrame(data, index=genomes).T
    return df

def accumulation_curves(presence_matrix, order):
    n_genomes = len(order)
    pangenome_sizes = []
    core_sizes = []
    current_genes = set()
    for i, strain in enumerate(order):
        strain_genes = set(presence_matrix.index[presence_matrix[strain] == 1])
        current_genes.update(strain_genes)
        pangenome_sizes.append(len(current_genes))
        if i == 0:
            core_genes = strain_genes
        else:
            core_genes = core_genes.intersection(strain_genes)
        core_sizes.append(len(core_genes))
    return pangenome_sizes, core_sizes

def create_panel_from_data(hist_data, pie_counts, stacked_df, accum_x, accum_pan, accum_core, n_strains, output_base):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    ax = axes[0,0]
    ax.hist(hist_data, bins=n_strains, histtype='stepfilled', alpha=0.7, color='steelblue')
    ax.set_xlabel('No. of genomes')
    ax.set_ylabel('No. of genes')
    sns.despine(ax=ax, left=True, bottom=True)
    ax.set_title('Gene frequency distribution')
    ax = axes[0,1]
    total = sum(pie_counts)
    def my_autopct(pct):
        val = int(round(pct * total / 100.0))
        return f'{val:d}'
    labels = [f'core\n(≥99%)', f'soft-core\n(95-99%)', f'shell\n(15-95%)', f'cloud\n(<15%)']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    ax.pie(pie_counts, labels=labels, autopct=my_autopct, startangle=90, colors=colors)
    ax.set_title('Pangenome categories')
    ax = axes[1,0]
    if stacked_df.shape[0] > 15:
        plot_df = stacked_df.iloc[:15]
    else:
        plot_df = stacked_df
    plot_df.plot(kind='bar', stacked=True, ax=ax, colormap='Blues', edgecolor='black')
    ax.set_xlabel('Genome')
    ax.set_ylabel('Number of genes')
    ax.set_title('Gene categories per genome')
    ax.legend(title='Category')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    ax = axes[1,1]
    ax.plot(accum_x, accum_pan, 'o-', label='Pangenome size')
    ax.plot(accum_x, accum_core, 's-', label='Core genome size')
    ax.set_xlabel('Number of genomes')
    ax.set_ylabel('Number of gene families')
    ax.legend()
    sns.despine(ax=ax)
    ax.set_title('Accumulation curves')
    plt.tight_layout()
    save_fig(output_base, dpi=300)

def main():
    parser = argparse.ArgumentParser(description='Create plots from PIRATE outputs')
    parser.add_argument('tree', help='Newick tree file (binary_presence_absence.nwk)')
    parser.add_argument('fasta', help='Binary FASTA file (binary_presence_absence.fasta)')
    parser.add_argument('--clades', help='Archivo TSV con columnas strain y clade para colorear ramas')
    parser.add_argument('--no-labels', action='store_true',
                        help='No mostrar los nombres de las ramas en el árbol (por defecto se muestran)')
    args = parser.parse_args()

    # Leer árbol
    t = Phylo.read(args.tree, 'newick')

    # Si se proporciona archivo de clados, asignar colores
    if args.clades:
        strain_to_clade = read_clade_mapping(args.clades)
        # Obtener lista única de clados en orden de aparición
        clade_list = list(dict.fromkeys(strain_to_clade.values()))
        clade_to_color = {clade: AAAS_COLORS[i % len(AAAS_COLORS)] for i, clade in enumerate(clade_list)}
        assign_colors_to_tree(t, strain_to_clade, clade_to_color)
        print(f"Se colorearon las ramas según {len(clade_list)} grupos")
    else:
        print("No se proporcionó archivo de clados; el árbol se dibujará sin colores")

    leaf_order = [x.name for x in t.get_terminals()]
    roary = read_binary_fasta(args.fasta)

    missing = [n for n in leaf_order if n not in roary.columns]
    if missing:
        print(f"Advertencia: algunos genomas del árbol no están en la matriz: {missing}")
        leaf_order = [n for n in leaf_order if n in roary.columns]
    roary = roary[leaf_order]

    n_strains = roary.shape[1]
    total_genes = roary.shape[0]
    gene_freq = roary.sum(axis=1)

    thresholds = {'core': 0.99, 'softcore': 0.95, 'shell': 0.15}
    core_count = ((gene_freq >= n_strains * thresholds['core']) & (gene_freq <= n_strains)).sum()
    softcore_count = ((gene_freq >= n_strains * thresholds['softcore']) & (gene_freq < n_strains * thresholds['core'])).sum()
    shell_count = ((gene_freq >= n_strains * thresholds['shell']) & (gene_freq < n_strains * thresholds['softcore'])).sum()
    cloud_count = (gene_freq < n_strains * thresholds['shell']).sum()

    # 1) Histograma de frecuencias
    plt.figure(figsize=(7,5))
    plt.hist(gene_freq, bins=n_strains, histtype='stepfilled', alpha=0.7, color='steelblue')
    plt.xlabel('No. of genomes')
    plt.ylabel('No. of genes')
    sns.despine(left=True, bottom=True)
    save_fig('pangenome_frequency')

    # 2) Matriz con árbol (SIEMPRE con nombres, y opcionalmente colores)
    idx = gene_freq.sort_values(ascending=False).index
    roary_sorted = roary.loc[idx]
    roary_sorted = roary_sorted[leaf_order]
    mdist = max([t.distance(t.root, x) for x in t.get_terminals()])

    with sns.axes_style('whitegrid'):
        fig = plt.figure(figsize=(17,10))
        ax1 = plt.subplot2grid((1,40), (0,10), colspan=30)
        ax1.matshow(roary_sorted.T, cmap=plt.cm.Blues, vmin=0, vmax=1, aspect='auto', interpolation='none')
        ax1.set_yticks([])
        ax1.set_xticks([])
        ax1.axis('off')
        ax = plt.subplot2grid((1,40), (0,0), colspan=10, facecolor='white')
        fig.subplots_adjust(wspace=0, hspace=0)
        ax1.set_title(f'PIRATE matrix\n({total_genes} gene clusters)')

        fsize = max(5, 12 - 0.25 * n_strains)
        with plt.rc_context({'font.size': fsize}):
            Phylo.draw(t, axes=ax, show_confidence=False,
                       label_func=(lambda x: None) if args.no_labels else (lambda x: str(x)),
                       xlim=(-mdist*0.1, mdist*1.8),
                       do_show=False)
        ax.set_axis_off()
        ax.set_title(f'Tree\n({n_strains} strains)')
        save_fig('pangenome_matrix')

    # 3) Gráfico circular
    def my_autopct(pct):
        val = int(round(pct * total_genes / 100.0))
        return f'{val:d}'
    plt.figure(figsize=(10,10))
    plt.pie([core_count, softcore_count, shell_count, cloud_count],
            labels=[f'core\n(≥ {int(n_strains*thresholds["core"])} strains)',
                    f'soft-core\n({int(n_strains*thresholds["softcore"])}–{int(n_strains*thresholds["core"])-1} strains)',
                    f'shell\n({int(n_strains*thresholds["shell"])}–{int(n_strains*thresholds["softcore"])-1} strains)',
                    f'cloud\n(< {int(n_strains*thresholds["shell"])} strains)'],
            explode=[0.1,0.05,0.02,0], radius=0.9,
            colors=[(0,0,1, x/total_genes) for x in (core_count, softcore_count, shell_count, cloud_count)],
            autopct=my_autopct)
    save_fig('pangenome_pie')

    # 4) Curvas de acumulación
    pangenome_sizes, core_sizes = accumulation_curves(roary, leaf_order)
    x = range(1, n_strains+1)
    plt.figure(figsize=(7,5))
    plt.plot(x, pangenome_sizes, 'o-', label='Pangenome size')
    plt.plot(x, core_sizes, 's-', label='Core genome size')
    plt.xlabel('Number of genomes')
    plt.ylabel('Number of gene families')
    plt.legend()
    sns.despine()
    save_fig('pangenome_accumulation')

    # 5) Barras apiladas absolutas
    core_genes = roary.index[(gene_freq >= n_strains * thresholds['core']) & (gene_freq <= n_strains)]
    softcore_genes = roary.index[(gene_freq >= n_strains * thresholds['softcore']) & (gene_freq < n_strains * thresholds['core'])]
    shell_genes = roary.index[(gene_freq >= n_strains * thresholds['shell']) & (gene_freq < n_strains * thresholds['softcore'])]
    cloud_genes = roary.index[gene_freq < n_strains * thresholds['shell']]

    def count_category(genes_list):
        return roary.loc[genes_list].sum(axis=0) if len(genes_list) > 0 else pd.Series(0, index=roary.columns)

    counts = pd.DataFrame({
        'core': count_category(core_genes),
        'softcore': count_category(softcore_genes),
        'shell': count_category(shell_genes),
        'cloud': count_category(cloud_genes)
    })
    counts = counts.loc[leaf_order]

    ax = counts.plot(kind='bar', stacked=True, figsize=(8,6), colormap='Blues', edgecolor='black')
    ax.set_xlabel('Genome')
    ax.set_ylabel('Number of genes')
    ax.set_title('Gene categories per genome')
    ax.legend(title='Category')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    save_fig('pangenome_stacked_bars')

    # 6) Barras apiladas relativas (%)
    totals = counts.sum(axis=1)
    counts_pct = counts.div(totals, axis=0) * 100
    ax = counts_pct.plot(kind='bar', stacked=True, figsize=(8,6), colormap='Blues', edgecolor='black')
    ax.set_xlabel('Genome')
    ax.set_ylabel('Percentage of genes')
    ax.set_title('Gene categories per genome (relative)')
    ax.legend(title='Category')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    save_fig('pangenome_stacked_bars_percent')

    # 7) Panel resumen
    create_panel_from_data(gene_freq, (core_count, softcore_count, shell_count, cloud_count),
                           counts, x, pangenome_sizes, core_sizes, n_strains, 'pangenome_panel')

    # 8) Reporte de texto
    with open('pangenome_report.txt', 'w') as f:
        f.write("PANGENOME ANALYSIS REPORT\n")
        f.write("=========================\n\n")
        f.write(f"Number of genomes analyzed: {n_strains}\n")
        f.write(f"Total gene families (pangenome size): {total_genes}\n")
        f.write(f"Core genome size (≥99% strains): {core_count}\n")
        f.write(f"Soft-core genome size (95-99%): {softcore_count}\n")
        f.write(f"Shell genome size (15-95%): {shell_count}\n")
        f.write(f"Cloud genome size (<15%): {cloud_count}\n\n")
        f.write("Interpretation:\n")
        if n_strains >= 5:
            if len(pangenome_sizes) >= 3:
                last_slope = (pangenome_sizes[-1] - pangenome_sizes[-2]) / pangenome_sizes[-2]
                if last_slope < 0.05:
                    f.write("- The pangenome accumulation curve is approaching a plateau, suggesting a closed pangenome.\n")
                else:
                    f.write("- The pangenome accumulation curve continues to rise, indicating an open pangenome.\n")
        else:
            f.write("- With fewer than 5 genomes, the pangenome status cannot be reliably assessed. More strains are needed.\n")
        f.write(f"- The core genome comprises {core_count} genes ({100*core_count/total_genes:.1f}% of the total).\n")
        f.write(f"- The accessory genome (shell + cloud) comprises {shell_count+cloud_count} genes ({(shell_count+cloud_count)*100/total_genes:.1f}%).\n")
        f.write("\nPer-genome statistics:\n")
        for strain in leaf_order:
            total_in_strain = roary[strain].sum()
            f.write(f"  {strain}: {total_in_strain} genes\n")
        f.write("\nRecommendations:\n")
        f.write("- Consider adding more genomes to better characterise the pangenome openness.\n")
        f.write("- Investigate the functional annotation of shell and cloud genes; many may be associated with mobile genetic elements.\n")
        f.write("- If you have phenotypic data, perform a pan-GWAS to link specific genes to traits.\n")

    print("Report written to pangenome_report.txt")
    print("Panel saved as pangenome_panel.png and pangenome_panel.svg")
    print("All figures were saved in PNG and SVG format.")

if __name__ == "__main__":
    main()
