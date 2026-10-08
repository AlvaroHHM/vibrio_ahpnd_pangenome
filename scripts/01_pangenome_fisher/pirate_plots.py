#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Figure 1 for the AHPND pangenome manuscript.

Ensamblado final con PIL para evitar el downsampling de matplotlib
en composites multi-panel.

Output: outputs/figures/Figure_1_pangenome.png / .pdf
"""

import argparse
import os
import tempfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio import Phylo
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.image import imread
from PIL import Image

plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.grid'] = False
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'


def read_binary_fasta(fasta_file):
    with open(fasta_file) as f:
        lines = [l.strip() for l in f if l.strip()]
    genomes, seqs = [], []
    cur_g, cur_s = None, []
    for line in lines:
        if line.startswith('>'):
            if cur_g is not None:
                seqs.append(''.join(cur_s))
            cur_g = line[1:]
            genomes.append(cur_g)
            cur_s = []
        else:
            cur_s.append(line)
    if cur_g is not None:
        seqs.append(''.join(cur_s))
    lengths = [len(s) for s in seqs]
    if len(set(lengths)) != 1:
        raise ValueError("Las secuencias binarias no tienen la misma longitud.")
    chars = set(seqs[0])
    if chars.issubset({'0', '1'}):
        mapping = {'0': 0, '1': 1}
    elif chars.issubset({'A', 'C'}):
        mapping = {'A': 1, 'C': 0}
    elif chars.issubset({'a', 'c'}):
        mapping = {'a': 1, 'c': 0}
    else:
        mapping = {c: 0 if c == '0' else 1 for c in chars}
    data = [[mapping.get(ch, 1) for ch in seq] for seq in seqs]
    return pd.DataFrame(data, index=genomes).T


def accumulation_curves(presence_matrix, order):
    pangenome_sizes, core_sizes = [], []
    current_genes = set()
    for i, strain in enumerate(order):
        strain_genes = set(presence_matrix.index[presence_matrix[strain] == 1])
        current_genes.update(strain_genes)
        pangenome_sizes.append(len(current_genes))
        core_genes = strain_genes if i == 0 else core_genes.intersection(strain_genes)
        core_sizes.append(len(core_genes))
    return pangenome_sizes, core_sizes


def render_panel_A(counts, outpath, figsize=(8, 6)):
    """Panel A: barras apiladas de categorías génicas."""
    fig, ax = plt.subplots(figsize=figsize, dpi=200)
    counts.plot(kind='bar', stacked=True, ax=ax, colormap='Blues', edgecolor='black')
    ax.set_xlabel('Genome', fontsize=12)
    ax.set_ylabel('Number of genes', fontsize=12)
    ax.set_title('Gene categories per genome', pad=12, fontsize=13)
    ax.legend(title='Category', loc='upper right', fontsize=9)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.text(-0.18, 1.08, 'A)', transform=ax.transAxes,
            fontsize=22, fontweight='bold', va='top', ha='left')
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()


def render_panel_B(x, pangenome_sizes, core_sizes, outpath, figsize=(8, 6)):
    """Panel B: curvas de acumulación."""
    fig, ax = plt.subplots(figsize=figsize, dpi=200)
    ax.plot(x, pangenome_sizes, 'o-', label='Pangenome size', color='#3B4992', linewidth=2)
    ax.plot(x, core_sizes, 's-', label='Core genome size', color='#EE0000', linewidth=2)
    ax.set_xlabel('Number of genomes', fontsize=12)
    ax.set_ylabel('Number of gene families', fontsize=12)
    ax.set_title('Accumulation curves', pad=12, fontsize=13)
    ax.legend(fontsize=11)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.text(-0.18, 1.08, 'B)', transform=ax.transAxes,
            fontsize=22, fontweight='bold', va='top', ha='left')
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()


def render_panel_C_tree(tree, mdist, n_strains, outpath, figsize=(3, 8), show_labels=True):
    """Panel C: árbol filogenético de genes accesorios."""
    fig, ax = plt.subplots(figsize=figsize, dpi=200)
    fsize = max(12, 16 - 0.15 * n_strains)
    with plt.rc_context({'font.size': fsize}):
        Phylo.draw(
            tree, axes=ax, show_confidence=False,
            label_func=(lambda x: None) if not show_labels else (lambda x: str(x)),
            xlim=(-mdist * 0.1, mdist * 2.5),
            do_show=False,
        )
    ax.set_axis_off()
    ax.text(-0.25, 1.05, 'C)', transform=ax.transAxes,
            fontsize=22, fontweight='bold', va='top', ha='left')
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()


def render_panel_C_matrix(matrix_binned, core_frac, outpath,
                          figsize=(14, 8), n_bins=500):
    """Panel C: matriz de presencia/ausencia."""
    contrast_cmap = LinearSegmentedColormap.from_list(
        'contrast_blues',
        ['#FFFFFF', '#D6E4F0', '#6BAED6', '#2171B5', '#08306B'],
        N=256,
    )
    fig, ax = plt.subplots(figsize=figsize, dpi=200)
    ax.imshow(matrix_binned, cmap=contrast_cmap, vmin=0, vmax=1,
              aspect='auto', interpolation='nearest')
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_title('Accessory-gene presence/absence matrix', pad=12, fontsize=13)
    ax.axvline(x=core_frac * n_bins - 0.5, color='red',
               linestyle='--', linewidth=2, alpha=0.9)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()


def compose_figure(panel_A, panel_B, panel_C_tree, panel_C_matrix, outpath_png,
                   outpath_pdf):
    """Ensambla los 4 paneles en una figura 2×2 con PIL."""
    img_A = Image.open(panel_A).convert('RGB')
    img_B = Image.open(panel_B).convert('RGB')
    img_C_tree = Image.open(panel_C_tree).convert('RGB')
    img_C_matrix = Image.open(panel_C_matrix).convert('RGB')

    # --- Alinear alturas de fila 1 (A y B) ---
    h_top = max(img_A.height, img_B.height)
    if img_A.height < h_top:
        new_w = int(img_A.width * h_top / img_A.height)
        img_A = img_A.resize((new_w, h_top), Image.LANCZOS)
    if img_B.height < h_top:
        new_w = int(img_B.width * h_top / img_B.height)
        img_B = img_B.resize((new_w, h_top), Image.LANCZOS)

    # --- Alinear alturas de fila 2 (árbol y matriz) ---
    h_bot = max(img_C_tree.height, img_C_matrix.height)
    if img_C_tree.height < h_bot:
        new_w = int(img_C_tree.width * h_bot / img_C_tree.height)
        img_C_tree = img_C_tree.resize((new_w, h_bot), Image.LANCZOS)
    if img_C_matrix.height < h_bot:
        new_w = int(img_C_matrix.width * h_bot / img_C_matrix.height)
        img_C_matrix = img_C_matrix.resize((new_w, h_bot), Image.LANCZOS)

    # --- Anchos de cada fila ---
    w_top = img_A.width + img_B.width
    w_bot = img_C_tree.width + img_C_matrix.width
    w_total = max(w_top, w_bot)

    # --- Padding ---
    pad = 40
    h_total = h_top + h_bot + 3 * pad

    composite = Image.new('RGB', (w_total + 2 * pad, h_total), 'white')

    # Fila 1: A | B (centrado)
    x_A = pad + (w_total - w_top) // 2
    composite.paste(img_A, (x_A, pad))
    composite.paste(img_B, (x_A + img_A.width, pad))

    # Fila 2: C_tree | C_matrix (centrado)
    x_bot = pad + (w_total - w_bot) // 2
    y_bot = pad + h_top + pad
    composite.paste(img_C_tree, (x_bot, y_bot))
    composite.paste(img_C_matrix, (x_bot + img_C_tree.width, y_bot))

    # --- Guardar ---
    composite.save(outpath_png, 'PNG')
    composite.save(outpath_pdf, 'PDF', resolution=200.0)
    return composite.size


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('tree')
    parser.add_argument('fasta')
    parser.add_argument('--outdir', default='outputs/figures')
    parser.add_argument('--no-labels', action='store_true')
    parser.add_argument('--n-bins', type=int, default=500)
    args = parser.parse_args()

    # ---------- Cargar datos ----------
    t = Phylo.read(args.tree, 'newick')
    leaf_order = [x.name for x in t.get_terminals()]
    roary = read_binary_fasta(args.fasta)

    missing = [n for n in leaf_order if n not in roary.columns]
    if missing:
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

    # ---------- Datos Panel A ----------
    core_genes = roary.index[(gene_freq >= n_strains * thresholds['core']) & (gene_freq <= n_strains)]
    softcore_genes = roary.index[(gene_freq >= n_strains * thresholds['softcore']) & (gene_freq < n_strains * thresholds['core'])]
    shell_genes = roary.index[(gene_freq >= n_strains * thresholds['shell']) & (gene_freq < n_strains * thresholds['softcore'])]
    cloud_genes = roary.index[gene_freq < n_strains * thresholds['shell']]

    def count_category(genes_list):
        return roary.loc[genes_list].sum(axis=0) if len(genes_list) > 0 else pd.Series(0, index=roary.columns)

    counts = pd.DataFrame({
        'core':     count_category(core_genes),
        'softcore': count_category(softcore_genes),
        'shell':    count_category(shell_genes),
        'cloud':    count_category(cloud_genes),
    }).loc[leaf_order]

    # ---------- Datos Panel B ----------
    pangenome_sizes, core_sizes = accumulation_curves(roary, leaf_order)
    x = list(range(1, n_strains + 1))

    # ---------- Datos Panel C ----------
    idx = gene_freq.sort_values(ascending=False).index
    roary_sorted = roary.loc[idx][leaf_order]
    mdist = max([t.distance(t.root, x) for x in t.get_terminals()])

    matrix_full = np.ascontiguousarray(roary_sorted.T.values.astype(float))
    n_rows, n_cols = matrix_full.shape
    n_bins = args.n_bins
    bin_size = n_cols // n_bins
    matrix_binned = matrix_full[:, :bin_size * n_bins].reshape(
        n_rows, n_bins, bin_size
    ).mean(axis=2)

    print(f"Matriz binned: shape={matrix_binned.shape}, "
          f"min={matrix_binned.min():.2f}, max={matrix_binned.max():.2f}, "
          f"mean={matrix_binned.mean():.2f}")

    # ---------- Renderizar cada panel a PNG temporal ----------
    tmp = tempfile.gettempdir()
    pA = os.path.join(tmp, 'fig1_A.png')
    pB = os.path.join(tmp, 'fig1_B.png')
    pC_tree = os.path.join(tmp, 'fig1_Ctree.png')
    pC_matrix = os.path.join(tmp, 'fig1_Cmatrix.png')

    render_panel_A(counts, pA, figsize=(8, 6))
    render_panel_B(x, pangenome_sizes, core_sizes, pB, figsize=(8, 6))
    render_panel_C_tree(t, mdist, n_strains, pC_tree, figsize=(3, 8),
                        show_labels=not args.no_labels)
    render_panel_C_matrix(matrix_binned, core_count / total_genes,
                          pC_matrix, figsize=(14, 8), n_bins=n_bins)

    # ---------- Componer con PIL ----------
    os.makedirs(args.outdir, exist_ok=True)
    outpath_png = os.path.join(args.outdir, 'Figure_1_pangenome.png')
    outpath_pdf = os.path.join(args.outdir, 'Figure_1_pangenome.pdf')

    size = compose_figure(pA, pB, pC_tree, pC_matrix,
                          outpath_png, outpath_pdf)

    # ---------- Reporte ----------
    with open('outputs/pangenome/pangenome_report.txt', 'w') as f:
        f.write("PANGENOME ANALYSIS REPORT\n")
        f.write("=========================\n\n")
        f.write(f"Number of genomes analyzed: {n_strains}\n")
        f.write(f"Total gene families: {total_genes}\n")
        f.write(f"Core genome size (>=99%): {core_count}\n")
        f.write(f"Soft-core (95-99%): {softcore_count}\n")
        f.write(f"Shell (15-95%): {shell_count}\n")
        f.write(f"Cloud (<15%): {cloud_count}\n")

    print("")
    print("=========================================")
    print("Figure 1 generated (composed with PIL)")
    print("=========================================")
    print(f"PNG : {outpath_png}")
    print(f"PDF : {outpath_pdf}")
    print(f"Size: {size[0]} x {size[1]} px")
    print(f"Genomes    : {n_strains}")
    print(f"Gene fams  : {total_genes}")
    print(f"Core       : {core_count} ({100*core_count/total_genes:.1f}%)")
    print(f"Accessory  : {shell_count + cloud_count} ({(shell_count+cloud_count)*100/total_genes:.1f}%)")


if __name__ == "__main__":
    main()
