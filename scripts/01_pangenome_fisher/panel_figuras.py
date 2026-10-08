#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genera un panel combinando tres figuras del pangenoma:
  - superior izquierda: pangenome_stacked_bars_percent
  - superior derecha:   pangenome_accumulation
  - inferior (ancho completo): pangenome_matrix

Se asume que existen las versiones PNG (o se convierten desde SVG).
El panel se guarda como panel_figuras.png y panel_figuras.svg
"""

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

# ---------- Configuración ----------
fig1 = "pangenome_stacked_bars_percent.png"   # arriba izquierda
fig2 = "pangenome_accumulation.png"           # arriba derecha
fig3 = "pangenome_matrix.png"                 # abajo ancho completo
output_png = "panel_figuras.png"
output_svg = "panel_figuras.svg"

# Si no existen los PNG, intentar convertir desde SVG con ImageMagick
def ensure_png(png_file):
    if os.path.exists(png_file):
        return png_file
    svg_file = png_file.replace('.png', '.svg')
    if os.path.exists(svg_file):
        print(f"Convirtiendo {svg_file} a PNG...")
        os.system(f"convert {svg_file} {png_file}")   # requiere ImageMagick
        if os.path.exists(png_file):
            return png_file
    return None

f1 = ensure_png(fig1)
f2 = ensure_png(fig2)
f3 = ensure_png(fig3)

if not all([f1, f2, f3]):
    print("Error: faltan archivos de figuras. Verifica que existan los PNG o SVG.")
    exit(1)

# ---------- Crear panel ----------
fig, axes = plt.subplots(2, 2, figsize=(18, 14),
                         gridspec_kw={'height_ratios': [1, 1.6], 'width_ratios': [1, 1]})

# Ocultar el eje vacío inferior derecho
axes[1, 1].axis('off')

# Cargar imágenes
img1 = mpimg.imread(f1)
img2 = mpimg.imread(f2)
img3 = mpimg.imread(f3)

# Superior izquierda
ax1 = axes[0, 0]
ax1.imshow(img1)
ax1.axis('off')
ax1.set_title('Gene categories per genome (relative)', fontsize=12)

# Superior derecha
ax2 = axes[0, 1]
ax2.imshow(img2)
ax2.axis('off')
ax2.set_title('Pangenome accumulation curves', fontsize=12)

# Inferior (ocupará todo el ancho)
ax3 = axes[1, 0]
ax3.imshow(img3)
ax3.axis('off')
ax3.set_title('Pangenome matrix with accessory-genome tree', fontsize=12)

# Ajustar para que la imagen inferior ocupe también el espacio derecho
# Fusionamos el eje inferior derecho con el izquierdo
gs = axes[1, 0].get_gridspec()
for ax in axes[1, :]:
    ax.remove()
ax_bottom = fig.add_subplot(gs[1, :])
ax_bottom.imshow(img3)
ax_bottom.axis('off')
ax_bottom.set_title('Pangenome matrix with accessory-genome tree', fontsize=12)

plt.tight_layout()
plt.savefig(output_png, dpi=300, bbox_inches='tight')
plt.savefig(output_svg, dpi=300, bbox_inches='tight')
print(f"Panel guardado como {output_png} y {output_svg}")
