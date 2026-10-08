# Vibrio parahaemolyticus AHPND pangenome

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23225902.svg)](https://doi.org/10.5281/zenodo.23225902)

Reproducible analysis code and key data for:

**"Comparative genomics and experimental virulence of Mexican *Vibrio parahaemolyticus* isolates identify candidate gene-content differences associated with AHPND"**

Hernández-Montiel ÁH, Mata-Torres F, Millán-Aguiñaga N, Torres-Beltrán M, Giffard-Mena I.

## Overview

Comparative genomic analysis of 19 *Vibrio parahaemolyticus* strains:

- 3 Mexican isolates sequenced in this study (11_VM, 15_CESAIBC, 6_VM)
- 16 public reference genomes from Asia and America

Focus: gene-content differences associated with AHPND (Acute Hepatopancreatic Necrosis Disease), driven by the acquisition of a pVA1-like plasmid encoding PirAB toxins.

## Repository structure

```
.
├── README.md
├── LICENSE
├── scripts/
│   ├── 01_pangenome_fisher/     # Pangenome + Fisher exact test + COG
│   ├── 02_huella_digital/       # Gene fingerprint (Roary → PIRATE)
│   ├── 03_filogenia/            # Core-SNP phylogeny, ancestral state
│   ├── 04_bioensayo/            # Survival analysis (Kaplan-Meier)
│   └── 05_pipeline/             # Genome download/annotation scripts
├── input_data/
│   ├── bioensayo/               # Raw survival data (360 individuals)
│   ├── eggnog/                  # eggNOG annotations of differential genes
│   ├── pangenome/               # Exclusive/absent gene lists
│   ├── phylogeny/               # Core-SNP treefile
│   ├── pirAB_references/        # pirA / pirB reference sequences
│   ├── pirate/                  # PIRATE presence/absence matrix
│   └── strain_metadata/         # phenotype.csv + metadata of 19 strains
└── outputs/
    ├── bioensayo/               # Kaplan-Meier PDF + summary tables
    ├── cog/                     # COG enrichment + barplot
    ├── colocalization/          # Co-localization of pirAB + neighbors
    ├── figures/                 # Supplementary Figures S1-S4
    ├── filogenia/               # Core-SNP tree, ancestral reconstruction
    ├── fisher/                  # Fisher results, FDR, sensitivity analyses
    ├── huella/                  # Top 50 fingerprint + heatmaps
    ├── pangenome/               # PIRATE summary, alpha Heaps, accession tree
    └── tables/                  # Supplementary Tables S1-S5
```

## Requirements

| Software | Version |
|----------|---------|
| Prokka | 1.14.6 |
| PIRATE | 1.0.4 |
| Gubbins | 2.4.1 |
| IQ-TREE | 3.1.3 |
| MAFFT | 7.310 |
| R | 4.3.3 |
| Python | 3.7.12 |

### R packages

- ape 5.8.1
- phangorn
- treeio 1.37.0.1
- survival
- rstatix

### Python

- Biopython 1.81
- SciPy
- pandas
- matplotlib
- seaborn

## Execution order

1. `scripts/05_pipeline/download_genomes.sh` — Download reference genomes from NCBI
2. `scripts/05_pipeline/process_genomes.sh` — Prokka annotation
3. `scripts/05_pipeline/run_pirate_all.sh` — Build PIRATE pangenome
4. `scripts/01_pangenome_fisher/fisher_from_pirate.py` — Fisher's exact test
5. `scripts/01_pangenome_fisher/process_fisher.sh` — FDR correction + candidate annotation
6. `scripts/01_pangenome_fisher/02_cog_analysis_roary.R` — COG enrichment
7. `scripts/02_huella_digital/04_huella_pirate_filtrada.R` — Gene fingerprint
8. `scripts/03_filogenia/05_core_snp_phylogeny.sh` — Recombination-filtered core-SNP phylogeny
9. `scripts/03_filogenia/07_ancestral_reconstruction.R` — Ancestral state reconstruction
10. `scripts/03_filogenia/08_panel_ancestral.R` — Panel of ancestral trees
11. `scripts/03_filogenia/10_curvas_acumulacion.R` — Randomized accumulation curves
12. `scripts/04_bioensayo/05_kaplan_meier.R` — Survival analysis

## Data availability

### Newly sequenced genomes

| Strain | BioProject | BioSample | GenBank assembly |
|---|---|---|---|
| 11_VM | PRJNA1499920 | SAMN61896695 | JCBNQA000000000 |
| 15_CESAIBC | PRJNA1499939 | SAMN61897767 | JCBNTU000000000 |
| 6_VM | PRJNA1499907 | SAMN61896394 | JCBNPX000000000 |

Reference genome accessions are listed in `input_data/strain_metadata/phenotype.csv`.

### Archived version

This repository is archived on Zenodo: [10.5281/zenodo.23225902](https://doi.org/10.5281/zenodo.23225902)

## Citation

If you use this code, please cite the manuscript above and the Zenodo release:

```
Hernández-Montiel ÁH, Mata-Torres F, Millán-Aguiñaga N, Torres-Beltrán M, Giffard-Mena I.
Comparative genomics and experimental virulence of Mexican Vibrio parahaemolyticus isolates
identify candidate gene-content differences associated with AHPND.
[Journal TBD], [Year]. DOI: [pending]

Code and data: https://doi.org/10.5281/zenodo.23225902
```

## License

MIT — see `LICENSE`.
