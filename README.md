# Vibrio parahaemolyticus AHPND pangenome

Reproducible analysis code and key data for:

**"Comparative genomics and experimental virulence of Mexican *Vibrio parahaemolyticus* isolates identify candidate gene-content differences associated with AHPND"**

Hernández-Montiel ÁH, Mata-Torres F, Millán-Aguiñaga N, Torres-Beltrán M, Giffard-Mena I.

## Overview

Comparative genomic analysis of 19 *Vibrio parahaemolyticus* strains:

- 3 Mexican isolates sequenced in this study (11_VM, 15_CESAIBC, 6_VM)
- 16 public reference genomes from Asia and America

Focus: gene-content differences associated with AHPND (Acute Hepatopancreatic Necrosis Disease), driven by the acquisition of a pVA1-like plasmid encoding PirAB toxins.

## Repository structure

    .
    ├── README.md
    ├── LICENSE
    ├── scripts/
    │   ├── 01_pangenome_fisher/     # Pangenome + Fisher exact test + COG
    │   ├── 02_huella_digital/       # Gene fingerprint (Roary → PIRATE)
    │   ├── 03_filogenia/            # Core-SNP phylogeny, ancestral state
    │   ├── 04_bioensayo/            # Survival analysis
    │   └── 05_pipeline/             # Genome download/annotation scripts
    ├── input_data/
    │   ├── strain_metadata/         # phenotype.csv
    │   ├── gffs_mexican_isolates/   # Prokka GFFs (3 Mexican strains)
    │   └── pirAB_references/        # pirAB reference sequences
    └── outputs/
        ├── pangenome/               # PIRATE output, accession tree
        ├── fisher/                  # Fisher results, FDR, sensitivity
        ├── cog/                     # COG enrichment + barplot
        ├── huella/                  # Top 50 fingerprint + heatmaps
        ├── filogenia/               # Core-SNP tree, Gubbins output
        └── colocalization/          # Co-localization of pirAB + neighbors

## Requirements

| Software | Version |
|----------|---------|
| Prokka | 1.14.6 |
| PIRATE | 1.0.4 |
| Gubbins | 2.4.1 |
| IQ-TREE | 3.1.3 |
| MAFFT | 7.310 |
| R | 4.3.3 |

### R packages

- ape 5.8.1
- phangorn
- treeio 1.37.0.1
- survival
- rstatix

### Python

- Python 3.7.12
- Biopython 1.81

## Execution order

1. `scripts/05_pipeline/00_prokka_annotate.sh` — Annotate genomes
2. `scripts/05_pipeline/00_pirate_pangenome.sh` — Build pangenome
3. `scripts/01_pangenome_fisher/01_pangenome_fisher.R` — Fisher's exact test
4. `scripts/01_pangenome_fisher/02_cog_analysis_roary.R` — COG enrichment
5. `scripts/01_pangenome_fisher/03_anotar_candidatos_fisher.R` — Annotate candidates
6. `scripts/02_huella_digital/04_huella_pirate_filtrada.R` — Gene fingerprint
7. `scripts/03_filogenia/05_core_snp_phylogeny.sh` — Recombination-filtered core-SNP phylogeny
8. `scripts/03_filogenia/07_ancestral_reconstruction.R` — Ancestral state reconstruction
9. `scripts/03_filogenia/08_panel_ancestral.R` — Panel of ancestral trees
10. `scripts/03_filogenia/10_curvas_acumulacion.R` — Randomized accumulation curves

## Data availability

Genome sequences:

- 11_VM: BioProject PRJNA1499920; BioSample SAMN61896695; GenBank JCBNAQ000000000
- 15_CESAIBC: BioProject PRJNA1499939; BioSample SAMN61897767; GenBank JCBNTU000000000
- 6_VM: BioProject PRJNA1499907; BioSample SAMN61896394; GenBank JCBNPX000000000

Reference genome accessions are listed in `input_data/strain_metadata/`.

## Citation

If you use this code, please cite the manuscript above.

## License

MIT — see `LICENSE`.
