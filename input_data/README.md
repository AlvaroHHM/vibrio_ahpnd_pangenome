# Input data for analysis scripts

This folder contains all input files required to run the analysis scripts (01–10).

## Structure

- `pirate/` — Main PIRATE output files: `PIRATE.gene_families.ordered.tsv`, `phenotype.csv`, and `feature_sequences/`.
- `pangenome/` — Roary pangenome files for the three Mexican strains: exclusive and absent gene lists (TSV and FASTA).
- `phylogeny/` — Core-genome alignment and core-SNP tree: `core_alignment.fasta`, `core_snp.treefile`.
- `eggnog/` — eggNOG annotation tables and Fisher results for annotation and sensitivity analyses.

## Notes

- The `feature_sequences/` folder contains one FASTA file per PIRATE gene family.
- All files were automatically located and copied from the original project directory.
- If a file is missing, the script will report it during execution.

## Usage

Set the `PROJECT_DIR` environment variable to the root of your project before running the scripts:

```bash
export PROJECT_DIR=/path/to/project
Rscript scripts/01_pangenome_fisher.R
