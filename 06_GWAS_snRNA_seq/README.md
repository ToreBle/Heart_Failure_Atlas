# 06 – GWAS integration with snRNA-seq

Scripts for the combined analysis of different genome-wide association studies (GWAS) with the snRNA-seq data from the human heart failure atlas.

## Scripts

- **`<GWAS_studies_summary>.txt`** – Overview of different GWAS studies and their accession number from NHGRI-EBI Catalog.
- **`<GWAS_summary_statistics_heart_atlas>.R`** – Meta data summary of the GWAS summary statistics.
- **`<scPagwas_main>.R`** - Main script for the analysis with scPagwas.
- **`<scPagwas_PCC_genes_traits>.R`** - PCC analysis from scPagwas.
- **`<Plotting_scPagwas_Overview>.R`** - Plotting of the heatmaps for patients and GWAS.
- **`<scPagwas_summary_by_celltype>.R`** - Summary of the GWAS by cell type.
- **`<scPagwas_cell_type_pathways>.R`** - Summary of the genetically-associated pathways by cell type.
- **`<eQTL_GWAS_mapping>.R`** - Plotting of eQTL data for selected loci with eQTpLot.
