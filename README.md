# Building of the human single-nuclei heart failure atlas
Heart failure encompasses a diverse group of cardiomyopathies, including myocardial infarction, hypertrophic, dilated, and arrhythmogenic forms, each defined by distinct etiologies. By integrating single-nuclei transcriptomic data from human heart tissue across these conditions, we constructed a unified atlas containing 1.8 million nuclei from 195 individuals.
The atlas reveals disease-specific cellular clusters, transcriptional changes and altered ligand-receptor interactions. We found cell states specific to the ischemic zone of myocardial infarction, explored the influence of cytokines on fibroblast cell-states, and identified etiological pathways in different cell types. The integration of summary statistics of 52 genome-wide association studies with atlas-wide gene expression highlighted genetically associated pathways involving metabolic dysregulation and ion channel dysfunction.
Implementation of an LLM-based workflow led to the identification of ZLN005, a small molecule that boosts mitochondrial biogenesis via PGC-1α, whose cardioprotective effect we validated experimentally, underscoring the utility of the atlas in early therapeutic target discovery for heart failure.

<center>
  <img src="https://github.com/ToreBle/Heart_Failure_Atlas/blob/main/figure/Overview_Human_Heart_Atlas.png" width="900">
  <figcaption>Figure 1: Overview of the building of the human heart atlas and subsequent analysis including a UMAP of the main cell types from human heart tissue and an schematic overview of the calcium handling within a cardiomyocyte.</figcaption>
</center>

<br>
<br>

## Workflow of the multi-modal human heart atlas:
The analysis of the integrated single-cell RNA-seq data and GWAS were performed as described in the following sections: 

-	01 – QC and snRNA-seq data integration
-	02 – Ischemic cluster
-	03 – Cell- and patient clustering

-	04 – CCI and bulkRNA-seq analysis
-	05 – Reactome pathway analysis
-	06 – GWAS analysis with scPagwas
- 07 - Setup of the LLM workflow for drug target discovery

## Data
- Raw sequencing data: [GSE350002] (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE350002)
- Human heart failure atlas: [h5ad](https://doi.org/10.5281/zenodo.23216491)

## Reference
The manuscript has been submitted for peer-review. A preprint is available on [medRxiv](https://www.medrxiv.org/content/10.1101/2025.10.09.25337646v1).
