library(Seurat)
library(tidyverse)
library(msigdbr)
library(stringr)
library(DESeq2)
library(openxlsx)

output_dir <- "./heart-atlas/Pseudobulk/"

pathwaysDF <- msigdbr("human", category="C2", subcategory = "CP:REACTOME")
pathway_Reactome <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

pseudobulk_DESeq2 <- function(subset, Entity, reference, compare_vs_all = FALSE) {
  # Aggregate counts
  cts <- AggregateExpression(subset, group.by = Entity, slot = "counts", return.seurat = FALSE)
  cts <- cts$RNA
  colData <- data.frame(samples = colnames(cts))
  
  # Extract "condition" column
  colData <- colData %>% mutate(condition = str_split(samples, "_", simplify = TRUE)[, 1])
  colData$condition <- sub("_[^_]+$", "", colData$samples)
  rownames(colData) <- colData$samples
  colData$samples <- NULL
  
  # Adjust condition grouping if comparing one condition vs all
  if (compare_vs_all) {
    colData$condition <- ifelse(colData$condition == reference, reference, "Other")
  }
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = colData, design = ~condition)
  dds$condition <- relevel(dds$condition, ref = reference)  # Set reference condition
  if (compare_vs_all) {
    dds$condition <- relevel(dds$condition, ref = "Other")  
  }
  
  
  # Filter low-count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  
  # Extract results
  result_col <- resultsNames(dds)[2]
  
  res <- results(dds, name = result_col)
  res <- data.frame(res)
  res$Gene <- rownames(res)
  
  return(list(results = res))
}

scRNA <- readRDS("/data/uv525372/heart-atlas/multinichenet/Heart_atlas_reduced.rds")
scRNA <- subset(scRNA, subset = (Annotation_1 != "Ischemic cells (MI)"))
scRNA <- subset(scRNA, subset = (Condition != "Ischemic"))

wb <- createWorkbook()
#cell_types <- unique(scRNA$Annotation_1)
cell_types <- c("Lymphocytes", "Cardiomyocytes", "Endocardial", "Endothelial", "Fibroblast", "Myeloid","Pericytes","VSMC")

for(celltype in cell_types){
  tmp <- subset(scRNA, subset = (Annotation_1 == celltype))
  print(celltype)
  for(subcluster in unique(tmp$Subclustering_2)){
    print(subcluster)
    results <- pseudobulk_DESeq2(tmp, c("Subclustering_2", "Donor"), subcluster, compare_vs_all=TRUE)
    
    addWorksheet(wb, sheetName = substr(subcluster, 1, 31))
    writeData(wb, sheet = substr(subcluster, 1, 31), results$results %>% arrange(desc(log2FoldChange)))
    
  }
}
saveWorkbook(wb, file = paste0(output_dir, "DEG_Subcluster.xlsx"), overwrite = TRUE)