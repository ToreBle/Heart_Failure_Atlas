library(Seurat)
library(tidyverse)
library(fgsea)
library(msigdbr)
library(stringr)
library(DESeq2)
library(EnhancedVolcano)
library(openxlsx)
library(biomaRt)
library(clusterProfiler)
library(ReactomePA)

output_dir <- "./heart-atlas/Pseudobulk/"

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
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
  #if (compare_vs_all) {
  #  dds$condition <- relevel(dds$condition, ref = "Other")  
  #}
  
  
  # Filter low-count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  
  # Extract results
  result_col <- resultsNames(dds)[2]
  
  res <- results(dds, name = result_col)
  res <- data.frame(res)
  res$Gene <- rownames(res)
  
  #results <- results %>% arrange(desc(log2FoldChange))
  
  # Get gene annotations
  gene_annotations <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_accession", "entrezgene_id"),
                            filters = "external_gene_name",
                            values = res$Gene,
                            mart = mart)
  tmp <- merge(res, gene_annotations, by.x = "Gene", by.y = "external_gene_name", all.x = TRUE)
  
  # Filter out NA values
  tmp <- tmp[!is.na(tmp$log2FoldChange) & !is.na(tmp$entrezgene_id), ]
  tmp <- tmp %>% arrange(desc(log2FoldChange))
  
  # Prepare gene list for pathway analysis
  geneList <- setNames(as.vector(tmp$log2FoldChange), tmp$entrezgene_id)
  
  return(list(results = res, geneList = geneList))
}

scRNA <- readRDS("/data/uv525372/heart-atlas/multinichenet/Heart_atlas_reduced.rds")
scRNA <- subset(scRNA, subset = (Annotation_1 != "Ischemic cells (MI)"))
scRNA <- subset(scRNA, subset = (Condition != "Ischemic"))

wb <- createWorkbook()

cell_types <- unique(scRNA$Annotation_1)

for(celltype in cell_types){
  tmp <- subset(scRNA, subset = (Annotation_1 == celltype))
  geneLists <- list()
  
  for(Condition in c("MI", "HCM", "DCM")){
    if(Condition == "HCM"){
      HCM <- subset(tmp, subset = (Condition %in% c("HCM","Control")))
      results <- pseudobulk_DESeq2(HCM, c("Condition", "Donor"), "Control")}
    if(Condition == "DCM"){
      DCM <- subset(tmp, subset = (Condition %in% c("DCM","Control")))
      results <- pseudobulk_DESeq2(DCM, c("Condition", "Donor"), "Control")}
    if(Condition == "MI"){
      Ischemic <- subset(tmp, subset = (Condition %in% c("Post MI","Control")))
      results <- pseudobulk_DESeq2(Ischemic, c("Condition", "Donor"), "Control")
    }
    
    geneLists[[Condition]] <- results$geneList
    print(head(results$geneList))
    react_analysis <- gsePathway(results$geneList, pvalueCutoff = 0.2, pAdjustMethod = "BH", verbose = FALSE)
    
    # Save to Excel
    addWorksheet(wb, paste0("Reactome ", celltype, " - ", Condition))
    writeData(wb, sheet = paste0("Reactome ", celltype, " - ", Condition), react_analysis@result)
    
    #FGSEA
    #tmp <- results$results %>% arrange(desc(log2FoldChange)) %>% filter(!str_starts(Gene, "LINC")) %>% filter(!str_starts(Gene, "AC")) %>% dplyr::select(Gene, log2FoldChange)
    #ranks <- tmp$log2FoldChange
    #names(ranks) <- tmp$Gene
    #ranks <- na.omit(ranks)
    
    #fgseaRes <- fgsea(pathways = pathway_Reactome, stats = ranks, minSize = 15, maxSize = 500)
    #addWorksheet(wb, paste0("Reactome ", Condition))
    #writeData(wb, sheet = paste0("Reactome ", Condition), fgseaRes)
    
    
    addWorksheet(wb, paste0("DEG ", celltype, " - ", Condition))
    writeData(wb, sheet = paste0("DEG ", celltype, " - ", Condition), results$results %>% arrange(desc(log2FoldChange)))
    
  }
}
saveWorkbook(wb, file = paste0(output_dir, "DEG_cell_types_heart_failure.xlsx"), overwrite = TRUE)