library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(DOSE)

library(ggplot2)
library(enrichplot)
library(openxlsx)

OrgDb <- org.Hs.eg.db
Conditions <- c("Hypoxia", "IL1b", "TGFb")
#Conditions <- c("Hypoxia", "Doxorubicin")

geneLists <- list()

for (condition in Conditions) {
  # Load results for the current condition
  results <- readRDS(paste0(condition, "_res.rds"))
  
  # Filter out rows with NA values for log2FoldChange and entrezgene_id
  results <- results[!is.na(results$log2FoldChange), ]
  results <- results[!is.na(results$entrezgene_id), ]
  
  # Create gene list for the current condition, sorted by log2FoldChange
  results <- results %>% arrange(desc(log2FoldChange))
  geneList <- as.vector(results$log2FoldChange)
  names(geneList) <- results$entrezgene_id
  
  # Store all Conditions in the geneList
  geneLists[[condition]] <- geneList
}


for (condition in Conditions){
  pdf(paste0(condition, "_pathway_analysis.pdf"))
  wb <- createWorkbook()
  
  # GO analysis
  go <- gseGO(
    geneList     = geneLists[[condition]],
    OrgDb        = org.Hs.eg.db,
    ont          = "CC",
    minGSSize    = 100,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    verbose      = FALSE
  )

  go_readable <- setReadable(go, 'org.Hs.eg.db', 'ENTREZID')
  addWorksheet(wb, "Gene Ontology")
  writeData(wb, sheet = "Gene Ontology", go_readable@result)
  
  print(dotplot(go, showCategory = 20) + ggtitle(paste("GO dotplot for ", condition)))
  print(heatplot(go_readable, foldChange = geneLists[[condition]], showCategory = 3) +
    coord_flip() +
    theme(axis.text.y = element_text(size = 3), axis.text.x = element_text(size = 8)) +
    ggtitle(paste("GO heatplot for ", condition)))
  go_pdo <- pairwise_termsim(go)
  print(emapplot(go_pdo) + ggtitle(paste("GO for ", condition)))
  
  # Reactome pathway analysis
  react <- gsePathway(
    geneList = geneLists[[condition]], 
    pvalueCutoff = 0.2,
    pAdjustMethod = "BH", 
    verbose = FALSE
  )
  
  react_readable <- setReadable(react, 'org.Hs.eg.db', 'ENTREZID')
  addWorksheet(wb, "Reactome pathways")
  writeData(wb, sheet = "Reactome pathways", react_readable@result)
  
  print(dotplot(react, showCategory = 20) + ggtitle(paste("Reactome pathways for ", condition)))
  print(heatplot(react_readable, foldChange = geneLists[[condition]], showCategory = 3) +
    coord_flip() +
    theme(axis.text.y = element_text(size = 3), axis.text.x = element_text(size = 8)) +
    ggtitle(paste("Reactome heatplot for ", condition)))
  react_pdo <- pairwise_termsim(react)
  print(emapplot(react_pdo) + ggtitle(paste("Reactome pathways for ", condition)))
  
  # Disease Ontology (DO) analysis
  do <- gseDO(
    geneList       = geneLists[[condition]],
    minGSSize      = 120,
    pvalueCutoff   = 0.2,
    pAdjustMethod  = "BH",
    verbose        = FALSE
  )

  do_readable <- setReadable(do, 'org.Hs.eg.db', 'ENTREZID')
  addWorksheet(wb, "Disease Ontology")
  writeData(wb, sheet = "Disease Ontology", do_readable@result)

  print(dotplot(do, showCategory = 20) + ggtitle(paste("DO dotplot for ", condition)))
  print(heatplot(do_readable, foldChange = geneLists[[condition]], showCategory = 3) +
    coord_flip() +
    theme(axis.text.y = element_text(size = 3), axis.text.x = element_text(size = 8)) +
    ggtitle(paste("DO heatmap ", condition)))
  do_pdo <- pairwise_termsim(do)
  print(emapplot(do_pdo) + ggtitle(paste("Disease ontology for ", condition)))
  
  # DisGeNET (DGN) analysis
  dgn <- gseDGN(
    geneList      = geneLists[[condition]],
    pvalueCutoff  = 0.2,
    pAdjustMethod = "BH",
    verbose       = FALSE
  )
  
  dgn_readable <- setReadable(dgn, 'org.Hs.eg.db', 'ENTREZID')
  addWorksheet(wb, "DisGeNET")
  writeData(wb, sheet = "DisGeNET", dgn_readable@result)
  
  print(dotplot(dgn, showCategory = 20) + ggtitle(paste("DGN dotplot for ", condition)))
  print(heatplot(dgn_readable, foldChange = geneLists[[condition]], showCategory = 3) +
    coord_flip() +
    theme(axis.text.y = element_text(size = 3), axis.text.x = element_text(size = 8)) +
    ggtitle(paste("DGN heatplot for ", condition)))
  dgn_pdo <- pairwise_termsim(dgn)
  print(emapplot(dgn_pdo) + ggtitle(paste("DisGeNET analysis for ", condition)))
  
  saveWorkbook(wb, file = paste0(condition, "_pathway_analysis.xlsx"), overwrite = TRUE)
  dev.off()
}
