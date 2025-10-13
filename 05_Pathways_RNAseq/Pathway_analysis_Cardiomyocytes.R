library(SeuratDisk)
library(Seurat)
library(tidyverse)
library(fgsea)
library(msigdbr)
library(stringr)
library(DESeq2)
library(openxlsx)

cell_state <- "CM"
output_dir <- "./heart-atlas/Pseudobulk/Cardiomyocytes/"

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
  result_col <- resultsNames(dds)[2]  # Ensure correct contrast selection
  
  res <- results(dds, name = result_col)
  res <- data.frame(res)
  res$Gene <- rownames(res)
  
  return(list(results = res))
}

# Load the single-cell RNA-seq samples for HCM/DCM/ICM
#Convert(paste0(output_dir, "DCM_", cell_state, "_DEG.h5ad"), paste0(output_dir, "DCM_", cell_state, "_DEG.h5Seurat"), assay = "RNA", overwrite = T)
CM <- LoadH5Seurat(paste0(output_dir, "DCM_", cell_state, "_DEG.h5Seurat"), meta.data=FALSE, misc = FALSE)
meta_data_CM = read.csv(paste0(output_dir, "meta_DCM_", cell_state ,"_DEG.csv"), row.names=1)
CM <- AddMetaData(CM, metadata = meta_data_CM)

HCM <- subset(CM, subset = (Condition != "DCM"))
DCM <- subset(CM, subset = (Condition != "HCM"))

#Ischemic
#Convert(paste0(output_dir, "Ischemic_", cell_state, "_DEG.h5ad"), paste0(output_dir, "Ischemic_", cell_state, "_DEG.h5Seurat"), assay = "RNA", overwrite = T)
Ischemic <- LoadH5Seurat(paste0(output_dir, "Ischemic_", cell_state ,"_DEG.h5Seurat"), meta.data=FALSE, misc = FALSE)
meta_data_Ischemic = read.csv(paste0(output_dir,"meta_Ischemic_", cell_state ,"_DEG.csv"), row.names=1)
Ischemic <- AddMetaData(Ischemic, metadata = meta_data_Ischemic)
Ischemic <- subset(Ischemic, subset = Study != "Simonson 2023")

wb <- createWorkbook()

for(Condition in c("MI", "HCM", "DCM")){
  if(Condition == "HCM"){
    results <- pseudobulk_DESeq2(HCM, c("Condition", "Donor"), "Control")}
  if(Condition == "DCM"){
    results <- pseudobulk_DESeq2(DCM, c("Condition", "Donor"), "Control")}
  if(Condition == "MI"){
    results <- pseudobulk_DESeq2(Ischemic, c("Condition", "Donor"), "Control")
  }
  
  tmp <- results$results %>% arrange(desc(log2FoldChange)) %>% filter(!str_starts(Gene, "LINC")) %>% filter(!str_starts(Gene, "AC")) %>% dplyr::select(Gene, log2FoldChange)
  
  ranks <- tmp$log2FoldChange
  names(ranks) <- tmp$Gene
  ranks <- na.omit(ranks)
  
  fgseaRes <- fgsea(pathways = pathway_Reactome, stats = ranks, minSize = 15, maxSize = 500)
  addWorksheet(wb, paste0("Reactome ", Condition))
  writeData(wb, sheet = paste0("Reactome ", Condition), fgseaRes)
  
  
  addWorksheet(wb, paste0("DEG pseudobulk ", Condition))
  writeData(wb, sheet = paste0("DEG pseudobulk ", Condition), results$results %>% arrange(desc(log2FoldChange)))
  
}

saveWorkbook(wb, file = paste0(output_dir, "Main_", cell_state ,"_CTRL.xlsx"), overwrite = TRUE)


CM_marker <- read.csv(paste0(output_dir, cell_state, "_to_all.csv"))
CM_marker <- CM_marker %>% filter(pvals_adj<0.05) %>% filter(logfoldchanges>0.5)

wb <- createWorkbook()
pdf(paste0(output_dir, "Main_", cell_state, ".pdf"), width = 15)
for(Condition in c("MI", "HCM", "DCM")){
  CM_H <- read.xlsx(paste0(output_dir, "Main_", cell_state, "_CTRL.xlsx"), sheet = paste0("DEG pseudobulk ", Condition))

CM_H$names = CM_H$Gene
tmp <- left_join(CM_marker, CM_H, by="names")

tmp <- tmp %>% arrange(desc(scores)) %>% arrange(desc(log2FoldChange)) %>% filter(!str_starts(names, "LINC")) %>% filter(!str_starts(names, "AC")) %>% dplyr::select(names, log2FoldChange)

ranks <- tmp$log2FoldChange
names(ranks) <- tmp$names ##
ranks <- na.omit(ranks)

fgseaRes <- fgsea(pathways = pathway_Reactome, stats = ranks, minSize = 15, maxSize = 500)
addWorksheet(wb, paste0("Reactome ", Condition))
writeData(wb, sheet = paste0("Reactome ", Condition), fgseaRes)


print(ggplot(fgseaRes %>% arrange(NES) %>% filter(padj<0.1) , aes(x=fct_inorder(pathway), y=padj, color=NES, size=size))+geom_point()+coord_flip()+
  scale_color_distiller(palette = "RdBu") +scale_size(range = c(0.0005, 5), breaks = c( 0.0005, 0.005, 0.05, 0.5), trans = 'reverse')+
  theme_classic()+theme(plot.title = element_text(colour = "black", size =12), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, colour = "black"), axis.text.y = element_text(size=6, colour = "black"))+
  labs(x="", y="", title = Condition))

}
saveWorkbook(wb, file = paste0(output_dir, "Pathways_", cell_state , "_CTRL.xlsx"), overwrite = TRUE)
dev.off()



pdf(paste0(output_dir,"Reactome_pathways_", cell_state, ".pdf"), width = 8)
Condition <- c("MI", "HCM", "DCM")
all_filtered <- list()
pathway_tracker <- list()

# 1. Load and filter all sheets first
for (c in seq_along(Condition)) {
  cell <- Condition[c]
  path <- read.xlsx(paste0(output_dir, "Pathways_", cell_state ,"_CTRL.xlsx"), sheet = paste0("Reactome ", cell))
  
  path$cell_type <- cell
  all_filtered[[cell]] <- path
  pathway_tracker[[cell]] <- unique(path$pathway)
}

# 2. Count in how many cell types each pathway appears
pathway_counts <- bind_rows(all_filtered) %>%
  distinct(pathway, cell_type) %>%
  group_by(pathway) %>%
  summarise(celltype_count = n(), .groups = "drop")

# 3. Filter overlaps within each cell type
final_filtered <- list()

for (cell in Condition) {
  path <- all_filtered[[cell]]
  gene_lists <- lapply(path$leadingEdge, function(x) unlist(strsplit(x, "/")))
  n <- length(gene_lists)
  keep <- rep(TRUE, n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      genes_i <- gene_lists[[i]]
      genes_j <- gene_lists[[j]]
      
      overlap_i <- length(intersect(genes_i, genes_j)) / length(genes_i)
      overlap_j <- length(intersect(genes_i, genes_j)) / length(genes_j)
      
      if (overlap_i >= 0.75 || overlap_j >= 0.75) {
        desc_i <- path$pathway[i]
        desc_j <- path$pathway[j]
        
        # Look up how many cell types they appear in
        ct_i <- pathway_counts$celltype_count[pathway_counts$pathway == desc_i]
        ct_j <- pathway_counts$celltype_count[pathway_counts$pathway == desc_j]
        
        # Keep the one more commonly seen
        if (ct_i > ct_j) {
          keep[j] <- FALSE
        } else if (ct_j > ct_i) {
          keep[i] <- FALSE
        } else {
          # If equal, keep the one with the shorter name
          if (nchar(desc_i) <= nchar(desc_j)) {
            keep[j] <- FALSE
          } else {
            keep[i] <- FALSE
          }
        }
      }
    }
  }
  
  final_filtered[[cell]] <- path[keep, ]
}

# 4. Combine the final results
final_result <- bind_rows(final_filtered)


# Step 1: Get union of top 15 genes from each cell type
top_and_bottom_pathways <- final_result %>%
  group_by(cell_type) %>%
  filter(padj<0.05) %>%
  slice_min(order_by = NES, n = 10, with_ties = FALSE) %>%
  bind_rows(
    final_result %>%
      group_by(cell_type) %>%
      filter(padj<0.05) %>%
      slice_max(order_by = NES, n = 10, with_ties = FALSE)
  ) %>%
  pull(pathway) %>%
  unique()


# Step 2: Keep those genes across all cell types
pathways_filtered <- final_result %>%
  filter(pathway %in% top_and_bottom_pathways)


pathway_mat <- pathways_filtered %>%
  select(pathway, cell_type, NES) %>%
  tidyr::pivot_wider(names_from = cell_type, values_from = NES, values_fill = 0) %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()

# Calculate distance and set the number of clusters
hc <- hclust(dist(pathway_mat))
k <- 6
pathway_clusters <- cutree(hc, k = k)
cluster_df <- data.frame(pathway = names(pathway_clusters), cluster = factor(pathway_clusters))

# Reorder levels to match dendrogram order
ordered_pathways <- names(pathway_clusters)[hc$order]
cluster_df$pathway <- factor(cluster_df$pathway, levels = ordered_pathways)

# Merge the cluster info into the main data
pathways_filtered <- merge(pathways_filtered, cluster_df, by = "pathway")
pathways_filtered$pathway <- gsub("^REACTOME_", "", pathways_filtered$pathway)

# Define a truncation function
truncate_pathway <- function(pathway, max_length = 50) {
  if (nchar(pathway) <= max_length) {
    return(pathway)
  }
  # Find underscores within the limit
  underscores <- gregexpr("_", pathway)[[1]]
  cut_pos <- max(underscores[underscores <= max_length])
  if (is.finite(cut_pos)) {
    return(paste0(substr(pathway, 1, cut_pos - 1), "..."))
  } else {
    return(paste0(substr(pathway, 1, max_length), "..."))
  }
}

# Apply to the column
pathways_filtered$pathway <- sapply(pathways_filtered$pathway, truncate_pathway)

ggplot(pathways_filtered, aes(x = cell_type, y = pathway, size = -log10(padj), color = NES)) +
  geom_point() +
  labs(y = "", x = "", size = "adj. p-value") +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    strip.text.y = element_text(angle = 0),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

dev.off()