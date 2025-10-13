library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(qs)
library(openxlsx)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(biomaRt)
library(stringr)
library(EnhancedVolcano)
library(msigdbr)
library(fgsea)

OrgDb <- org.Hs.eg.db
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = "99")

pathwaysDF <- msigdbr("human", category="C2", subcategory = "CP:REACTOME")
pathway_Reactome <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

sample_list <- list.files(
  path = "./Cellbender/RNA",
  recursive = FALSE, 
  full.names = FALSE
)

directory_list = paste(
  sample_list, "/", sep=""
)

cell_ranger_merged <- Read10X_h5_Multi_Directory(base_path = "./Cellbender/RAW/",
                                                 default_10X_path = TRUE, h5_filename = "filtered_feature_bc_matrix.h5", sample_list = sample_list, merge = TRUE, parallel = FALSE)

cell_bender_merged <- Read_CellBender_h5_Multi_File(data_dir = "./Cellbender/RNA",
                                                    custom_name = "output_file_filtered.h5", sample_names = sample_list, sample_list = directory_list, merge = TRUE, )

dual_seurat <- Create_CellBender_Merged_Seurat(raw_cell_bender_matrix = cell_bender_merged, raw_counts_matrix = cell_ranger_merged,
                                               raw_assay_name = "RAW")


#dual_seurat <- readRDS("./Cellbender/RNA/CellBender_dual_seurat_MI.rds")
Idents(dual_seurat) <- "Annotation_1"

# Cell types with low abundancy and quality
dual_seurat <- subset(dual_seurat, subset = (Subclustering != "EC, IFN"))
dual_seurat <- subset(dual_seurat, subset = (Subclustering != "Discard"))
dual_seurat <- subset(dual_seurat, subset = (Subclustering != "Fibroblast (APOD/E)"))
dual_seurat <- subset(dual_seurat, subset = (Subclustering != "Cardiomyocytes (LINGO1)"))
print(table(dual_seurat$Condition))

wb <- createWorkbook()
pdf("Volcano_plot_Ischemic_Cluster_Counts.pdf")

RAW_marker <- FindAllMarkers(dual_seurat, logfc.threshold = 0.25, only.pos = TRUE, assay = "RAW")
addWorksheet(wb, "Main_Cell_Types_Counts")
writeData(wb, sheet = "Main_Cell_Types_Counts", RAW_marker)

RNA_marker <- FindAllMarkers(dual_seurat, logfc.threshold = 0.25, only.pos = TRUE, assay = "RNA")
addWorksheet(wb, "Main_Cell_Types_CellBender")
writeData(wb, sheet = "Main_Cell_Types_CellBender", RNA_marker)


ischemic_cluster <- subset(dual_seurat, subset = (Annotation_1 == "Ischemic cells (MI)"))
Idents(ischemic_cluster) <- "Subclustering"

RAW_marker <- FindAllMarkers(ischemic_cluster, logfc.threshold = 0.15, only.pos = TRUE, assay = "RAW")
addWorksheet(wb, "Ischemic_Subcluster_Counts")
writeData(wb, sheet = "Ischemic_Subcluster_Counts", RAW_marker)

RNA_marker <- FindAllMarkers(ischemic_cluster , logfc.threshold = 0.15, only.pos = TRUE, assay = "RNA")
addWorksheet(wb, "Ischemic_Subcluster_CellBender")
writeData(wb, sheet = "Ischemic_Subcluster_CellBender", RNA_marker)


Cell_types <- c("Endothelial", "Fibroblast", "Cardiomyocytes", "Myeloid")

for(i in 1:length(Cell_types)){
  
  if(Cell_types[i] == "Endothelial"){
    subcluster <- subset(dual_seurat, subset = (Annotation_1 == "Endothelial" | Subclustering == "Ischemic Endothelial cells" ))
  }
  
  if(Cell_types[i] == "Fibroblast"){
    subcluster <- subset(dual_seurat, subset = (Annotation_1 == "Fibroblast" | Subclustering == "Ischemic Fibroblasts" ))
  }
  

  if(Cell_types[i] == "Cardiomyocytes"){
    subcluster <- subset(dual_seurat, subset = (Annotation_1 == "Cardiomyocytes" | Subclustering == "Ischemic Cardiomyocytes #1" | Subclustering == "Ischemic Cardiomyocytes #2" | Subclustering == "Ischemic Cardiomyocytes #3" ))
  }
  if(Cell_types[i] == "Myeloid"){
    subcluster <- subset(dual_seurat, subset = (Annotation_1 == "Monocytes" | Subclustering == "Ischemic Mononuclear Phagocytes"))
  }
  
  print(table(subcluster$Subclustering))
  
  
  marker_RAW <- FindMarkers(subcluster, ident.1 = "Ischemic cells (MI)", assay = "RAW")
  marker_RAW$Gene <- rownames(marker_RAW)
  
  marker_RNA <- FindMarkers(subcluster, ident.1 = "Ischemic cells (MI)", assay = "RNA")
  marker_RNA$Gene <- rownames(marker_RNA)
  
  ###
  gene_annotations <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "entrezgene_accession", "entrezgene_id"),
                            filters = "external_gene_name",
                            values = marker_RAW$Gene,
                            mart = mart)
  tmp <- merge(marker_RAW, gene_annotations, by.x = "Gene", by.y = "external_gene_name", all.x = TRUE)
  tmp <- tmp[!is.na(tmp$avg_log2FC) & !is.na(tmp$entrezgene_id), ]
  tmp <- tmp[tmp$pct.1 > 0.1 | tmp$pct.2 > 0.1, ]
  tmp <- tmp %>% arrange(desc(avg_log2FC))
  geneList <- setNames(as.vector(tmp$avg_log2FC), tmp$entrezgene_id)
  
  react <- gsePathway(
    geneList = geneList,
    minGSSize = 15,
    maxGSSize = 500,
    pvalueCutoff = 0.5,
    pAdjustMethod = "BH", 
    verbose = FALSE
  )
  
  react_readable <- setReadable(react, 'org.Hs.eg.db', 'ENTREZID')
  
  
  biotype <- gene_annotations %>% dplyr::select(external_gene_name, ensembl_gene_id, gene_biotype)
  print(head(biotype))
  marker_RAW <- marker_RAW %>%
    left_join(biotype, by = c("Gene" = "external_gene_name"))

  
  addWorksheet(wb, paste(Cell_types[i], "_Counts"))
  writeData(wb, sheet = paste(Cell_types[i], "_Counts"), marker_RAW)
  
  addWorksheet(wb, paste(Cell_types[i], "_CellBender"))
  writeData(wb, sheet = paste(Cell_types[i], "_CellBender"), marker_RNA)
  
  addWorksheet(wb, paste("CP_Reactome ", Cell_types[i]))
  writeData(wb, sheet = paste("CP_Reactome ", Cell_types[i]), react_readable@result)
  
  print(head(tmp))
  ranks <- tmp$avg_log2FC
  names(ranks) <- tmp$Gene
  ranks <- na.omit(ranks)
  
  fgseaRes <- fgsea(pathways = pathway_Reactome, stats = ranks, minSize = 15, maxSize = 500)
  addWorksheet(wb, paste0("Reactome ", Cell_types[i]))
  writeData(wb, sheet = paste0("Reactome ", Cell_types[i]), fgseaRes)
  
  
  if(Cell_types[i] == "Cardiomyocytes"){
    shape_map <- c(
    "protein_coding" = 17,
    "lncRNA" = 3)
  
    keyvals.shape <- shape_map[marker_RAW$gene_biotype]
    keyvals.shape[!marker_RAW$gene_biotype %in% names(shape_map)] <- NA  # Set all others to NA
    names(keyvals.shape) <- ifelse(marker_RAW$gene_biotype %in% names(shape_map),
                                 marker_RAW$gene_biotype, NA)
  
    selectGenes <- marker_RAW$Gene[marker_RAW$gene_biotype %in% c('protein_coding', 'lncRNA')]
  
    print(EnhancedVolcano(marker_RAW,
                  lab = marker_RAW$Gene,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  title = Cell_types[i],labSize = 4.0,
                  subtitle = "Ischemic zone vs control", col=c('gray', 'gray', 'gray', 'red3'),
                  shapeCustom = keyvals.shape, legendPosition = 'bottom',
                  legendLabSize = 10,
                  legendIconSize = 4,)+coord_flip())
    }
  
}

saveWorkbook(wb, file = "DEG_Ambient_RNA_Ischemic_Cluster.xlsx", overwrite = TRUE)
dev.off()


### Combined plotting of pathways by cell type
pdf("Reactome_pathways_Ischemic_Cells.pdf")
Cell_types <- c("Endothelial", "Fibroblast", "Cardiomyocytes", "Myeloid")
all_filtered <- list()
pathway_tracker <- list()

# Load and filter all sheets first
for (c in seq_along(Cell_types)) {
  cell <- Cell_types[c]
  path <- read_excel("DEG_Ambient_RNA_Ischemic_Cluster.xlsx", sheet = paste0("CP_Reactome  ", cell))
  #path <- path %>% filter(p.adjust < 0.05)
  path$cell_type <- cell
  all_filtered[[cell]] <- path
  pathway_tracker[[cell]] <- unique(path$Description)
}

# Count in how many cell types each pathway appears
pathway_counts <- bind_rows(all_filtered) %>%
  distinct(Description, cell_type) %>%
  group_by(Description) %>%
  summarise(celltype_count = n(), .groups = "drop")

# Filter overlaps within each cell type
final_filtered <- list()

for (cell in Cell_types) {
  path <- all_filtered[[cell]]
  gene_lists <- lapply(path$core_enrichment, function(x) unlist(strsplit(x, "/")))
  n <- length(gene_lists)
  keep <- rep(TRUE, n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      genes_i <- gene_lists[[i]]
      genes_j <- gene_lists[[j]]
      
      overlap_i <- length(intersect(genes_i, genes_j)) / length(genes_i)
      overlap_j <- length(intersect(genes_i, genes_j)) / length(genes_j)
      
      if (overlap_i >= 0.75 || overlap_j >= 0.75) {
        desc_i <- path$Description[i]
        desc_j <- path$Description[j]
        
        # Look up how many cell types they appear in
        ct_i <- pathway_counts$celltype_count[pathway_counts$Description == desc_i]
        ct_j <- pathway_counts$celltype_count[pathway_counts$Description == desc_j]
        
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

# Combine the final results
final_result <- bind_rows(final_filtered)


# Get union of top 15 genes from each cell type
top_pathways <- final_result %>%
  filter(Description != "Innate Immune System") %>%
  group_by(cell_type) %>%
  slice_min(order_by = p.adjust, n = 35, with_ties = FALSE) %>%
  pull(Description) %>%
  unique()

# Keep those genes across all cell types
pathways_filtered <- final_result %>%
  filter(Description %in% top_pathways)

print(ggplot(pathways_filtered, aes(x=cell_type, y=reorder(Description, NES), size=p.adjust, color=NES))+geom_point()+
        labs(y="", x="", size = "adj. p_val")+
        scale_color_distiller(palette = "RdBu", direction = 1) +
        theme_minimal() +
        scale_size(trans = "reverse", range = c(0, 4))+ 
        theme(
          #plot.background = element_rect(fill = "transparent"),
          #panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          axis.title = element_text(size=14),
          axis.text = element_text(size=11),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=12),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
        ))
dev.off()


### Plotting of pathway with hierarchical clustering 
pdf(paste0("Reactome_pathways_Ischemic_Cells2.pdf"), width = 8, height = 11)
Cell_types <- c("Endothelial", "Fibroblast", "Cardiomyocytes", "Myeloid")
all_filtered <- list()
pathway_tracker <- list()

# Load and filter all sheets first
for (c in seq_along(Cell_types)) {
  cell <- Cell_types[c]
  path <- read_excel("DEG_Ambient_RNA_Ischemic_Cluster.xlsx", sheet = paste0("Reactome ", cell))
  #path <- path %>% filter(padj < 0.05)
  path$cell_type <- cell
  all_filtered[[cell]] <- path
  pathway_tracker[[cell]] <- unique(path$pathway)
}

# Count in how many cell types each pathway appears
pathway_counts <- bind_rows(all_filtered) %>%
  distinct(pathway, cell_type) %>%
  group_by(pathway) %>%
  summarise(celltype_count = n(), .groups = "drop")

# Filter overlaps within each cell type
final_filtered <- list()

for (cell in Cell_types) {
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

# Combine the final results
final_result <- bind_rows(final_filtered)


# Get union of top 10 genes from each cell type
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


# Keep those genes across all cell types
pathways_filtered <- final_result %>%
  filter(pathway %in% top_and_bottom_pathways)


pathway_mat <- pathways_filtered %>%
  dplyr::select(pathway, cell_type, NES) %>%
  tidyr::pivot_wider(names_from = cell_type, values_from = NES, values_fill = 0) %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()

# Calculate distance and clustering
hc <- hclust(dist(pathway_mat))

# Cut the dendrogram into k clusters
k <- 8
pathway_clusters <- cutree(hc, k = k)

# Store as data frame for merging
cluster_df <- data.frame(pathway = names(pathway_clusters), cluster = factor(pathway_clusters))

# Optional: reorder levels to match dendrogram order
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

# Apply to pathways that actually contain the gene
pathways_filtered$pathway <- sapply(pathways_filtered$pathway, truncate_pathway)

# Flag pathways that contained the gene once
pathways_filtered <- pathways_filtered %>%
  group_by(pathway) %>%
  mutate(
    mark_star = any(grepl("STAB1|SERPINE1", leadingEdge))
  ) %>%
  ungroup() %>%
  mutate(
    pathway = ifelse(mark_star, paste0("*", pathway), pathway)
  )

ggplot(pathways_filtered, aes(x = cell_type, y = pathway, size = -log10(padj), color = NES)) +
  geom_point() +
  labs(y = "", x = "", size = "-log10(adj. p-value)") +
  scale_color_distiller(palette = "RdBu", direction = 1, guide = guide_colorbar(
    title.position = "top", title.hjust = 0.5, reverse = TRUE)) +
  scale_size(
    range = c(1, 4),
    guide = guide_legend(
      title.position = "top",
      title.hjust = 0.5, reverse = TRUE))+
  facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    strip.text.y = element_text(angle = 0),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
  )
dev.off()