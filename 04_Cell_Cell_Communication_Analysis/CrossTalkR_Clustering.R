library(CrossTalkeR)
library(igraph)
library(ggraph)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(dplyr)
library(proxy)
library(RColorBrewer)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)

setwd("./CCI/")

# Define paths and dataset names
dataset_paths <- list(
  Chaffin_DCM   = "Chaffin_DCM_LR_data_final.Rds",
  Chaffin_HCM   = "Chaffin_HCM_LR_data_final.Rds",
  Simonson_PostMI = "Simonson_PostMI_LR_data_final.Rds",
  Kuppe_Ischemic = "Kuppe_Ischemic_LR_data_final.Rds"
)

Cell_type_colors <- c("Adipocytes" = '#1f77b4', "Cardiomyocytes" = '#ff7f0e', "Endocardial" = '#279e68', "Endothelial" = '#d62728', 
                      "Epicardium" = '#aa40fc', "Fibroblast" = '#8c564b',"Ischemic cells (MI)" = '#e377c2', "Lymphatic EC" = '#b5bd61', 
                      "Lymphocytes" = '#17becf', "Mast cells" = '#aec7e8', "Myeloid" = '#ffbb78', "Neuronal" = '#98df8a',
                      "Pericytes" = '#ff9896', "VSMC" = '#c5b0d5')


cci_data <- list()
table_list <- list()

# Load all RDS files and assign colors
for (name in names(dataset_paths)) {
  cat("Loading", name, "\n")
  obj <- readRDS(dataset_paths[[name]])
  obj@colors <- Cell_type_colors
  cci_data[[name]] <- obj
}

# Plot only graphs named 'Disease_X_Control'
pdf("CrosstalkR.pdf")

for (name in names(cci_data)) {
  cci <- cci_data[[name]]
  
  # Only keep graphs that match the pattern "Disease_x_Control"
  graph_names <- names(cci@graphs)
  graph_names <- grep("_x_Control$", graph_names, value = TRUE)  # <- filter here
  
  for (gname in graph_names) {
    plot_cci(
      graph = cci@graphs[[gname]],
      colors = cci@colors,
      plt_name = paste(gname, "-", name),
      coords = cci@coords[V(cci@graphs[[gname]])$name, ],
      col_pallet = c("#fdbb2d", "white", "#1a2a6c"),
      emax = NULL,
      leg = TRUE,
      low = 0,
      high = 0,
      ignore_alpha = FALSE,
      log = FALSE,
      efactor = 8,
      vfactor = 12,
      pg = cci@rankings[[gname]]$Pagerank
    )
  }
}
dev.off()


###
cluster_crosstalk_matrix <- function(cci_data, k = 8,
                                     include_unclustered = FALSE,
                                     zscore_by = "column",
                                     distance_method = "spearman") {
  require(dplyr)
  require(reshape2)
  require(matrixStats)
  require(proxy)
  
  table_list <- list()
  
  for (name in names(cci_data)) {
    cci <- cci_data[[name]]
    table_names <- grep("_x_Control$", names(cci@tables), value = TRUE)
    
    for (tname in table_names) {
      disease_label <- gsub(".*_", "", name)
      tbl <- cci@tables[[tname]] %>% mutate(Condition = disease_label)
      table_list[[paste(name, tname, sep = "_")]] <- tbl
    }
  }
  
  # Combine all into one data frame
  table_all <- bind_rows(table_list)
  
  # Create wide matrix
  mat <- table_all %>%
    group_by(allpair) %>%
    filter(any(LRScore > 0)) %>%
    ungroup() %>%
    select(allpair, Condition, LRScore) %>%
    distinct() %>%
    reshape2::dcast(allpair ~ Condition, value.var = "LRScore", fill = 0)
  
  rownames(mat) <- mat$allpair
  mat <- as.matrix(mat[, -1])
  
  # Z-score scaling
  if (tolower(zscore_by) == "row") {
    mat_scaled <- t(scale(t(mat)))  # Z-score by row
  } else {
    mat_scaled <- scale(mat)        # Z-score by column (default)
  }
  
  # Distance calculation
  if (tolower(distance_method) == "cosine") {
    dist_mat <- proxy::dist(mat_scaled, method = "cosine")
  } else if (tolower(distance_method) == "spearman") {
    dist_mat <- as.dist(1 - cor(t(mat_scaled), method = "spearman"))
  } else {
    stop("Unsupported distance_method. Use 'spearman' or 'cosine'.")
  }
  
  # Hierarchical clustering
  hc <- hclust(dist_mat)
  clusters <- cutree(hc, k = k)
  
  # Join cluster info
  cluster_df <- data.frame(allpair = names(clusters), cluster = factor(clusters))
  
  if (!include_unclustered) {
    clustered_crosstalkr_table <- left_join(table_all, cluster_df, by = "allpair") %>%
      filter(!is.na(cluster))
  } else {
    clustered_crosstalkr_table <- left_join(table_all, cluster_df, by = "allpair")
  }
  
  # Order matrix by cluster assignment
  mat_ordered <- mat_scaled[names(clusters)[order(clusters)], ]
  
  return(list(
    clustered_crosstalkr_table = clustered_crosstalkr_table,
    mat = mat_ordered,
    clusters = clusters
  ))
}


res <- cluster_crosstalk_matrix(
  cci_data,
  k = 8,
  include_unclustered = FALSE,
  zscore_by = "column",
  distance_method = "spearman"
)

clustered_df <- res$clustered_crosstalkr_table
mat <- res$mat
clusters <- res$clusters

###

pdf("CrossTalkR_downstream_anaylsis.pdf", width = 3, height = 5)
# Cluster annotation
cluster_levels <- sort(unique(clusters))
cluster_colors <- setNames(
  brewer.pal(n = length(cluster_levels), name = "Paired"),
  cluster_levels
)

annotation_row <- data.frame(cluster = factor(clusters, levels = cluster_levels))
rownames(annotation_row) <- names(clusters)

ann_colors <- list(cluster = cluster_colors)

# Plot
pheatmap(
  mat_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#fdbb2d", "white", "#1a2a6c"))(100),
  legend = TRUE
)



# Transpose matrix: now rows = conditions, columns = allpairs
mat_condition <- t(mat)

# PCA on conditions
pca_cond <- prcomp(mat_condition, center = TRUE, scale. = TRUE)
pca_df_cond <- as.data.frame(pca_cond$x)
pca_df_cond$Condition <- rownames(pca_df_cond)

# Variance explained
percent_var_cond <- round(100 * pca_cond$sdev^2 / sum(pca_cond$sdev^2), 1)

# Plot PCA
ggplot(pca_df_cond, aes(PC1, PC2, label = Condition)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  labs(
    title = "PCA of Conditions by LRScore",
    x = paste0("PC1 (", percent_var_cond[1], "%)"),
    y = paste0("PC2 (", percent_var_cond[2], "%)")
  ) +
  theme_minimal()
dev.off()



clustered_df <- clustered_df %>%
  separate(allpair, into = c("ligand_part", "receptor_part"), sep = "@", remove = FALSE) %>%
  mutate(
    ligand_gene = sub(".*/", "", sub("\\|L", "", ligand_part)),
    receptor_gene = sub(".*/", "", sub("\\|R", "", receptor_part)),
    ligand_receptor_pair = paste(ligand_gene, receptor_gene, sep = "→")
  )


top_LR_per_cluster <- clustered_df %>%
  group_by(cluster, ligand_receptor_pair) %>%
  summarise(mean_score = mean(LRScore, na.rm = TRUE), .groups = "drop") %>%
  group_by(cluster) #%>%
  #slice_max(mean_score, n = 10) 

unique_LR <- top_LR_per_cluster %>%
  group_by(ligand_receptor_pair) %>%
  summarise(n_clusters = n_distinct(cluster), .groups = "drop") %>%
  filter(n_clusters == 1)

# Get the actual unique pairs with their assigned cluster
unique_LR_per_cluster <- top_LR_per_cluster %>%
  semi_join(unique_LR, by = "ligand_receptor_pair")


# Filter and transform
top_LRs <- unique_LR_per_cluster %>%
  group_by(cluster) %>%
  slice_max(order_by = mean_score, n = 5) %>%
  ungroup() %>%
  filter(mean_score > 0.02) #%>%
  #mutate(label_wrapped = gsub("→", "→\n", ligand_receptor_pair))

# Generate Paired color palette
n_clusters <- length(unique(top_LRs$cluster))
cluster_colors <- setNames(brewer.pal(n = max(n_clusters, 3), name = "Paired")[1:n_clusters],
                           sort(unique(top_LRs$cluster)))

pdf("top_and_unique_cluster_LRscores.pdf", width = 8, height = 11)
ggplot(top_LRs, aes(x = mean_score, y = reorder(ligand_receptor_pair, mean_score), fill = cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  labs(
    x = "Mean LRScore",
    y = "Ligand\n→Receptor",
    fill = "Cluster"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = ""
  )


top_LR_gene_pairs <- clustered_df %>%
  group_by(cluster, ligand_receptor_pair) %>%
  summarise(mean_score = mean(LRScore, na.rm = TRUE), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = mean_score, n = 10) %>%
  mutate(label_wrapped = gsub("→", "\n→", ligand_receptor_pair))

ggplot(top_LR_gene_pairs, aes(x = reorder(label_wrapped, mean_score), y = mean_score, fill = cluster)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(~ cluster, scales = "free_y") +
  scale_fill_brewer(palette = "Paired")+
  coord_flip() +
  theme(axis.text.y = element_text(size = 12), strip.text = element_text(size = 10))
  labs(
    x = "Ligand → Receptor",
    y = "Mean LRScore"
  ) +
  theme_minimal(base_size = 11)
dev.off()