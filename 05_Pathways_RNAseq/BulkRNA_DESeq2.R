library(Rsubread)
library(DESeq2)
library(biomaRt)

library(dplyr)
library(openxlsx)
library(optparse)

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(stringr)


#name <- "Bulk_RNA_Heart_Atlas"
#base_directory = "./heart-atlas/bulkRNA"
#meta_data_csv = "BulkRNA_samples_Heart_Atlas.csv"


option_list <- list(
  make_option(c("-n", "--name"), type = "character", default = NULL, 
              help = "Experiment name", metavar = "name"),
  make_option(c("-d", "--directory"), type = "character", default = NULL, 
              help = "Path to the base directory", metavar = "directory"),
  make_option(c("-g", "--gtffile"), type = "character", default = NULL,
              help = "Path to the GTF file", metavar = "gtffile"),
  make_option(c("-m", "--meta_data"), type = "character", default = NULL,
              help = "Path to the meta data csv", metavar = "meta_data")
)

# Parse the options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the options are provided
if(is.null(opt$name)) {
  stop("No experiment name for bulk samples provided. Use -n or --name to specify.")
}
if(is.null(opt$directory)) {
  stop("No main directory for bulk samples provided. Use -d or --directory to specify the base directory.")
}
if(is.null(opt$gtffile)) {
  stop("No GTF file provided. Use -gtf or --gtffile to specify the GTF file path that has been used for the alignment.")
}
if(is.null(opt$meta_data)) {
  stop("No meta data of the samples provided. Use -meta or --meta_data to specify the meta_data file path that contains the sample name and patient information.")
}

# Assign the arguments to variables
base_directory <- opt$directory
gtf_file_path <- opt$gtffile
meta_data_csv <- opt$meta_data
name <- opt$name

perform_DESeq2_analysis <- function(base_directory, gtf_file_path, meta_data, BAM_extension_pattern = "_sorted.bam$") {
  # List all bam files
  file_list <- list.files(
    path = base_directory, 
    pattern = BAM_extension_pattern, 
    recursive = TRUE, 
    full.names = TRUE
  )
  
  # FeatureCounts
  fc <- featureCounts(
    files = file_list, 
    annot.ext = gtf_file_path,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    isPairedEnd = TRUE
  )
  
  # Process counts
  cts <- fc$counts
  colnames(cts) <- gsub(BAM_extension_pattern, "", colnames(cts))
  
  # Prepare metadata
  colData <- meta_data %>% dplyr::select(Sample, Replicate, Condition)
  rownames(colData) <- colData$Sample
  #colData$Sample <- NULL
  
  #check
  cts <- cts[, rownames(colData)]
  all(rownames(colData) == colnames(cts))
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = colData, design = ~ Condition)
  dds <- DESeq(dds)
  
  # Differential expression analysis
  #res <- results(dds)
  #res <- res[order(res$padj, na.last = NA), ]
  
  # Variance stabilizing transformation
  vsd <- vst(dds, blind=TRUE)
  
  return(list(cts = cts, vsd = vsd, dds = dds))
}


# Metadata
meta_data <- read.table(file.path(base_directory, "analysis", meta_data_csv), header = TRUE, sep = ";")

# Perform DESeq2 analysis
bulk_RNA_results <- perform_DESeq2_analysis(base_directory, gtf_file_path, meta_data)

cts <- bulk_RNA_results$cts
vsd <- bulk_RNA_results$vsd
dds <- bulk_RNA_results$dds


#saveRDS(cts, file.path(base_directory, "analysis", paste0(name, "_counts.rds")))

mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 98)
gene_annotations <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_accession", "entrezgene_id"),
                          filters = "ensembl_gene_id",
                          values = rownames(cts),
                          mart = mart)

# Sample Plotting
pdf(file.path(base_directory, "analysis", paste0(name, "_PCA.pdf")))

wb <- createWorkbook()

Conditions <- c("Hypoxia", "IL1b", "TGFb")

for(i in 1:length(Conditions)){
  res <- results(dds, contrast=c("Condition",Conditions[i],"Control"))
  res$ensembl_gene_id = rownames(res)
  res <- merge(data.frame(res), gene_annotations, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)
  res <- data.frame(res) %>% arrange(padj)
  
  #Plot DEG genes
  results = as.data.frame(dplyr::mutate(as.data.frame(res), sig=ifelse(res$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res))
  DEgenes_DESeq <- results[which(abs(results$log2FoldChange) > log2(1.5) & results$padj < 0.05),]
  
  p = ggplot2::ggplot(results, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
    ggplot2::geom_point(ggplot2::aes(col = sig)) +
    ggplot2::scale_color_manual(values = c("red", "black")) +
    ggplot2::ggtitle(paste0("Volcano Plot of DESeq2 analysis - ", Conditions[i]))
  
  DEG_plot <- p + ggrepel::geom_text_repel(data=results[1:10, ], ggplot2::aes(label=results$hgnc_symbol[1:10]))
  print(DEG_plot)
  
  #Save the DEG
  saveRDS(res, file.path(base_directory, "analysis", paste0(Conditions[i], "_res.rds")))
  addWorksheet(wb, Conditions[i])
  writeData(wb, sheet = Conditions[i], res)
}

saveWorkbook(wb, file = file.path(base_directory, "analysis", paste0(name, "_DESeq2.xlsx")), overwrite = TRUE)


#Heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Condition, vsd$Replicate, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
plot_heatmap <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
print(plot_heatmap)

# Plot PCA
pcaData <- plotPCA(vsd, intgroup=c("Condition", "Replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$Replicate <- as.character(pcaData$Replicate)
plot_PCA <- ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Replicate)) +
  geom_point() +
  theme_minimal() +
  theme(
    aspect.ratio = 1,
    plot.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
  )
print(plot_PCA)
dev.off()
