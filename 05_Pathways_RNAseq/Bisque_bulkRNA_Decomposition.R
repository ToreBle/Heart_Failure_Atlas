library(Biobase)
library(BisqueRNA)
library(Seurat)
library(biomaRt)
library(tidyverse)
library(scales)
library(ggsignif)
library(ggpubr)

setwd("./heart-atlas/bulkRNA/analysis/")

counts <- readRDS("Bulk_RNA_Heart_Atlas_counts.rds")
counts <- data.frame(counts)
counts$ensembl_gene_id <- rownames(counts)

mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 99) #98
gene_annotations <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                          filters = "ensembl_gene_id",
                          values = counts$ensembl_gene_id,
                          mart = mart)
bulk <- merge(counts, gene_annotations, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)

bulk <- bulk %>%
  add_count(external_gene_name) %>%
  filter(n == 1) %>%
  select(-n)

rownames(bulk) <- bulk$external_gene_name
bulk$ensembl_gene_id <- NULL
bulk$external_gene_name <- NULL

bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(bulk))


Seurat <- readRDS("./heart-atlas/multinichenet/Heart_atlas_reduced.rds")

Seurat <- subset(
  Seurat,
  subset = Annotation_2 %in% c('Fibroblast')
)
Single_data <- subset(Seurat, subset = Subclustering_2 != "Fibroblast, Perivascular")


#print(table(Single_data$Subclustering))
print(table(Single_data$Subclustering_2))
sample.ids <- colnames(Single_data)
sc_counts <- as.matrix(GetAssayData(Single_data, layer = "counts"))


# individual.ids and cell.types should be in the same order as in sample.ids
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=sample.ids,
                       SubjectName=Single_data$Donor,
                       cellType=Single_data$Subclustering_2)
sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.pdata <- new("AnnotatedDataFrame",
                data=sc.pheno,
                varMetadata=sc.meta)
sc.eset <- Biobase::ExpressionSet(assayData=sc_counts,
                                  phenoData=sc.pdata)

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=F)
write.csv(res$bulk.props, file = "Proportions_BisqueRNA.csv")

proportions <- data.frame(res$bulk.props)
proportions$Cell_type <- rownames(proportions)
df_long <- proportions %>%
  pivot_longer(
    cols = starts_with("ds"),
    names_to = "Sample",
    values_to = "Proportion"
  )

Samples <- read.csv("BulkRNA_samples_Heart_Atlas.csv", sep=";")
Samples$Replicate <- NULL
combined <- merge(df_long, Samples, by = "Sample", all.x = TRUE)
#combined$Cell_types <- gsub(" \\(", "\n(", combined$Cell_type)
combined$Cell_types <- gsub(" ", "\n", combined$Cell_type)


sample_colors <- c("darkgray", "#2e86c1", "#ff7f0e", "#d62728")

pdf("Proportions_Fibroblast.pdf")#, width = 7, height = 6)
ggplot(combined, aes(y = Condition, x = Proportion, fill=Condition)) +
  geom_boxplot() +
  theme_classic()+
  facet_wrap(Cell_types ~ ., scales = "free", ncol = 2)+theme(strip.text.y.left = element_text(angle = 0, size = 22))+
  theme(line = element_blank(), axis.title = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 16), legend.text = element_text(colour = "black", size = 16), axis.text.x = element_text(colour = "black", size = 16), axis.text.y =  element_blank()) +
  labs(x="Proportion", y="", fill="Treatment")+
  scale_fill_manual(values = sample_colors)+
  geom_signif(
    comparisons = list(c("Hypoxia", "Control"), c("IL1b", "Control"), c("TGFb", "Control")), test = "t.test",
    map_signif_level = TRUE, textsize = 4, step_increase = 0.2, margin_top = 0.1)+
  scale_x_continuous(labels = scales::percent, breaks=c(0, 0.2, 0.4, 0.6))


ggplot(combined, aes(y = Cell_types, x = Proportion, fill=Cell_types)) +
  geom_boxplot() +
  theme_classic()+
  facet_wrap(Condition ~ ., scales = "free", ncol = 2)+theme(strip.text.y.left = element_text(angle = 0, size = 22))+
  theme(line = element_blank(), axis.title = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 16), legend.text = element_text(colour = "black", size = 16), axis.text.x = element_text(colour = "black", size = 16), axis.text.y =  element_blank()) +
  labs(x="Proportion", y="", fill="Treatment")+
  #scale_fill_manual(values = sample_colors)+
  geom_signif(
    comparisons = list(c("Hypoxia", "Control"), c("IL1b", "Control"), c("TGFb", "Control")), test = "t.test",
    map_signif_level = TRUE, textsize = 4, step_increase = 0.2, margin_top = 0.1)+
  scale_x_continuous(labels = scales::percent, breaks=c(0, 0.2, 0.4, 0.6))


ggplot(combined, aes(y = Condition, x = Proportion, fill=Condition)) +
  geom_boxplot() +
  theme_minimal()+
  facet_grid(Cell_types ~ ., scales = "free")+theme(strip.text = element_text(angle = 0, size = 12))+
  theme(line = element_blank(), legend.position = "right", axis.title = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 16), legend.text = element_text(colour = "black", size = 16), axis.text.x = element_text(colour = "black", size = 16), axis.text.y =  element_text(colour = "black", size = 16)) +
  labs(x="Proportion", y="", fill="Treatment")+
  scale_fill_manual(values = sample_colors)+
  geom_signif(
    comparisons = list(c("Hypoxia", "Control"), c("IL1b", "Control"), c("TGFb", "Control")), test = "wilcox.test",
    map_signif_level = TRUE, textsize = 4, step_increase = 0.1, margin_top = 0.01)+
  scale_x_continuous(labels = scales::percent, breaks=c(0, 0.2, 0.4, 0.6))


ggplot(combined, aes(y = Proportion, x = Condition, fill = Condition)) +
  geom_boxplot() +
  theme_minimal()+
  facet_grid(. ~ Cell_types, scales = "free")+theme(strip.text = element_text(angle = 90, size = 12))+
  theme(line = element_blank(), legend.position = "bottom", axis.title = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 16), legend.text = element_text(colour = "black", size = 16), axis.text.x = element_text(colour = "black", size = 14), axis.text.y =  element_text(colour = "black", size = 16)) +
  labs(x="", y="Proportion", fill="Treatment")+
  scale_fill_manual(values = sample_colors)+
  RotatedAxis()+
  geom_signif(
    comparisons = list(c("Hypoxia", "Control"), c("IL1b", "Control"), c("TGFb", "Control")), test = "wilcox.test",
    map_signif_level = TRUE, textsize = 4, step_increase = 0.1, margin_top = 0.01)+
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6))

dev.off()