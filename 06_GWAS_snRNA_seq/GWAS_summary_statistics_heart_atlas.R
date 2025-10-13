library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(openxlsx)

setwd("./scPAGWAS/GWAS_stats/")

# All ancestry data v1.0 from https://www.ebi.ac.uk/gwas/docs/file-downloads
# Direct download: https://www.ebi.ac.uk/gwas/api/search/downloads/ancestry
GWAS_stats <- read.csv("gwas_catalog-ancestry_2025_heart_atlas.tsv", sep="\t")

GWAS_overview <- read.table("GWAS_studies_summary.txt", sep="\t")
colnames(GWAS_overview) <- c("STUDY.ACCESSION", "Name", "Category")
GWAS_overview$STUDY.ACCESSION <- sub("-.*", "", GWAS_overview$STUDY.ACCESSION)
GWAS_Heart_atlas <- left_join(GWAS_overview, GWAS_stats, by="STUDY.ACCESSION")

palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(GWAS_Heart_atlas$BROAD.ANCESTRAL.CATEGORY))+1)

pdf("GWAS_summary_statistics.pdf", width = 8, height = 8)
ggplot(GWAS_Heart_atlas, aes(y=Name, x=NUMBER.OF.INDIVDUALS/1000000, fill=BROAD.ANCESTRAL.CATEGORY))+geom_bar(position="stack", stat="identity")+
  facet_grid(Category ~., space = "free", scales = "free")+
  scale_fill_manual(values = palette, labels = str_replace_all(levels(as.factor(GWAS_Heart_atlas$BROAD.ANCESTRAL.CATEGORY)), ", ", "\n"))+
  theme_bw()+labs(y="", x="Number of individuals [in million]", fill="Ancestral Category")+
  scale_x_continuous(expand = c(0,0))+
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text = element_text(size=10, colour = "black"),
    legend.key.height = unit(0.8, "cm"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
dev.off()

#save as SI table in excel
wb <- createWorkbook()
addWorksheet(wb, "GWAS studies")
writeData(wb, sheet = "GWAS studies", GWAS_Heart_atlas)
saveWorkbook(wb, file = "Supplementary Table 10 - Overview GWAS summary statistics.xlsx", overwrite = TRUE)

