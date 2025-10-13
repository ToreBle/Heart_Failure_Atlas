library(dplyr)
library(Seurat)
library(scPagwas)
library(tidyverse)
library(ggplot2)
library(grDevices)
library(stats)
library(FactoMineR)
library(scales)
library(gridExtra)
library(grid)
library(ComplexHeatmap)
library(reshape2)
library(circlize)
library(ggsignif)
library(ggpubr)

source(system.file("extdata", "plot_scpathway_contri_dot.R", package = "scPagwas"))

# Set in parallel
path = "./heart-atlas/scPagwas/outs/"
fig_out = "./heart-atlas/scPagwas/figures/"

#Read only studies
GWAS_overview <- read.table("./heart-atlas/scPagwas/GWAS_studies_summary.txt", sep="\t")
colnames(GWAS_overview) <- c("GWAS_Study", "Name", "Category")

GWAS_overview <- GWAS_overview %>%
  filter(Name %in% c("Coronary atherosclerosis", "Hypertension", "Type 2 diabetes", "Type 1 diabetes", "Obesity", "Myocarditis", "Myocardial infarction", "Arrythmia", "Dilated Cardiomyopathy 1", "Dilated Cardiomyopathy 2", "Hyperthropic Cardiomyopathy", "Ischemic stroke", "Left ventricle ejection fraction"))

clinical_meta_data <- read.csv("./heart-atlas/scPagwas/meta_data_all_Heart_atlas.csv", row.names = 2)
clinical_meta_data$X.1 <- NULL
clinical_meta_data$X_scvi_batch <- NULL
clinical_meta_data$X_scvi_labels <- NULL
clinical_meta_data$leiden02 <- NULL
clinical_meta_data$leiden03 <- NULL
clinical_meta_data$leiden04 <- NULL

pdf(paste0(fig_out, "TRS_Patient_GWAS.pdf"))
Patients <- c("3_DCM_Re_DS3",  "3_DCM_Re_H27", "3_DCM_Re_H01", "3_DCM_Re_H28","3_DCM_Re_H03", "3_DCM_Re_H33",
              "3_DCM_Re_H06", "3_DCM_Re_H35", "3_DCM_Re_H13", "3_DCM_Re_H36","3_DCM_Re_H18", "3_DCM_Re_H37",
              "3_DCM_Re_H19", "3_DCM_Re_H40", "3_DCM_Re_H20", "3_DCM_Re_H41", "3_DCM_Re_H22", "3_DCM_Re_H56", "3_DCM_Re_H24","3_DCM_Re_H78",
              
              "1_DCM_Kö_H_ZC-LVAD",   "1_DCM_Kö_TWCM-11-3", "1_DCM_Kö_TWCM-13-102", "1_DCM_Kö_TWCM-13-17",  
              "1_DCM_Kö_TWCM-13-208", "1_DCM_Kö_TWCM-13-285", "1_DCM_Kö_TWCM-13-47",  "1_DCM_Kö_TWCM-13-84",  "1_DCM_Kö_TWCM-LVAD2",  
              "1_DCM_Kö_TWCM-LVAD3",  "1_DCM_Re_DL3", "1_DCM_Re_DP1", "1_DCM_Re_DP2", "1_DCM_Re_H02", "1_DCM_Re_H15",         
              "1_DCM_Re_H21", "1_DCM_Re_H32", "1_DCM_Re_H57", "1_DCM_Re_H58", "1_DCM_Re_IC_H03", "1_DCM_Re_IC_H04",
              
              "2_DCM_Ch_P1290", "2_DCM_Ch_P1430", "2_DCM_Ch_P1617", "2_DCM_Ch_P1300", "2_DCM_Ch_P1437", "2_DCM_Re_DT4",
              "2_DCM_Ch_P1304", "2_DCM_Ch_P1472", "2_DCM_Re_H04", "2_DCM_Ch_P1358", "2_DCM_Ch_P1504", "2_DCM_Re_H16", "2_DCM_Ch_P1371", "2_DCM_Ch_P1606",
              
              "0_DCM_Re_H12", "0_DCM_Re_H14", "0_DCM_Re_H25", "0_DCM_Re_H26", "0_DCM_Re_H31", "0_DCM_Re_H34", 
              "0_DCM_Re_H39", "4_DCM_Re_DL2", "4_DCM_Re_H17", "4_DCM_Re_IC_H01",
              "HCM_Ch_P1422", "HCM_Ch_P1425", "HCM_Ch_P1447", "HCM_Ch_P1462", "HCM_Ch_P1479", "HCM_Ch_P1508",
              "HCM_Ch_P1510", "HCM_Ch_P1602", "HCM_Ch_P1630", "HCM_Ch_P1631", "HCM_Ch_P1685", "HCM_Ch_P1707", "HCM_Ch_P1722", "HCM_Ch_P1726", "HCM_Ch_P1735",
              
              "0_Control_Ch_P1515", "0_Control_Ch_P1516", "0_Control_Ch_P1539", "0_Control_Ch_P1540", "0_Control_Ch_P1547", 
              "0_Control_Ch_P1558", "0_Control_Ch_P1582", "0_Control_Ch_P1600", "0_Control_Ch_P1603", "0_Control_Ch_P1610", "0_Control_Ch_P1678",
              "0_Control_Re_H2", "0_Control_Re_H3", "0_Control_Re_H4", "0_Control_Re_H49", "0_Control_Re_H5",
              "0_Control_Re_H53", "0_Control_Re_H6", "0_Control_Re_H7")
CM <- data.frame()
MYL <- data.frame()
Cells <- data.frame()
for(Patient in Patients){
  Patient_GWAS <- data.frame()
  Single_data <- readRDS(paste0("./heart-atlas/scPagwas/data/", Patient, ".rds"))
  
  for (GWAS in unique(GWAS_overview$GWAS_Study)) {
    print(paste("Processing GWAS study:", GWAS))
    
    # Load meta_data_TRS file for the current GWAS study
    trs_file <- paste0(path, "reactome_", Patient, "_", GWAS, "/reactome_", Patient, "_", GWAS, "_singlecell_scPagwas_score_pvalue.Result.csv")
    
    if (file.exists(trs_file)) {
      meta_data_TRS <- read.csv(trs_file, header = TRUE, row.names = 1)
      Pagwas <- AddMetaData(Single_data, metadata = meta_data_TRS)
      Pagwas <- AddMetaData(Pagwas, metadata = clinical_meta_data)
      
      # Calculate summary statistics for scPagwas.TRS.Score
      St <- Pagwas@meta.data %>%
        group_by(Annotation_1) %>%
        summarise(
          mean_TRS = mean(scPagwas.TRS.Score, na.rm = TRUE),
          sd_TRS = sd(scPagwas.TRS.Score, na.rm = TRUE),
          BMI = paste(unique(BMI), collapse = ","),
          Genotype = paste(unique(Genotype), collapse = ","),
          Sex = paste(unique(Sex), collapse = ","),
          Age = paste(unique(Age), collapse = ","),
          Hypertension = paste(unique(Hypertension), collapse = ","),
          nCount_RNA = mean(nCount_RNA, na.rm = TRUE),
          nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
          Study = paste(unique(Study), collapse = ","),
          Diabetes = paste(unique(Diabetes), collapse = ","),
          LVEF = paste(unique(LVEF), collapse = ","),
        )
      # Add additional columns
      St$Patient <- Patient
      St$Condition <- paste(strsplit(Patient, "_")[[1]][1:2], collapse = "_") #Condition
      St$Cluster <-  paste(strsplit(Patient, "_")[[1]][1], collapse = "_")
      St$GWAS_Study <- GWAS
      
      # Join with GWAS_overview for additional metadata
      St <- St %>% left_join(GWAS_overview, by = "GWAS_Study")
      
      # Combine results into the final dataframe
      Patient_GWAS <- bind_rows(Patient_GWAS, St)
      
    } else {
      print(paste("GWAS", GWAS, "for Patient", Patient, "does not exist."))
    }
    Cell <- Pagwas@meta.data %>% select(Donor, Subclustering, Condition2, Cluster_Condition, scPagwas.TRS.Score, scPagwas.downTRS.Score, scPagwas.gPAS.score, Annotation_1)
    Cell$GWAS_Study <- GWAS
    Cell <- Cell %>% left_join(GWAS_overview, by = "GWAS_Study")
    Cells <- bind_rows(Cell, Cells)
  }
  CM_Pat <- Patient_GWAS %>% filter(Annotation_1=="Cardiomyocytes")
  CM <- bind_rows(CM_Pat, CM)
  MYL_Pat <- Patient_GWAS %>% filter(Annotation_1=="Monocytes")
  MYL <- bind_rows(MYL_Pat, MYL)
}
dev.off()


Cells$Annotation_1[Cells$Annotation_1 == "Monocytes"] <- "Myeloid"
Cells$Condition <- gsub("^[01234]_", "", Cells$Condition)
saveRDS(Cells, "./heart-atlas/scPagwas/Cells_overview.rds")


Statistics <- Cells %>%
  group_by(Annotation_1,Name,Category) %>%
  summarise(
    mean_TRS = mean(scPagwas.TRS.Score, na.rm = TRUE),
    sd_TRS = sd(scPagwas.TRS.Score, na.rm = TRUE),
    mean_gPAS = mean(scPagwas.TRS.Score, na.rm = TRUE),
    sd_gPAS = sd(scPagwas.TRS.Score, na.rm = TRUE),
    mean_TRS_down = mean(scPagwas.downTRS.Score, na.rm = TRUE),
    sd_TRS_down = sd(scPagwas.downTRS.Score, na.rm = TRUE))

pdf(paste0(fig_out, "gPAS_TRS_cell_types.pdf"), width = 8, height = 7)
summary_TRS_plot <- ggplot(Statistics, aes(x=Annotation_1, y=Name, color=mean_TRS, size =sd_TRS))+
  geom_point() +
  scale_color_gradientn(colors = rev(hcl.colors(20, "Reds")))+
  #scale_fill_gradientn(colors = hcl.colors(20, "Reds"))+
  theme_classic()+
  theme(
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    strip.text = element_text(size = 13))+
  facet_grid(Category ~ ., scales = "free", space = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14))
print(summary_TRS_plot)
summary_TRS_plot_down <- ggplot(Statistics, aes(x=Annotation_1, y=Name, color=mean_TRS_down, size =sd_TRS_down))+
  geom_point() +
  scale_color_gradientn(colors = rev(hcl.colors(20, "Blues")))+
  #scale_fill_gradientn(colors = hcl.colors(20, "Reds"))+
  theme_classic()+
  theme(
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    strip.text = element_text(size = 13))+
  facet_grid(Category ~ ., scales = "free", space = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14))
print(summary_TRS_plot_down)
summary_gPAS_plot <- ggplot(Statistics, aes(x=Annotation_1, y=Name, color=mean_gPAS, size =sd_gPAS))+
  geom_point() +
  scale_color_gradientn(colors = rev(hcl.colors(20, "Reds")))+
  #scale_fill_gradientn(colors = hcl.colors(20, "Reds"))+
  theme_classic()+
  theme(
    axis.text.y = element_text(size = 14, colour = "black"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    strip.text = element_text(size = 13))+
  facet_grid(Category ~ ., scales = "free", space = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14, colour = "black"))
print(summary_gPAS_plot)
dev.off()  

saveRDS(CM, "./heart-atlas/scPagwas/All_CM_overview.rds")


pdf(paste0(fig_out, "Heatmap_subset.pdf"), width = 12, height = 6)
All <- CM
All$Condition <- gsub("^[01234]_", "", All$Condition)
All <- All %>% filter(Name != "Hyperthropic Cardiomyopathy")
All <- All %>% filter(Name != "Left ventricle ejection fraction")
All <- All %>% filter(Condition != "HCM_Ch")

data_matrix <- reshape2::acast(All, Name ~ Patient, value.var = "mean_TRS")
print(head(data_matrix))
col_annotation_values <- All[match(colnames(data_matrix), All$Patient), c("Condition", "BMI", "Sex"), drop = FALSE]
col_annotation_values$BMI <- as.numeric(as.character(col_annotation_values$BMI))

# Column annotations only (using exact values from `df`)
col_annotation <- HeatmapAnnotation(
  Condition = col_annotation_values$Condition,
  BMI = col_annotation_values$BMI,
  Sex = col_annotation_values$Sex,
  col = list(
    Condition = c("DCM" = "orange", "Control"="#2e86c1"),  # Define colors for conditions
    Sex = c("Female" = "red", "Male"='#8e44ad'),
    BMI = colorRamp2(c(min(col_annotation_values$BMI, na.rm = TRUE), mean(col_annotation_values$BMI, na.rm = TRUE), max(col_annotation_values$BMI, na.rm = TRUE)), 
                     c("darkred", "orange", "white"))
  )
)

# Create the heatmap
Heatmap(
  data_matrix,
  name = "mean TRS for\nCardiomyocytes",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  show_column_names = FALSE,
  row_names_max_width = unit(12, "cm"),
  row_names_side = "left",
  column_names_side = "top",
  col = colorRampPalette(c("blue", "white", "red"))(50),
  top_annotation = col_annotation,
  show_heatmap_legend = TRUE
)
dev.off()


pdf(paste0(fig_out, "Correlation_Cardiomyocytes.pdf"))
All <- CM
All$Condition <- gsub("^[01234]_", "", All$Condition)

All$LVEF <- as.numeric(as.character(All$LVEF))
All$BMI <- as.numeric(as.character(All$BMI))


print(ggplot(All %>% filter(GWAS_Study =="GCST90018834-EFO_0000407"), aes(x=mean_TRS, y=Condition))+
        geom_boxplot()+geom_jitter(aes(color=Condition))+
        theme_classic()+theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14))+
        ggtitle("DCM 1"))
print(ggplot(All %>% filter(GWAS_Study =="Jurgens_DCM"), aes(x=mean_TRS, y=Condition))+
        geom_boxplot(outlier.size = 0)+geom_jitter(aes(color=Condition, size=4))+
        theme_classic()+theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))+
        ggtitle("DCM 2"))
print(ggplot(All %>% filter(GWAS_Study =="GCST90435254"), aes(x=mean_TRS, y=Condition))+
        geom_boxplot(outlier.size = 0)+geom_jitter(aes(color=Condition, size=4))+
        theme_classic()+theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))+
        ggtitle("HCM"))
print(ggplot(All %>% filter(GWAS_Study =="GCST90043957-MONDO_0021661"), aes(x=mean_TRS, y=Condition))+
        geom_boxplot(outlier.size = 0)+geom_jitter(aes(color=Condition, size=4))+
        theme_classic()+theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))+
        ggtitle("Coronary Atherosclerosis"))
print(ggplot(All %>% filter(GWAS_Study =="GCST90043949-EFO_0000537"), aes(x=mean_TRS, y=Condition))+
        geom_boxplot(outlier.size = 0)+geom_jitter(aes(color=Hypertension, size=4))+
        theme_classic()+theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))+
        ggtitle("Hypertension"))
print(ggplot(All %>% filter(GWAS_Study =="GCST90043672-EFO_0001073"), aes(x=mean_TRS, y=BMI, color=Condition))+
        geom_point(size=4)+geom_smooth(aes(color = NULL), method = "lm", se = FALSE)+
        theme_classic()+theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))+
        ggtitle("Obesity"))

print(ggplot(All %>% filter(GWAS_Study == "GCST90043672-EFO_0001073"), 
             aes(x = mean_TRS, y = BMI, color = Condition)) +
        geom_point(size = 4) +
        geom_smooth(aes(color = NULL), method = "lm", se = FALSE) +
        stat_cor(
          aes(x = mean_TRS, y = BMI),
          method = "pearson",
          label.x.npc = "left",
          label.y.npc = "top",
          size = 6,
          inherit.aes = FALSE
        ) +
        theme_classic() +
        theme(
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 24)
        ) +
        ggtitle("Obesity - pearson")
)

print(ggplot(All %>% filter(GWAS_Study == "GCST90043672-EFO_0001073"), 
             aes(x = mean_TRS, y = BMI, color = Condition)) +
        geom_point(size = 4) +
        geom_smooth(aes(color = NULL), method = "lm", se = FALSE) +
        stat_cor(
          aes(x = mean_TRS, y = BMI),
          method = "spearman",
          label.x.npc = "left",
          label.y.npc = "top",
          size = 6,
          inherit.aes = FALSE
        ) +
        theme_classic() +
        theme(
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 24)
        ) +
        ggtitle("Obesity")
)


print(ggplot(All %>% filter(GWAS_Study == "GCST90268128"), 
             aes(x = mean_TRS, y = LVEF, color = Condition)) +
        geom_point(size = 4) +
        geom_smooth(aes(color = NULL), method = "lm", se = FALSE) +
        stat_cor(
          aes(x = mean_TRS, y = LVEF),
          method = "pearson",
          label.x.npc = "left",
          label.y.npc = "top",
          size = 6,
          inherit.aes = FALSE
        ) +
        theme_classic() +
        theme(
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 24)
        ) +
        ggtitle("LVEF - pearson")
)

print(ggplot(All %>% filter(GWAS_Study == "GCST90268128"), 
             aes(x = mean_TRS, y = LVEF, color = Condition)) +
        geom_point(size = 4) +
        geom_smooth(aes(color = NULL), method = "lm", se = FALSE) +
        stat_cor(
          aes(x = mean_TRS, y = LVEF),
          method = "spearman",
          label.x.npc = "left",
          label.y.npc = "top",
          size = 6,
          inherit.aes = FALSE
        ) +
        theme_classic() +
        theme(
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 24)
        ) +
        ggtitle("LVEF")
)

dev.off()