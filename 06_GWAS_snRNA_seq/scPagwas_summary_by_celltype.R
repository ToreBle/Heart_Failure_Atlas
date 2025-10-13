library(Seurat)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyverse)
library(broom)
library(effsize)
library(ggpubr)
library(scales)
library(reshape2)
library(ggdendro)
library(gridExtra)
library(grid)
library(ComplexHeatmap)
library(reshape2)
library(circlize)
library(ggsignif)


#-------------------------------------------------------------------------------
# get_summary_by_celltype()
#-------------------------------------------------------------------------------

#' @param patients       character vector of patient IDs
#' @param gwas_overview  data.frame with GWAS_Study etc.
#' @param clinical_meta  data.frame for AddMetaData()
#' @param data_dir       directory of Single‐cell RDS files
#' @param result_dir     root path to reactome_<Patient>_<GWAS> folders
#' @param celltypes      optional character vector; if NULL, uses all Annotation_1 levels
#'
#' @return data.frame of one row per (Patient × GWAS_Study × Annotation_1) with summary stats
get_summary_by_celltype <- function(patients,
                                    gwas_overview,
                                    clinical_meta,
                                    data_dir,
                                    result_dir,
                                    celltypes = NULL) {
  
  map_dfr(patients, function(patient) {
    single_data <- readRDS(file.path(data_dir, paste0(patient, ".rds")))
    
    map_dfr(unique(gwas_overview$GWAS_Study), function(gwas) {
      trs_csv <- file.path(
        result_dir,
        paste0("reactome_", patient, "_", gwas),
        paste0("reactome_", patient, "_", gwas, "_singlecell_scPagwas_score_pvalue.Result.csv")
      )
      
      if (!file.exists(trs_csv)) {
        warning(sprintf("Skipping missing file for %s / %s", patient, gwas))
        return(NULL)
      }
      
      # attach metadata
      meta_trs <- read.csv(trs_csv, row.names = 1)
      pagwas   <- AddMetaData(single_data, metadata = meta_trs)
      pagwas   <- AddMetaData(pagwas,   metadata = clinical_meta)
      md       <- pagwas@meta.data
      
      # optionally filter to a subset of celltypes
      if (!is.null(celltypes)) {
        md <- filter(md, Annotation_1 %in% celltypes)
      }
      
      # summarise per‐celltype
      md %>%
        group_by(Annotation_1) %>%
        summarise(
          mean_TRS     = mean(scPagwas.TRS.Score, na.rm = TRUE),
          sd_TRS       = sd(  scPagwas.TRS.Score, na.rm = TRUE),
          BMI          = paste(unique(BMI),        collapse = ","),
          Genotype     = paste(unique(Genotype),   collapse = ","),
          Sex          = paste(unique(Sex),        collapse = ","),
          Age          = paste(unique(Age),        collapse = ","),
          Hypertension = paste(unique(Hypertension),collapse = ","),
          nCount_RNA   = mean(nCount_RNA,            na.rm = TRUE),
          nFeature_RNA = mean(nFeature_RNA,          na.rm = TRUE),
          Study        = paste(unique(Study),       collapse = ","),
          Diabetes     = paste(unique(Diabetes),    collapse = ",")
        ) %>%
        ungroup() %>%
        mutate(
          Patient     = patient,
          Condition   = paste(strsplit(patient, "_")[[1]][1:2], collapse = "_"),
          Cluster     = strsplit(patient, "_")[[1]][1],
          GWAS_Study  = gwas
        ) %>%
        left_join(gwas_overview, by = "GWAS_Study")
    })
  })
}


#-------------------------------------------------------------------------------
# Prepare inputs
#-------------------------------------------------------------------------------

# Define your inputs
GWAS_overview <- read.table("./heart-atlas/scPagwas/GWAS_studies_summary.txt", sep="\t")
colnames(GWAS_overview) <- c("GWAS_Study", "Name", "Category")
GWAS_overview <- GWAS_overview %>%
  filter(Name %in% c("Coronary atherosclerosis", "Hypertension", "Type 2 diabetes", "Type 1 diabetes", "Obesity", "Myocarditis", "Myocardial infarction", "Arrythmia", "Dilated Cardiomyopathy 1", "Dilated Cardiomyopathy 2", "Hyperthropic Cardiomyopathy", "Ischemic stroke"))



clinical_meta_data <- read.csv("./heart-atlas/scPagwas/meta_data_all_Heart_atlas.csv", row.names = 2)


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


data_dir <- "./heart-atlas/scPagwas/data/"
result_dir <- "./heart-atlas/scPagwas/outs/"

#-------------------------------------------------------------------------------
# Run analysis aggregated by cell types
#-------------------------------------------------------------------------------

summary_cell_types <- get_summary_by_celltype(
  patients      = Patients,
  gwas_overview = GWAS_overview,
  clinical_meta = clinical_meta_data,
  data_dir      = data_dir,
  result_dir    = result_dir
)

print(head(summary_cell_types))
summary_cell_types$Annotation_1[summary_cell_types$Annotation_1 == "Monocytes"] <- "Myeloid"
summary_cell_types$Condition <- gsub("^[01234]_", "", summary_cell_types$Condition)
summary_cell_types$Condition[summary_cell_types$Condition == "HCM_Ch"] <- "HCM"

#Heatmap

pdf("./heart-atlas/scPagwas/figures/Heatmap_subset.pdf", width = 12, height = 6)
All <- summary_cell_types %>% filter(Annotation_1 == "Cardiomyocytes")
All <- All %>% filter(Name != "Hyperthropic Cardiomyopathy")
All <- All %>% filter(Condition != "HCM")

data_matrix <- reshape2::acast(All, Name ~ Patient, value.var = "mean_TRS")
print(head(data_matrix))
col_annotation_values <- All[match(colnames(data_matrix), All$Patient), c("Condition", "BMI", "Sex"), drop = FALSE]
col_annotation_values$BMI <- as.numeric(as.character(col_annotation_values$BMI))


# Column annotations only (using exact values from `df`)
col_annotation <- HeatmapAnnotation(
  Condition = col_annotation_values$Condition,
  BMI = col_annotation_values$BMI,
  Sex = col_annotation_values$Sex,
  #Hypertension = col_annotation_values$Hypertension,
  col = list(
    Condition = c("DCM" = "orange", "Control"="#2e86c1"),  # Define colors for conditions
    Sex = c("Female" = "red", "Male"='#8e44ad'),
    #Hypertension = c("Hypertension"='red', "normal" = "blue", "PAH" = "orange"),
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
  #show_column_names = FALSE,
  col = colorRampPalette(c("blue", "white", "red"))(50),
  top_annotation = col_annotation,
  show_heatmap_legend = TRUE
)
dev.off()


pdf("./heart-atlas/scPagwas/figures/Correlation_Cardiomyocytes.pdf")
All <- summary_cell_types %>% filter(Annotation_1 == "Cardiomyocytes")
All$BMI <- as.numeric(as.character(All$BMI))

print(ggplot(All %>% filter(GWAS_Study =="GCST90018834-EFO_0000407"), aes(x=mean_TRS, y=Condition))+
        geom_boxplot(outlier.size = 0)+geom_jitter(aes(color=Condition, size=4))+
        scale_color_manual(values = c("gray70", "orange", "#D11D36","#940031"))+
        theme_classic()+theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))+
        ggtitle("DCM 1"))
print(ggplot(All %>% filter(GWAS_Study =="GCST90018834-EFO_0000407"), aes(x=mean_TRS, y=Condition))+
        geom_boxplot(outlier.size = 0)+geom_jitter(aes(color=Genotype, size=4))+
        theme_classic()+theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))+
        ggtitle("DCM 1"))
print(ggplot(All %>% filter(GWAS_Study =="Jurgens_DCM"), aes(x=mean_TRS, y=Condition))+
        geom_boxplot(outlier.size = 0)+geom_jitter(aes(color=Condition, size=4))+
        scale_color_manual(values = c("gray70", "orange", "#D11D36","#940031"))+
        theme_classic()+theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))+
        ggtitle("DCM 2"))
print(ggplot(All %>% filter(GWAS_Study =="GCST90435254"), aes(x=mean_TRS, y=Condition))+
        geom_boxplot(outlier.size = 0)+geom_jitter(aes(color=Condition, size=4))+
        scale_color_manual(values = c("gray70", "orange", "#D11D36","#940031"))+
        theme_classic()+theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))+
        ggtitle("HCM"))
print(ggplot(All %>% filter(GWAS_Study =="GCST90043957-MONDO_0021661"), aes(x=mean_TRS, y=Condition))+
        geom_boxplot(outlier.size = 0)+geom_jitter(aes(color=Condition, size=4))+
        scale_color_manual(values = c("gray70", "orange", "#D11D36"))+
        theme_classic()+theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))+
        ggtitle("Coronary Atherosclerosis"))
print(ggplot(All %>% filter(GWAS_Study =="GCST90043949-EFO_0000537"), aes(x=mean_TRS, y=Condition))+
        geom_boxplot(outlier.size = 0)+geom_jitter(aes(color=Hypertension, size=3))+
        guides(colour = guide_legend(override.aes = list(size=6)))+
        scale_color_manual(values = c("gray70", "orange", "#D11D36","blue"))+
        theme_classic()+theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24), legend.title = element_text(size = 21, color="black"), legend.text  = element_text(size = 18, color="black"),)+
        ggtitle("Hypertension"))

print(ggplot(All %>% filter(GWAS_Study == "GCST90043672-EFO_0001073"), 
             aes(x = mean_TRS, y = BMI, color = Condition)) +
        geom_point(size = 4) +
        scale_color_manual(values = c("gray70", "orange", "#D11D36"))+
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

genotyped <- All %>% filter(Name == "Dilated Cardiomyopathy 1" | Name == "Dilated Cardiomyopathy 2") %>%
  filter(Condition != "HCM") %>%
  filter(Patient != "0_Control_Ch_P1582") # non HF donor with NEBL mutation that would be single point in ctrl
  ggplot(genotyped, aes(x = mean_TRS, y = Condition, color = Genotype)) +
  geom_boxplot(outlier.size = 0, position = position_dodge(width = 0.7))+geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7))+
  scale_color_manual(values = c("black", "#005184", "#0070AD", "#0087BB", "#3C9FC8", "#7CB5D5", "#A7CAE0", "orange", "#E9A79D", "gray70", "#E18A7B", "#D76A53", "#C94816", "#AE2D00", "#872000", "#611300"))+
  facet_grid(Name ~ ., scales = "free", space = "free") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(size = 16, color="black"),
    axis.text.y  = element_text(size = 16, color="black"),
    legend.title = element_text(size = 14, color="black"),
    legend.text  = element_text(size = 13, color="black"),
    strip.text   = element_text(size = 14, color="black")
  )
  

  HCM_DCM <- All %>% filter(Name == "Hyperthropic Cardiomyopathy" | Name == "Dilated Cardiomyopathy 2")
  ggplot(HCM_DCM, aes(x = mean_TRS, y = Condition, color = Condition)) +
    geom_boxplot(outlier.size = 0, position = position_dodge(width = 0.7))+geom_point(aes(size=3), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7))+
    scale_color_manual(values = c("gray70", "orange", "#D11D36"))+
    facet_grid(Name ~ ., scales = "free", space = "free") +
    theme_classic() +
    theme(
      axis.text.x  = element_text(size = 21, color="black"),
      axis.title  = element_text(size = 21, color="black"),
      axis.text.y  = element_text(size = 24, color="black"),
      strip.text   = element_text(size = 16, color="black"),
      legend.position = ""
    )+
    labs(y="")

dev.off()