#Load libraries
library(Seurat)
library(tidyverse)
library(dplyr)
require(biomaRt)
require(SingleCellExperiment)
require(liana)
require(purrr)
require(reticulate)
library(magrittr)
library(SeuratObject)
require(tibble)
library(patchwork)

# creating gpcr and consensus external resorce 
lr_df = read.csv('41586_2023_6311_MOESM5_ESM.xlsx - interaction_input.csv')
lr_clean <- lr_df %>%
  filter(!is.na(protein_name_a) & !is.na(protein_name_b) &
           protein_name_a != "" & protein_name_b != "")
lr_custom <- lr_clean %>%
  transmute(
    source_genesymbol = protein_name_a,    # ligand
    target_genesymbol = protein_name_b,    # receptor
    source  = partner_a,                   
    target  = partner_b,
    sources = source,                     
    references = NA_character_,            
    category_intercell_source = "ligand",  
    category_intercell_target = "receptor"
  )
write.csv(lr_custom, file = 'custom_gpcr.csv')

cpdb_res <- liana::select_resource("Consensus")[[1]]

combined_res <- bind_rows(cpdb_res, lr_custom) %>%
  distinct(source_genesymbol, target_genesymbol, .keep_all = TRUE) 

#write.csv(cpdb_res, file = 'cpdb_res.csv')
#write.csv(combined_res, file = 'combined_res_cons_gpcr.csv')

#Read file
SeuratObj = readRDS('path/to/seurat/object')
Idents(SeuratObj) = SeuratObj$Annotation_1
SeuratObj = subset(SeuratObj, idents ='Ischemic cells (MI)', invert = TRUE )

# For Post MI
Idents(SeuratObj) = SeuratObj$Study
Post_MI <- subset(SeuratObj, subset = Condition %in% c('Post MI','Control'))

Idents(SeuratObj) = SeuratObj$Annotation_1
SeuratObj = subset(SeuratObj, idents =c('Mast cells', 'Adipocytes', 'Lymphatic EC' ), invert = TRUE )

Idents(Post_MI) = Post_MI@meta.data$Annotation_1

data = as.SingleCellExperiment(Post_MI)
gp <-assays(data)[["counts"]][,colData(data)$Condition=='Control']
pp <-assays(data)[["counts"]][,colData(data)$Condition=='Post MI']

metadata_gp <-tibble::tibble(rownames(colData(data))[colData(data)$Condition=='Control'],
                             as.character(colData(data)$Annotation_1[colData(data)$Condition=='Control']))
metadata_pp<- tibble::tibble(rownames(colData(data))[colData(data)$Condition=='Post MI'],
                             as.character(colData(data)$Annotation_1[colData(data)$Condition=='Post MI']))
#create seurat object for running liana for conditions separately
good = CreateSeuratObject(counts = gp)
Idents(good) = metadata_gp$"as.character(...)"
poor = CreateSeuratObject(counts = pp)
Idents(poor) = metadata_pp$"as.character(...)"
combined_res = read.csv('read_gpcr_combined_LR/combined_res_cons_gpcr.csv')
#Run liana for both conditions
good = NormalizeData(good, normalization.method = "LogNormalize", scale.factor = 10000)
liana_test <-  liana_wrap(good, method = 'cellphonedb', resource = 'custom', external_resource = combined_res)
saveRDS(liana_test, file = 'savefile/control.rds')
poor = NormalizeData(poor, normalization.method = "LogNormalize", scale.factor = 10000)
liana_test2 <- liana_wrap(poor, method = 'cellphonedb', resource = 'custom', external_resource = combined_res)
saveRDS(liana_test2, file = 'savefile/condition.rds')

# For each study separately
study = unique(SeuratObj$Study)
Idents(SeuratObj) = SeuratObj$Study
for (s in study) {
  sub_seuratobj = subset(SeuratObj, idents == s)
  Idents(sub_seuratobj) = sub_seuratobj@meta.data$Annotation_1
  cond = unique(sub_seuratobj$Condition)
  data = as.SingleCellExperiment(sub_seuratobj)
  for (con in cond) {
    gp <-assays(data)[["counts"]][,colData(data)$Condition==cond]
    metadata_gp <-tibble::tibble(rownames(colData(data))[colData(data)$Condition==cond],
                                 as.character(colData(data)$Annotation_1[colData(data)$Condition==cond]))

    #create seurat object for running liana for conditions separately
    good = CreateSeuratObject(counts = gp)
    Idents(good) = metadata_gp$"as.character(...)"
    #Run liana for both conditions
    good = NormalizeData(good, normalization.method = "LogNormalize", scale.factor = 10000)
    liana_test <-  liana_wrap(good, method = 'cellphonedb', resource = 'custom', external_resource = combined_res)
    saveRDS(liana_test, file =paste0('save/file/', con, '.rds'))
    
  }
}

paths <- c('list/of/path/to/diff/study/output')
for (s in seq_along(paths)) {
  p <- paths[s]
  liana_ctrl <- readRDS(file.path(p, 'control.rds'))
  liana_cond <- readRDS(file.path(p, 'condition.rds'))
  
  ctrl <- liana_ctrl %>%
    filter(!is.na(pvalue), pvalue <= 0.05) %>%
    mutate(
      gene_A = gsub("_", "*", ligand.complex),
      gene_B = gsub("_", "*", receptor.complex),
      type_gene_A = "Ligand",
      type_gene_B = "Receptor"
    )
  
  cond <- liana_cond %>%
    filter(!is.na(pvalue), pvalue <= 0.05) %>%
    mutate(
      gene_A = gsub("_", "*", ligand.complex),
      gene_B = gsub("_", "*", receptor.complex),
      type_gene_A = "Ligand",
      type_gene_B = "Receptor"
    )
  
  lr_list <- list(control = ctrl, condition = cond)
  out_file <- paste0("condition_vs_control_", s, ".html")
  generate_report(
    lrpaths     = lr_list,
    genes       = NULL,
    out_path    = p,                    
    threshold   = 0,
    out_file    = out_file,
    output_fmt  = "html_document",
    report      = TRUE,
    sel_columns = c("source","target","gene_A","gene_B","type_gene_A","type_gene_B","lr.mean")
  )
}
