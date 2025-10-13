library(scPagwas)
library(Seurat)
library(readr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
GWAS = args[1]
directory <- "./heart-atlas/scPagwas/data"
files <- list.files(path = directory, pattern = "^[0123].*\\.rds$", full.names = TRUE)


for (data in files){
  Single_data <- readRDS(data)
  Patient <- sub("^.*/data/.*/data/(.*)\\.rds$", "\\1", data)
  Single_data <- NormalizeData(Single_data, normalization.method = "LogNormalize", scale.factor = 10000)
  Single_data <- ScaleData(Single_data)

  print(GWAS)
  GWAS_study <- read_table(file = paste0("./heart-atlas/scPagwas/GWAS/Studies/", GWAS, ".h.tsv.gz"))
  if (GWAS == "Jurgens_DCM"){
    GWAS_study <- GWAS_study %>% rename("chrom"=CHR, "pos"=POS, "REF"=EA, "ALT"=NEA, "beta"=BETA, "se"=SE, "maf"=EAFREQ, "p_value"=P, "rsid"=rsID)
  } else {
    GWAS_study <- GWAS_study %>% rename("chrom"=chromosome, "pos"=base_pair_location, "REF"=effect_allele, "ALT"=other_allele, "beta"=beta, "se"=standard_error, "maf"=effect_allele_frequency, "p_value"=p_value)
    if (!"rsid" %in% colnames(GWAS_study)) {
      GWAS_study <- GWAS_study %>%
        rename(rsid = variant_id)
  }
   }
  GWAS_study <- GWAS_study %>% select("chrom", "pos", "REF", "ALT", "beta", "se", "maf", "p_value", "rsid")
  GWAS_study <- data.frame(GWAS_study)

  print(GWAS)
  print(head(GWAS_study))

  Pagwas <- list()
  Pagwas <- Single_data_input(Pagwas=Pagwas,
                            assay="RNA",
                            Single_data=Single_data,
                            Pathway_list=reduce_genes.by.reactome.pathway)

  Pagwas <- scPagwas_main(Pagwas =Pagwas,
                                  gwas_data = GWAS_study,
                                  Single_data = Single_data,
                                  output.prefix=paste("reactome", Patient, GWAS, sep="_"),
                                  output.dirs=paste("reactome", Patient, GWAS, sep="_"),
                                  Pathway_list=reduce_genes.by.reactome.pathway,
                                  assay="RNA",
                                  block_annotation = block_annotation,
                                  iters_singlecell = 100,
                                  chrom_ld = chrom_ld,
                                  singlecell=T,
                                  celltype=F)
}