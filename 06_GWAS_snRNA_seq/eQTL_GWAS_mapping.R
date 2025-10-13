library(ggplot2)
library(eQTpLot)
library(tidyverse)
library(biomaRt)
library(rtracklayer)

setwd("./heart-atlas/scPagwas/figures")

#Matching reference genome for GTEX v8 (ENCODE v26)
gtf <- import("./heart-atlas/scPagwas/GWAS/Studies/gencode.v26.chr_patch_hapl_scaff.basic.annotation.gtf.gz")
genes_df <- as.data.frame(gtf[gtf$type == "gene", c("gene_id", "gene_name")])

#Prepare DCM 2 GWAS
GWAS.df <- read.table("./heart-atlas/scPagwas/GWAS/Studies/Jurgens_DCM.h.tsv.gz", header=T)

# same variant ID
GWAS.df$SNP <- paste0("chr", gsub(":", "_", GWAS.df$ID_B38), "_b38")

#GWAS_study <- GWAS_study %>% dplyr::rename("chrom"=CHR, "pos"=POS, "effect_allele"=EA, "other_allele"=NEA, "beta"=BETA, "se"=SE, "maf"=EAFREQ, "p"=P, "rsid"=rsID)
GWAS.df <- GWAS.df %>% dplyr::rename("BP"=POS) %>%
    dplyr::select("CHR", "BP", "SNP", "P", "BETA")
GWAS.df <- data.frame(GWAS.df)

Genes.df <- genes_df %>% dplyr::rename("CHR" = seqnames, "Start" = start, "Stop" = end, "Gene" = gene_name) %>%
  dplyr::select(CHR, Start, Stop, Gene)
Genes.df$CHR <- sub("^chr", "", Genes.df$CHR)
Genes.df$Build = "hg38"
Genes.df$CHR <- as.numeric(Genes.df$CHR)
Genes.df$Start <- as.numeric(Genes.df$Start)
Genes.df$Stop <- as.numeric(Genes.df$Stop)



# Prepare eQTL data for different tissues downloaded from: https://gtexportal.org/home/downloads/adult-gtex/qtl
eQTL.df <- read.table("./heart-atlas/scPagwas/GWAS/Studies/GTEx_Analysis_v8_eQTL/Artery_Tibial.v8.signif_variant_gene_pairs.txt.gz", header=T)

eqtl_annotated <- merge(
  eQTL.df,
  genes_df[, c("gene_id", "gene_name")],
  by = "gene_id",
  all.x = TRUE
)
eqtl_annotated$Tissue = "Artery - Tibial"

eQTL.df <- eqtl_annotated %>% dplyr::rename("Gene.Symbol"=gene_name, "SNP.Id"=variant_id, "P.Value" = pval_nominal, "NES" = slope) %>%
  dplyr::select("Gene.Symbol", "SNP.Id", "P.Value", "NES", "Tissue")
eQTL.df <- data.frame(eQTL.df)

eQTpLot(GWAS.df = GWAS.df, eQTL.df = eQTL.df, gene = c("MLIP","ACTN2"), 
        gbuild = "hg38", trait ="DCM", tissue = "Artery - Tibial", CollapseMethod = "min", 
        GeneList = T, Genes.df = Genes.df, getplot=FALSE) 


pdf("./heart-atlas/scPagwas/figures/MLIP_eQTpLot.pdf") #, width = 1000, height = 800, res = 150)
  eQTpLot(GWAS.df = GWAS.df,
          eQTL.df = eQTL.df,
          Genes.df = Genes.df,
          gbuild = "hg38",
          gene = c("MLIP"),
          tissue = "Artery - Tibial",
          trait = "DCM")
dev.off()


eQTL.df <- read.table("./heart-atlas/scPagwas/GWAS/Studies/GTEx_Analysis_v8_eQTL/Muscle_Skeletal.v8.signif_variant_gene_pairs.txt.gz", header=T)

eqtl_annotated <- merge(
  eQTL.df,
  genes_df[, c("gene_id", "gene_name")],
  by = "gene_id",
  all.x = TRUE
)
eqtl_annotated$Tissue = "Muscle_Skeletal"

eQTL.df <- eqtl_annotated %>% dplyr::rename("Gene.Symbol"=gene_name, "SNP.Id"=variant_id, "P.Value" = pval_nominal, "NES" = slope) %>%
  dplyr::select("Gene.Symbol", "SNP.Id", "P.Value", "NES", "Tissue")
eQTL.df <- data.frame(eQTL.df)

eQTpLot(GWAS.df = GWAS.df, eQTL.df = eQTL.df, gene = c("TMEM178B","ACTN2"),
        gbuild = "hg38", trait ="DCM", tissue = "Muscle_Skeletal", CollapseMethod = "min", 
        GeneList = T, Genes.df = Genes.df, getplot=FALSE) 


pdf("./heart-atlas/scPagwas/figures/TMEM178B_eQTpLot.pdf") #, width = 1000, height = 800, res = 150)
eQTpLot(GWAS.df = GWAS.df,
        eQTL.df = eQTL.df,
        Genes.df = Genes.df,
        gbuild = "hg38",
        gene = c("TMEM178B"),
        tissue = "Muscle_Skeletal",
        trait = "DCM")
dev.off()

eQTpLot(GWAS.df = GWAS.df, eQTL.df = eQTL.df, gene = c("FHOD3","ACTN2"), 
        gbuild = "hg38", trait ="DCM", tissue = "Muscle_Skeletal", CollapseMethod = "min", 
        GeneList = T, Genes.df = Genes.df, getplot=FALSE) 


pdf("./heart-atlas/scPagwas/figures/FHOD3_eQTpLot_muscle.pdf") #, width = 1000, height = 800, res = 150)
eQTpLot(GWAS.df = GWAS.df,
        eQTL.df = eQTL.df,
        Genes.df = Genes.df,
        gbuild = "hg38",
        gene = c("FHOD3"),
        tissue = "Muscle_Skeletal",
        trait = "DCM")
dev.off()


eQTL.df <- read.table("./heart-atlas/scPagwas/GWAS/Studies/GTEx_Analysis_v8_eQTL/Heart_Atrial_Appendage.v8.signif_variant_gene_pairs.txt.gz", header=T)

eqtl_annotated <- merge(
  eQTL.df,
  genes_df[, c("gene_id", "gene_name")],
  by = "gene_id",
  all.x = TRUE
)
eqtl_annotated$Tissue = "Heart_Atrial_Appendage"

eQTL.df <- eqtl_annotated %>% dplyr::rename("Gene.Symbol"=gene_name, "SNP.Id"=variant_id, "P.Value" = pval_nominal, "NES" = slope) %>%
  dplyr::select("Gene.Symbol", "SNP.Id", "P.Value", "NES", "Tissue")
eQTL.df <- data.frame(eQTL.df)

eQTpLot(GWAS.df = GWAS.df, eQTL.df = eQTL.df, gene = c("FHOD3","ACTN2"),
        gbuild = "hg38", trait ="DCM", tissue = "Heart_Atrial_Appendage", CollapseMethod = "min", 
        GeneList = T, Genes.df = Genes.df, getplot=FALSE) 


pdf("./heart-atlas/scPagwas/figures/FHOD3_eQTpLot_heart.pdf") #, width = 1000, height = 800, res = 150)
eQTpLot(GWAS.df = GWAS.df,
        eQTL.df = eQTL.df,
        Genes.df = Genes.df,
        gbuild = "hg38",
        gene = c("FHOD3"),
        tissue = "Heart_Atrial_Appendage",
        trait = "DCM")
dev.off()


eQTL.df <- read.table("./heart-atlas/scPagwas/GWAS/Studies/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz", header=T)


eqtl_annotated <- merge(
  eQTL.df,
  genes_df[, c("gene_id", "gene_name")],
  by = "gene_id",
  all.x = TRUE
)
eqtl_annotated$Tissue = "Whole_Blood"

eQTL.df <- eqtl_annotated %>% dplyr::rename("Gene.Symbol"=gene_name, "SNP.Id"=variant_id, "P.Value" = pval_nominal, "NES" = slope) %>%
  dplyr::select("Gene.Symbol", "SNP.Id", "P.Value", "NES", "Tissue")
eQTL.df <- data.frame(eQTL.df)

eQTpLot(GWAS.df = GWAS.df, eQTL.df = eQTL.df, gene = c("MYBPC3","ACTN2"),
        gbuild = "hg38", trait ="DCM", tissue = "Whole_Blood", CollapseMethod = "min", 
        GeneList = T, Genes.df = Genes.df, getplot=FALSE) 


pdf("./heart-atlas/scPagwas/figures/MYBPC3_eQTpLot_heart.pdf") #, width = 1000, height = 800, res = 150)
eQTpLot(GWAS.df = GWAS.df,
        eQTL.df = eQTL.df,
        range  = 200,
        Genes.df = Genes.df,
        gbuild = "hg38",
        gene = c("MYBPC3"),
        tissue = "Whole_Blood",
        trait = "DCM")
dev.off()
