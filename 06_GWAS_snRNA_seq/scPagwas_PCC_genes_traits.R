library(dplyr)
library(tidyverse)
library(ggplot2)
library(grDevices)
library(stats)
library(scales)
library(reshape2)
library(effsize)
library(openxlsx)
library(purrr)

# Read-in and read-out dir
path = "./heart-atlas/scPagwas/outs/"
fig_out = "./heart-atlas/scPagwas/figures/"

# PCC evaluation for different traits
Patients <- c("0_Control_Ch_P1515", "0_Control_Ch_P1516", "0_Control_Ch_P1539", "0_Control_Ch_P1540", "0_Control_Ch_P1547", 
              "0_Control_Ch_P1558", "0_Control_Ch_P1582", "0_Control_Ch_P1600", "0_Control_Ch_P1603", "0_Control_Ch_P1610", "0_Control_Ch_P1678",
              "0_Control_Re_H2", "0_Control_Re_H3", "0_Control_Re_H4", "0_Control_Re_H49", "0_Control_Re_H5",
              "0_Control_Re_H53", "0_Control_Re_H6", "0_Control_Re_H7")

GWAS_Studies <- c("GCST90018834-EFO_0000407", "Jurgens_DCM", "GCST90018861-EFO_0000538", "GCST90435254", "GCST90043672-EFO_0001073", "GCST90268128", "GCST90043949-EFO_0000537", "GCST011365-EFO_0000612")

for(GWAS in GWAS_Studies){

Patient_PCC <- data.frame()
for(Patient in Patients){
  # Load PCC file for the GWAS study
  PCC_file <- paste0(path, "reactome_", Patient, "_", GWAS, "/reactome_", Patient, "_", GWAS, "_gene_PCC.csv")
  
  if (file.exists(PCC_file)) {
    PCC <- read.csv(PCC_file, header = TRUE, row.names = 1)
    PCC$Condition <- paste(strsplit(Patient, "_")[[1]][1:2], collapse = "_")
    PCC$Gene <- row.names(PCC)
    rownames(PCC) <- NULL
    Patient_PCC <- rbind(Patient_PCC, PCC)
  } else {
    print(paste("GWAS", GWAS, "for Patient", Patient, "does not exist."))
  }
}

heritability_cor_scatterplot2 <- function(gene_heri_cor = NULL,
                                          topn_genes_label = 10,
                                          color_low = "#035397",
                                          color_high = "#F32424",
                                          color_mid = "white",
                                          text_size = 5,
                                          do_plot = TRUE,
                                          max.overlaps = 10,
                                          figurenames = NULL,
                                          width = 4,
                                          height = 4) {
  text <- NULL
  cor_df <- data.frame(genes = gene_heri_cor$Gene, PCC = gene_heri_cor$mean_PCC, text = rep(NA, nrow(gene_heri_cor)))
  
  cor_df <- cor_df[order(cor_df$PCC, decreasing = T), ]
  cor_df$text[1:topn_genes_label] <- cor_df$genes[1:topn_genes_label]
  
  cor_df$order <- seq_len(nrow(cor_df))
  
  p <- ggplot(data = cor_df) +
    geom_point(mapping = aes(x = order, y = PCC, color = PCC)) +
    scale_colour_gradient2(
      low = color_low,
      mid = color_mid,
      high = color_high,
      midpoint = 0
    ) +
    theme_classic() +
    theme(axis.text = element_text(size = 16, colour = "black"), axis.title = element_text(size = 21, colour = "black"))+
    ggrepel::geom_text_repel(aes(
      x = order, y = PCC,
      label = text
    ),
    max.overlaps = max.overlaps,
    na.rm = T, size = text_size
    )
  
  if (do_plot) print(p)
  if (!is.null(figurenames)) {
    grDevices::pdf(file = figurenames, width = width, height = height)
    print(p)
    grDevices::dev.off()
  }
}



results <- Patient_PCC %>% group_by(Gene) %>% summarise(mean_PCC = mean(PCC))
top_PCC <- results %>% arrange(desc(mean_PCC)) %>% slice_head(n=20) %>% pull(Gene)
Patient_top_PCC <- Patient_PCC %>% filter(Gene %in% top_PCC)

plot_Patient_top_PCC <- ggplot(Patient_top_PCC, aes(x=reorder(Gene, -PCC), y=PCC, color = PCC))+
        scale_colour_gradient2(
          low = "#035397",
          mid = "white",
          high = "#F32424",
          midpoint = 0,
          limits = range(results$mean_PCC, na.rm = TRUE),
          oob = scales::squish)+
        geom_point(aes(color=PCC))+
        theme_classic()+theme(axis.text = element_text(size = 21, colour = "black"), axis.title = element_text(size = 21, colour = "black"),
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf(paste0(fig_out, GWAS, "_control.pdf"))
heritability_cor_scatterplot2(gene_heri_cor = results, max.overlaps = 20,topn_genes_label = 2)
print(plot_Patient_top_PCC)
dev.off()
}


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
              
              "0_Control_Ch_P1515", "0_Control_Ch_P1516", "0_Control_Ch_P1539", "0_Control_Ch_P1540", "0_Control_Ch_P1547", 
              "0_Control_Ch_P1558", "0_Control_Ch_P1582", "0_Control_Ch_P1600", "0_Control_Ch_P1603", "0_Control_Ch_P1610", "0_Control_Ch_P1678",
              "0_Control_Re_H2", "0_Control_Re_H3", "0_Control_Re_H4", "0_Control_Re_H49", "0_Control_Re_H5",
              "0_Control_Re_H53", "0_Control_Re_H6", "0_Control_Re_H7")

GWAS_Studies <- c("GCST90018834-EFO_0000407", "Jurgens_DCM")
GWAS_union <- list()

for(GWAS in GWAS_Studies){
  
  Patient_PCC <- data.frame()
  for(Patient in Patients){
    # Load PCC file for the GWAS study
    PCC_file <- paste0(path, "reactome_", Patient, "_", GWAS, "/reactome_", Patient, "_", GWAS, "_gene_PCC.csv")
    
    if (file.exists(PCC_file)) {
      PCC <- read.csv(PCC_file, header = TRUE, row.names = 1)
      PCC$Condition <- paste(strsplit(Patient, "_")[[1]][1:2], collapse = "_")
      PCC$Gene <- row.names(PCC)
      rownames(PCC) <- NULL
      Patient_PCC <- rbind(Patient_PCC, PCC)
    } else {
      print(paste("GWAS", GWAS, "for Patient", Patient, "does not exist."))
    }
  }

  #GWAS specific
  Patient_PCC$Condition <- gsub("^[01234]_", "", Patient_PCC$Condition)

  # Split by gene, calculate everything per group
  results <- Patient_PCC %>%
    group_by(Gene) %>%
    group_split() %>%
    map_dfr(function(df_gene) {
      gene_name <- unique(df_gene$Gene)
      group1 <- df_gene$PCC[df_gene$Condition == "DCM"]
      group2 <- df_gene$PCC[df_gene$Condition == "Control"]
      
      # Run safely
      test <- tryCatch(t.test(group1, group2), error = function(e) NULL)
      d <- tryCatch(cohen.d(PCC ~ Condition, data = df_gene)$estimate, error = function(e) NA_real_)
      g <- tryCatch(hedges.g(PCC ~ Condition, data = df_gene)$estimate, error = function(e) NA_real_)
      glass <- tryCatch(glass.delta(PCC ~ Condition, data = df_gene)$estimate, error = function(e) NA_real_)
      
      tibble(
        Gene = gene_name,
        mean_DCM = mean(group1, na.rm = TRUE),
        mean_Control = mean(group2, na.rm = TRUE),
        log2FC = log2(mean(group1, na.rm = TRUE) / mean(group2, na.rm = TRUE)),
        p_value = if (!is.null(test)) test$p.value else NA_real_,
        cohen_d = d,
        hedges_g = g,
        glass_d = glass
      )
    }) %>%
    mutate(
      adj_p = p.adjust(p_value, method = "BH"),
      log10_adj_p = -log10(adj_p)
    )
  
  
  write.csv(results, file=paste0(path, "reactome_", GWAS, "summary_gene_PCC.csv"))
  
  top_GWAS_genes <- results %>% 
    arrange(desc(mean_Control)) %>% 
    #filter(adj_p < 0.08) %>%
    filter(log2FC < 0) %>%
    mutate(rank = "top", GWAS = GWAS) #%>%
    #slice_head(n = 30)
  
  bottom_GWAS_genes <- results %>%
    arrange(desc(mean_DCM)) %>%
    #filter(adj_p < 0.08) %>%
    filter(log2FC > 0) %>%
    mutate(rank = "bottom", GWAS = GWAS) #%>%
    #slice_head(n = 30)
  #top_GWAS_genes$Gene <- factor(top_GWAS_genes$Gene, levels = rev(top_GWAS_genes$Gene))

  GWAS_union[[GWAS]] <- bind_rows(top_GWAS_genes, bottom_GWAS_genes)

  pdf(paste0(fig_out, GWAS, "_PCC_genes_comparision.pdf"))
  print(ggplot(top_GWAS_genes, aes(y=reorder(Gene, mean_Control), x=mean_Control, color=log2FC, size=adj_p))+
    geom_point()+
    scale_color_gradient2(low = "darkred", mid = "red", high = "orange")+
    labs(y="", fill="log2FC", x="mean PCC - Control")+
    theme_minimal() +
    theme(
      axis.title = element_text(size=14),
      axis.text = element_text(size=10),
      axis.ticks.y = element_line(),
      panel.grid.minor = element_blank()
    ))
  
  print(ggplot(bottom_GWAS_genes, aes(y=reorder(Gene, mean_DCM), x=mean_DCM, color=log2FC, size=adj_p))+
          geom_point()+
          scale_color_gradient2(low = "darkred", mid = "red", high = "orange")+
          labs(y="", fill="log2FC", x="mean PCC - DCM")+
          theme_minimal() +
          theme(
            axis.title = element_text(size=14),
            axis.text = element_text(size=10),
            axis.ticks.y = element_line(),
            panel.grid.minor = element_blank()
          ))
  dev.off()
}


####
combined_GWAS <- bind_rows(GWAS_union)

gene_counts <- combined_GWAS %>%
  group_by(Gene, rank) %>%
  summarise(n_GWAS = n_distinct(GWAS), .groups = "drop") %>%
  filter(n_GWAS > 1) %>%
  mutate(label = "both") %>%
  select(Gene, rank, label)

overlapping_GWAS <- combined_GWAS %>%
  left_join(gene_counts, by = c("Gene", "rank")) %>%
  mutate(label = ifelse(is.na(label), GWAS, label))

top_GWAS_genes <- overlapping_GWAS %>%
  filter(label == "both") %>%
  group_by(Gene) %>%
  summarise(sorting_mean = mean(mean_Control)) %>%
  arrange(desc(sorting_mean)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  slice_head(n = 50) %>%
  pull(Gene)

top50_GWAS <- overlapping_GWAS %>%
  filter(Gene %in% top_GWAS_genes) %>%
  filter(GWAS == "Jurgens_DCM")


pdf(paste0(fig_out, "_DCM_GWAS_PCC_genes_comparision.pdf"), height = 8, width = 6)
print(ggplot(top50_GWAS, aes(y=reorder(Gene, mean_Control), x=mean_Control, color=log2FC, size=log10_adj_p))+
        geom_point()+#xlim(0.05,max(top50_GWAS$mean_Control))+
        scale_color_gradient2(low = "darkred", mid = "red", high = "orange", midpoint = -2)+
        labs(y="", color="log2FC\nDCM\nvs.\nCTRL", x="mean (PCC - Control)", size="-log10\n(adj.pval)")+
        theme_minimal() +
        theme(
          axis.title = element_text(size=14),
          axis.text = element_text(size=10, color = "black"),
          axis.ticks.y = element_line(),
          panel.grid.minor = element_blank()
        ))

print(ggplot(top50_GWAS %>% filter(adj_p <0.05), aes(y=reorder(Gene, mean_Control), x=mean_Control, color=log2FC, size=log10_adj_p))+
        geom_point()+xlim(0.05,max(top50_GWAS$mean_Control))+
        scale_color_gradient2(low = "darkred", mid = "red", high = "orange", midpoint = -1.5)+
        labs(y="", color="log2FC\nDCM\nvs.\nCTRL", x="mean (PCC - Control)", size="-log10\n(adj.pval)")+
        theme_minimal() +
        theme(
          axis.title = element_text(size=14),
          axis.text = element_text(size=10, color = "black"),
          axis.ticks.y = element_line(),
          panel.grid.minor = element_blank()
        ))

dev.off()

### HCM

Patients <- c("HCM_Ch_P1422", "HCM_Ch_P1425", "HCM_Ch_P1447", "HCM_Ch_P1462", "HCM_Ch_P1479", "HCM_Ch_P1508",
              "HCM_Ch_P1510", "HCM_Ch_P1602", "HCM_Ch_P1630", "HCM_Ch_P1631", "HCM_Ch_P1685", "HCM_Ch_P1707", "HCM_Ch_P1722", "HCM_Ch_P1726", "HCM_Ch_P1735",
              
              "0_Control_Ch_P1515", "0_Control_Ch_P1516", "0_Control_Ch_P1539", "0_Control_Ch_P1540", "0_Control_Ch_P1547", 
              "0_Control_Ch_P1558", "0_Control_Ch_P1582", "0_Control_Ch_P1600", "0_Control_Ch_P1603", "0_Control_Ch_P1610", "0_Control_Ch_P1678",
              "0_Control_Re_H2", "0_Control_Re_H3", "0_Control_Re_H4", "0_Control_Re_H49", "0_Control_Re_H5",
              "0_Control_Re_H53", "0_Control_Re_H6", "0_Control_Re_H7")

GWAS_Studies <- c("GCST90435254")

for(GWAS in GWAS_Studies){
  
  Patient_PCC <- data.frame()
  for(Patient in Patients){
    # Load PCC file for the GWAS study
    PCC_file <- paste0(path, "reactome_", Patient, "_", GWAS, "/reactome_", Patient, "_", GWAS, "_gene_PCC.csv")
    
    if (file.exists(PCC_file)) {
      PCC <- read.csv(PCC_file, header = TRUE, row.names = 1)
      PCC$Condition <- paste(strsplit(Patient, "_")[[1]][1:2], collapse = "_")
      PCC$Gene <- row.names(PCC)
      rownames(PCC) <- NULL
      Patient_PCC <- rbind(Patient_PCC, PCC)
    } else {
      print(paste("GWAS", GWAS, "for Patient", Patient, "does not exist."))
    }
  }
  
  #GWAS specific
  Patient_PCC$Condition <- gsub("^[01234]_", "", Patient_PCC$Condition)
  Patient_PCC$Condition[Patient_PCC$Condition == "HCM_Ch"] <- "HCM"
  
  results <- Patient_PCC %>%
    group_by(Gene) %>%
    summarise(
      p_value = t.test(PCC ~ Condition, data = cur_data())$p.value,
      mean_HCM = mean(PCC[Condition == "HCM"]),
      mean_Control = mean(PCC[Condition == "Control"]),
      cohen_d = tryCatch(
        cohen.d(PCC ~ Condition)$estimate,
        error = function(e) NA
      )
    ) %>%
    mutate(log2FC = log2(mean_HCM / mean_Control))
  results <- results %>%
    mutate(
      adj_p = p.adjust(p_value, method = "BH"),
      log10_adj_p = -log10(adj_p)
    )
  
  write.csv(results, file=paste0(path, "reactome_", GWAS, "summary_gene_PCC.csv"))
  
  top_GWAS_genes <- results %>% 
    arrange(desc(mean_Control)) %>% 
    filter(log2FC < 0) %>%
    mutate(rank = "top", GWAS = GWAS)
  
  bottom_GWAS_genes <- results %>%
    arrange(desc(mean_HCM)) %>%
    filter(log2FC > 0) %>%
    mutate(rank = "bottom", GWAS = GWAS)
  
  #top_GWAS_genes <- results %>% 
  #  arrange(desc(mean_Control)) %>% 
  #  filter(p_value < 0.05) %>%
  #  filter(log2FC > 0) %>%
  #  slice_head(n = 30) #[1:30, ]
  
  #bottom_GWAS_genes <- results %>%
  #  arrange(desc(mean_HCM)) %>%
  #  filter(p_value < 0.05) %>%
  #  filter(log2FC < 0) %>%
  #  slice_head(n = 30)
  #top_GWAS_genes$Gene <- factor(top_GWAS_genes$Gene, levels = rev(top_GWAS_genes$Gene))
  
  combined <- bind_rows(
    top_GWAS_genes %>% mutate(Category = "Top"),
    bottom_GWAS_genes %>% mutate(Category = "Bottom")
  )
  
  pdf(paste0(fig_out, "HCM_", GWAS, "_PCC_genes_comparision.pdf"))
  print(ggplot(top_GWAS_genes, aes(y=reorder(Gene, mean_Control), x=mean_Control, color=log2FC, size=adj_p))+
          geom_point()+
          scale_color_gradient2(low = "darkred", mid = "red", high = "orange", midpoint = (-1.5))+
          labs(y="", fill="log2FC", x="mean PCC - Control")+
          theme_minimal() +
          theme(
            axis.title = element_text(size=14),
            axis.text = element_text(size=10),
            axis.ticks.y = element_line(),
            panel.grid.minor = element_blank()
          ))
  
  print(ggplot(bottom_GWAS_genes, aes(y=reorder(Gene, mean_HCM), x=mean_HCM, color=log2FC, size=adj_p))+
          geom_point()+
          scale_color_gradient2(low = "darkred", mid = "red", high = "orange", midpoint = 1)+
          labs(y="", fill="log2FC", x="mean PCC - DCM")+
          theme_minimal() +
          theme(
            axis.title = element_text(size=14),
            axis.text = element_text(size=10),
            axis.ticks.y = element_line(),
            panel.grid.minor = element_blank()
          ))
  dev.off()
}

top_GWAS_genes <- combined %>%
  group_by(Gene) %>%
  filter(adj_p < 0.05) %>%
  summarise(sorting_mean = mean(mean_Control)) %>%
  arrange(desc(sorting_mean)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  slice_head(n = 40) %>%
  pull(Gene)

top50_GWAS <- combined %>%
  filter(Gene %in% top_GWAS_genes)


pdf(paste0(fig_out, "HCM_GWAS_PCC_genes_comparision.pdf"), height = 6, width = 4)
print(ggplot(top50_GWAS, aes(y=reorder(Gene, mean_Control), x=mean_Control, color=log2FC, size=log10_adj_p))+
        geom_point()+
        scale_color_gradient2(low = "darkred", mid = "red", high = "orange", midpoint = 2)+
        labs(y="", color="log2FC\nDCM\nvs.\nCTRL", x="mean (PCC - Control)", size="-log10\n(adj.pval)")+
        theme_minimal() +
        theme(
          axis.title = element_text(size=12),
          axis.text = element_text(size=10, color = "black"),
          axis.ticks.y = element_line(),
          panel.grid.minor = element_blank()
        ))
dev.off()


