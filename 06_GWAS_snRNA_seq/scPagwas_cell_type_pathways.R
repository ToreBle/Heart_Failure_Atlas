library(Seurat)
library(dplyr)
library(tidyr)
library(reshape2)
library(purrr)
library(ggplot2)
library(ggh4x)

prepare_pagwas_dotplot_data_from_df <- function(rds_df,
                                                pvalue_threshold = 0.05,
                                                max_log10_pvalue = 15,
                                                size_column_name = "CellpValue",
                                                color_column_name = "proportion") {
  
  add_specificity <- function(sc_rankp_matrix) {
    pathway_df <- as.data.frame(sc_rankp_matrix)
    pathway_df$path <- rownames(pathway_df)
    
    pathway_df_long <- pathway_df %>%
      add_count(path) %>%
      filter(n == 1) %>%
      select(-n) %>%
      pivot_longer(-path, names_to = "celltype", values_to = "rankp") %>%
      as_tibble() %>%
      group_by(celltype) %>%
      mutate(paras_sum_mean = (rankp * nrow(sc_rankp_matrix)) / sum(rankp)) %>%
      group_by(path) %>%
      mutate(specificity = (paras_sum_mean * nrow(sc_rankp_matrix)) / sum(paras_sum_mean)) %>%
      ungroup()
    
    return(pathway_df_long)
  }
  
  combined_plot_data <- list()
  
  for (i in seq_len(nrow(rds_df))) {
    seurat_path     <- rds_df$rds_path[i]
    condition_label <- rds_df$condition[i]
    
    pagwas_obj <- readRDS(seurat_path)
    celltypes   <- unique(Idents(pagwas_obj))
    
    proportion_by_celltype <- tapply(
      as.vector(Idents(pagwas_obj)),
      Idents(pagwas_obj),
      function(x) {
        her_mat <- t(GetAssayData(pagwas_obj, assay = "scPagwasPaHeritability"))
        apply(her_mat, 2, function(y) sum(y > 0) / length(y))
      }
    )
    
    proportion_df <- Reduce(cbind, proportion_by_celltype)
    colnames(proportion_df) <- names(proportion_by_celltype)
    rownames(proportion_df) <- gsub("-", "_", rownames(proportion_df))
    
    rankp_matrix <- -log10(pagwas_obj@misc$scPathways_rankPvalue + 1e-20)
    
    valid_pathways <- intersect(rownames(proportion_df), rownames(rankp_matrix))
    proportion_df <- proportion_df[valid_pathways, celltypes, drop = FALSE]
    rankp_matrix <- rankp_matrix[valid_pathways, celltypes, drop = FALSE]
    
    rankp_matrix[rankp_matrix > max_log10_pvalue]               <- max_log10_pvalue
    rankp_matrix[rankp_matrix < -log10(pvalue_threshold)]       <- 0
    
    proportion_df <- as.data.frame(proportion_df, stringsAsFactors = FALSE)
    proportion_df$pathways <- rownames(proportion_df)
    
    rankp_df <- as.data.frame(rankp_matrix, stringsAsFactors = FALSE)
    rankp_df$pathways    <- rownames(rankp_df)
    
    melted_proportion <- reshape2::melt(proportion_df,
                                        id.vars = "pathways",
                                        variable.name = "celltype",
                                        value.name    = color_column_name)
    
    melted_rankp <- reshape2::melt(rankp_df,
                                   id.vars = "pathways",
                                   variable.name = "celltype",
                                   value.name    = size_column_name)
    
    merged_df <- merge(melted_rankp, melted_proportion,
                       by = c("pathways", "celltype"))
    
    specificity_df <- add_specificity(rankp_matrix)
    merged_df <- left_join(merged_df, specificity_df,
                           by = c("pathways" = "path", "celltype"))
    
    merged_df$condition <- condition_label
    merged_df$patient   <- unique(pagwas_obj$Donor)
    
    combined_plot_data[[ length(combined_plot_data) + 1 ]] <- merged_df
  }
  
  final_dotplot_df <- bind_rows(combined_plot_data)
  return(final_dotplot_df)
}


truncate_pathway <- function(pathway, max_length = 50, prefix_to_remove = NULL) {
  # Remove prefix if provided
  if (!is.null(prefix_to_remove)) {
    pathway <- gsub(paste0("^", prefix_to_remove), "", pathway)
  }
  
  # If within length, return as-is
  if (nchar(pathway) <= max_length) {
    return(pathway)
  }
  
  # Truncate at nearest underscore
  underscores <- gregexpr("_", pathway)[[1]]
  cut_pos <- max(underscores[underscores <= max_length])
  
  if (is.finite(cut_pos)) {
    return(paste0(substr(pathway, 1, cut_pos - 1), "..."))
  } else {
    return(paste0(substr(pathway, 1, max_length), "..."))
  }
}

###

plot_pathways_for_celltype <- function(df,
                                       cell_state,
                                       preselect_condition = "ctrl",
                                       fc_conditions       = c("ctrl","dcm"),
                                       n_ctrl              = 120,
                                       n_up                = 20,
                                       n_down              = 20,
                                       max_label_len       = 50,
                                       prefix_to_remove    = "REACTOME_",
                                       effect_measure      = c("log2FC","cohens_d"),
                                       pseudocount         = 1e-6) {
  effect_measure <- match.arg(effect_measure)
  
  # 1) Top pathways by cell-level specificity
  top_ctrl <- df %>%
    filter(celltype == cell_state,
           condition == preselect_condition) %>%
    group_by(pathways) %>%
    summarise(mean_spec = mean(specificity, na.rm=TRUE),
              .groups = "drop") %>%
    arrange(desc(mean_spec)) %>%
    slice_head(n = n_ctrl) %>%
    pull(pathways)

  
  # 2) Compute log2FC, Cohen's d and p-values
  stats_df <- df %>%
    filter(celltype == cell_state,
           condition %in% fc_conditions,
           pathways  %in% top_ctrl) %>%
    group_by(pathways) %>%
    summarise(
      m1 = mean(proportion[condition == fc_conditions[1]], na.rm=TRUE),
      m2 = mean(proportion[condition == fc_conditions[2]], na.rm=TRUE),
      v1 = list(proportion[condition == fc_conditions[1]]),
      v2 = list(proportion[condition == fc_conditions[2]]),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      log2FC    = log2((m2 + pseudocount)/(m1 + pseudocount)),
      cohens_d  = tryCatch(cohen.d(unlist(v2), unlist(v1), pooled=TRUE)$estimate,
                           error=function(e) NA_real_),
      p_value   = tryCatch(wilcox.test(unlist(v1), unlist(v2), paired=FALSE)$p.value,
                           error=function(e) NA_real_),
      p_adj     = p.adjust(p_value, method="fdr"),
      log10_padj = -log10(p_adj),
      effect_size = if (effect_measure=="log2FC") log2FC else cohens_d,
      regulation  = ifelse(effect_size > 0, "Up", "Down")
    ) %>%
    ungroup()
  
  # 3) Top-Up and Top-Down
  shorten <- function(x) sapply(x,
                                truncate_pathway,
                                max_length=max_label_len,
                                prefix_to_remove=prefix_to_remove)
  
  top_up <- stats_df %>%
    filter(regulation=="Up") %>%
    arrange(desc(effect_size)) %>%
    slice_head(n = n_up) %>%
    mutate(short = shorten(pathways)) %>%
    filter(m1>0.5)
    #filter(p_adj <0.05)
  
  
  
  top_down <- stats_df %>%
    filter(regulation=="Down") %>%
    arrange(effect_size) %>%
    slice_head(n = n_down) %>%
    mutate(short = shorten(pathways)) %>%
    filter(m1>0.5)
    #filter(p_adj <0.05)
  
  both <- bind_rows(top_up, top_down) %>%
    filter(p_adj <0.05)
  
  # 4) “No data” placeholder
  no_data <- function(title) {
    ggplot() +
      annotate("text", x=0.5, y=0.5, label=title, size=6, inherit.aes=FALSE) +
      theme_void()
  }
  
  # 5) Three “ctrl-style” plots (RdBu palette)
  make_ctrl <- function(df_sub, title) {
    if (nrow(df_sub)==0) return(no_data(title))
    ggplot(df_sub,
           aes(y=reorder(short, m1),
               x=m1,
               size=-log10(p_adj),
               color=effect_size)) +
      geom_point(alpha=0.8, na.rm=TRUE) +
      scale_color_distiller(name=effect_measure,
                            palette="Reds",
                            direction=-1,
                            guide=guide_colorbar(reverse=TRUE)) +
      labs(title=title,
           x=paste0("Mean prop.(", fc_conditions[1], ")"),
           size="adj. p-value") +
      theme_minimal(base_size=12)
  }
  
  ctrl_up    <- make_ctrl(top_up,    paste(cell_state, "- top Up"))
  ctrl_down  <- make_ctrl(top_down,  paste(cell_state, "- top Down"))
  ctrl_facet <- {
    title <- paste(cell_state, "- Up & Down")
    if (nrow(both)==0) {
      no_data(title)
    } else {
      ggplot(both,
             aes(y=reorder(short, m1),
                 x=m1,
                 size=log10_padj,
                 color=effect_size)) +
        geom_point(alpha=0.8, na.rm=TRUE) +
        facet_grid(regulation ~ ., scales="free_y", space="free_y") +
        scale_color_distiller(name=effect_measure,
                              palette="RdBu",
                              direction=-1,
                              guide=guide_colorbar(reverse=TRUE)) +
        labs(title=title,
             x=paste0("Mean proportion (", fc_conditions[1], ")"),
             size="-log10\n(adj. p-value)",
             color="log2FC\nDCM\nvs.\nctrl",
             y="") +
        theme_minimal(base_size=12)
    }
  }
  
  # 6) Dot-plot data: mean per (celltype, pathway, condition)
  dot_base <- df %>%
    filter(pathways  %in% both$pathways,
           condition %in% fc_conditions) %>%
    group_by(celltype, pathways, condition) %>%
    summarise(
      mean_prop = mean(proportion, na.rm=TRUE),
      mean_pval = mean(CellpValue, na.rm=TRUE),
      .groups   = "drop"
    ) %>%
    left_join(both %>% select(pathways, regulation, short),
              by="pathways") %>%
    filter(!is.na(regulation))
  
  # 7) Dot-plot helper without jitter
  make_dot <- function(subdf, title, facet_by_reg=FALSE) {
    if (nrow(subdf)==0) return(no_data(title))
    p <- ggplot(subdf,
                aes(x=condition,
                    y=short,
                    size=mean_pval,
                    color=mean_prop)) +
      geom_point(alpha=0.8, na.rm=TRUE) +
      scale_color_distiller(name="mean prop.",
                            palette="RdBu",
                            direction=-1,
                            guide=guide_colorbar(reverse=TRUE)) +
      scale_size_continuous(range=c(1,4),
                            guide=guide_legend(title="mean Cell p-val",
                                               title.position="top",
                                               title.hjust=0.5,
                                               reverse=TRUE)) +
      labs(title=title, x=NULL, y=NULL) +
      theme_minimal(base_size=12) +
      theme(axis.text.x=element_text(angle=90, vjust=0.5),
            legend.position="top",
            legend.direction="horizontal")
    if (facet_by_reg) {
      Cell_type_colors <- c(
        "Adipocytes" = '#1f77b4', "Cardiomyocytes" = '#ff7f0e', "Endocardial" = '#279e68', "Endothelial" = '#d62728', 
          "Epicardium" = '#aa40fc', "Fibroblast" = '#8c564b',"Ischemic cells (MI)" = '#e377c2', "Lymphatic EC" = '#b5bd61', 
          "Lymphocytes" = '#17becf', "Mast cells" = '#aec7e8', "Myeloid" = '#ffbb78', "Neuronal" = '#98df8a',
          "Pericytes" = '#ff9896', "VSMC" = '#c5b0d5'
      )
      subdf$celltype <- factor(subdf$celltype, levels = names(Cell_type_colors))
      
      # Create named list of strip themes
      strip_colors <- lapply(Cell_type_colors[levels(subdf$celltype)], function(col) {
        element_rect(fill = col, color = "white", linewidth = 0.5)
      })
      names(strip_colors) <- levels(subdf$celltype)
      
      p + ggh4x::facet_grid2(regulation ~ celltype, scales = "free_y", space = "free_y",
                                 strip = ggh4x::strip_themed(
                                   background_x = strip_colors,
                                   text_x = list(element_text(color = strip_colors))
                                 ))
    } else {
      p + facet_grid(. ~ celltype,
                     scales="free_y", space="free_y")
    }
  }
  
  dot_up    <- make_dot(filter(dot_base, regulation=="Up"),
                        paste(cell_state, "- dot Up"),
                        facet_by_reg=FALSE)
  dot_down  <- make_dot(filter(dot_base, regulation=="Down"),
                        paste(cell_state, "- dot Down"),
                        facet_by_reg=FALSE)
  dot_facet <- make_dot(dot_base,
                        paste(cell_state, "- dot Up & Down"),
                        facet_by_reg=TRUE)
  
  # 8) Return six plots
  list(
    ctrl_up    = ctrl_up,
    ctrl_down  = ctrl_down,
    ctrl_facet = ctrl_facet,
    dot_up     = dot_up,
    dot_down   = dot_down,
    dot_facet  = dot_facet
  )
}

# ── Running the GWAS ──

for(GWAS in c("Jurgens_DCM", "GCST90018834-EFO_0000407", "GCST90043672-EFO_0001073", "GCST90043949-EFO_0000537")){ 
  rds_ctrl <- data.frame(
    rds_path = c(list.files(path = "./heart-atlas/scPagwas/outs/reactome", pattern = paste0("^0_Control.*\\", GWAS, ".rds$"), full.names = TRUE)),
    condition = c("ctrl"),
    stringsAsFactors = FALSE
  )
  
  rds_dcm <- data.frame(
    rds_path = c(list.files(path = "./heart-atlas/scPagwas/outs/reactome", pattern = paste0("^[01234]_DCM.*\\", GWAS, ".rds$"), full.names = TRUE)),
    condition = c("dcm"),
    stringsAsFactors = FALSE
  )
  
  rds_info <- rbind(rds_dcm, rds_ctrl)
  
  final_df <- prepare_pagwas_dotplot_data_from_df(rds_info)
  
  #print(final_df %>% filter(celltype == "Cardiomyocytes") %>% filter(pathways == "REACTOME_AUTOPHAGY"))
  #print(final_df %>% filter(celltype == "Cardiomyocytes") %>% filter(pathways =="REACTOME_PHASE_0_RAPID_DEPOLARISATION"))
  #print(final_df %>% filter(celltype == "Cardiomyocytes") %>% filter(pathways == "REACTOME_MUSCLE_CONTRACTION"))
  
  
  pdf(paste0("./heart-atlas/scPagwas/figures/scPagwas_pathways_", GWAS, ".pdf"), width = 10)
  all_plots <- map(
    c("Cardiomyocytes"),
    ~ plot_pathways_for_celltype(
      df                   = final_df,
      cell_state           = .x,
      preselect_condition  = "ctrl",
      fc_conditions        = c("ctrl","dcm"),
      n_ctrl               = 200,
      n_up                 = 30,
      n_down               = 30,
      effect_measure       = "log2FC"   # or "log2FC", "r"
    )
  )
  names(all_plots) <- c("Cardiomyocytes")
  
  for (ct in names(all_plots)) {
    cat("\n\n=== Cell type:", ct, "===\n")
    print(all_plots[[ct]]$ctrl_up)
    print(all_plots[[ct]]$ctrl_down)
    print(all_plots[[ct]]$ctrl_facet)
    print(all_plots[[ct]]$dot_up)
    print(all_plots[[ct]]$dot_down)
    print(all_plots[[ct]]$dot_facet)
  }
  dev.off()
}


GWAS = "GCST90435254"
rds_ctrl <- data.frame(
  rds_path = c(list.files(path = "./heart-atlas/scPagwas/outs/reactome", pattern = paste0("^0_Control.*\\", GWAS, ".rds$"), full.names = TRUE)),
  condition = c("ctrl"),
  stringsAsFactors = FALSE
)

rds_hcm <- data.frame(
  rds_path = c(list.files(path = "./heart-atlas/scPagwas/outs/reactome", pattern = paste0("^HCM.*\\", GWAS, ".rds$"), full.names = TRUE)),
  condition = c("hcm"),
  stringsAsFactors = FALSE
)

rds_info <- rbind(rds_hcm, rds_ctrl)

final_df <- prepare_pagwas_dotplot_data_from_df(rds_info)

print("HCM")

pdf(paste0("./heart-atlas/scPagwas/figures/scPagwas_pathways_HCM_", GWAS, ".pdf"), width = 10)
    all_plots <- map(
      c("Cardiomyocytes"),
      ~ plot_pathways_for_celltype(
        df                   = final_df,
        cell_state           = .x,
        preselect_condition  = "ctrl",
        fc_conditions        = c("ctrl","hcm"),
        n_ctrl               = 200,
        n_up                 = 30,
        n_down               = 30,
        effect_measure       = "log2FC"   # or "log2FC", "r"
      )
    )
    names(all_plots) <- c("Cardiomyocytes")
    
    for (ct in names(all_plots)) {
      cat("\n\n=== Cell type:", ct, "===\n")
      print(all_plots[[ct]]$ctrl_up)
      print(all_plots[[ct]]$ctrl_down)
      print(all_plots[[ct]]$ctrl_facet)
      print(all_plots[[ct]]$dot_up)
      print(all_plots[[ct]]$dot_down)
      print(all_plots[[ct]]$dot_facet)
    }
    dev.off()
