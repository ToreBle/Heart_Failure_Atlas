# ============================================================
# Combined pathway PubMed analysis:
# Agonist versus Inhibitor
# ============================================================

library(jsonlite)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggforce)
library(ggwordcloud)
library(stringi)
library(stringr)

`%||%` <- function(x, y) if (!is.null(x)) x else y
clamp01 <- function(x) pmin(pmax(x, 0), 1)


# ============================================================
# 1. Input files
# ============================================================

agonist_file <- paste0(
  "./Heart atlas/Claude_MCP/",
  "pubmed_pathway_analysis_2026-08-05_agonist.json"
)

inhibitor_file <- paste0(
  "./Heart atlas/Claude_MCP/",
  "pubmed_pathway_analysis_2026-08-05_inhibitor.json"
)


# ============================================================
# 2. Read both JSONs
# ============================================================

read_analysis_json <- function(file, context) {
  
  dat <- fromJSON(
    file,
    simplifyVector = FALSE
  )
  
  paths <- dat$pathways %||% list()
  
  if (is.data.frame(paths)) {
    paths <- purrr::transpose(paths)
  }
  
  if (length(paths) == 0) {
    stop(
      paste(
        "No pathways found in",
        file
      )
    )
  }
  
  list(
    paths = paths,
    metadata = dat$analysis_metadata %||% list(),
    context = context
  )
}


agonist_dat <- read_analysis_json(
  agonist_file,
  "Agonist"
)

inhibitor_dat <- read_analysis_json(
  inhibitor_file,
  "Inhibitor"
)


analyses <- list(
  agonist_dat,
  inhibitor_dat
)


# ============================================================
# 3. Extract pathway-level metrics
# ============================================================

extract_metrics <- function(analysis) {
  
  map_dfr(
    analysis$paths,
    
    function(p) {
      
      pi <- p$pathway_info %||% list()
      la <- p$literature_analysis %||% list()
      pq <- la$pubmed_query_results %||% list()
      totals <- p$total_occurrences %||% list()
      
      tibble(
        
        context = analysis$context,
        
        pathway =
          (pi$name %||% pi$id) %||%
          "Unknown pathway",
        
        NES_score =
          as.numeric(
            pi$nes_score %||% NA_real_
          ),
        
        gene_count =
          as.numeric(
            pi$gene_count %||% 0
          ),
        
        citations =
          as.numeric(
            pq$total_citations %||% 0
          ),
        
        abstracts =
          as.numeric(
            pq$abstracts_processed %||% 0
          ),
        
        disease =
          as.numeric(
            totals$disease %||% 0
          ),
        
        heart_failure =
          as.numeric(
            totals$heart_failure %||% 0
          ),
        
        compounds =
          as.numeric(
            totals$compounds %||% 0
          ),
        
        genetics =
          as.numeric(
            totals$genetics %||% 0
          )
      )
    }
  )
}


metrics_raw <- map_dfr(
  analyses,
  extract_metrics
)


# ============================================================
# 4. Factors / ordering
# ============================================================

all_pathways <- unique(
  metrics_raw$pathway
)

metrics_raw <- metrics_raw %>%
  
  mutate(
    
    pathway = factor(
      pathway,
      levels = all_pathways
    ),
    
    context = factor(
      context,
      levels = c(
        "Agonist",
        "Inhibitor"
      )
    )
  )


# ============================================================
# 5. Metrics for radar plot
# ============================================================

var_levels <- c(
  "NES_score",
  "citations",
  "disease",
  "heart_failure",
  "compounds",
  "genetics"
)


metrics_long <- metrics_raw %>%
  
  select(
    context,
    pathway,
    all_of(var_levels)
  ) %>%
  
  pivot_longer(
    
    cols = all_of(var_levels),
    
    names_to = "variable",
    
    values_to = "value"
  ) %>%
  
  mutate(
    
    variable = factor(
      variable,
      levels = var_levels
    )
  )


# ============================================================
# 6. NES normalization
#
# The pathway with the highest ABSOLUTE NES = 100%.
#
# Example:
#
# NES:
#   -1.30
#   -1.82
#   -2.27
#
# becomes:
#    57%
#    80%
#   100%
#
# Sign is intentionally ignored ONLY for radar magnitude.
# Original NES_score is retained in metrics_raw.
# ============================================================

max_abs_nes <- max(
  abs(metrics_raw$NES_score),
  na.rm = TRUE
)


cat(
  "\nMaximum absolute NES:",
  max_abs_nes,
  "\n"
)


# ============================================================
# 7. Scaling for the other radar variables
#
# Keep global scaling across BOTH agonist and inhibitor so the
# two contexts remain directly comparable.
# ============================================================

scale_tbl <- metrics_long %>%
  
  filter(
    variable != "NES_score"
  ) %>%
  
  group_by(
    variable
  ) %>%
  
  summarise(
    
    minv = min(
      value,
      na.rm = TRUE
    ),
    
    maxv = max(
      value,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


metrics_scaled <- metrics_long %>%
  
  left_join(
    scale_tbl,
    by = "variable"
  ) %>%
  
  mutate(
    
    scaled = case_when(
      
      # ------------------------------------------------------
      # NES:
      # absolute strongest NES across all pathways = 100%
      # ------------------------------------------------------
      variable == "NES_score" ~
        
        ifelse(
          is.na(value),
          0,
          abs(value) / max_abs_nes
        ),
      
      # ------------------------------------------------------
      # Other metrics:
      # global min-max scaling across agonist + inhibitor
      # ------------------------------------------------------
      TRUE ~ {
        
        v <- ifelse(
          is.na(value),
          0,
          value
        )
        
        ifelse(
          
          is.finite(maxv - minv) &
            (maxv - minv) > 0,
          
          (v - minv) /
            (maxv - minv),
          
          0.5
        )
      }
    ),
    
    scaled = clamp01(
      scaled
    )
  ) %>%
  
  select(
    context,
    pathway,
    variable,
    value,
    scaled
  )


# ============================================================
# 8. Radar coordinates
# ============================================================

n_axes <- length(
  var_levels
)


axis_layout <- tibble(
  
  variable = factor(
    var_levels,
    levels = var_levels
  ),
  
  idx = seq_len(
    n_axes
  )
) %>%
  
  mutate(
    
    theta =
      pi / 2 -
      2 * pi *
      (idx - 1) /
      n_axes,
    
    x_unit = cos(theta),
    
    y_unit = sin(theta)
  )


df_xy <- metrics_scaled %>%
  
  left_join(
    axis_layout,
    by = "variable"
  ) %>%
  
  arrange(
    pathway,
    context,
    idx
  ) %>%
  
  mutate(
    
    x = scaled * x_unit,
    
    y = scaled * y_unit
  )


# ============================================================
# 9. Radar background
# ============================================================

ring_vals <- c(
  0.25,
  0.50,
  0.75,
  1.00
)


rings <- tibble(
  r = ring_vals
)


spokes <- axis_layout %>%
  
  transmute(
    
    x0 = 0,
    y0 = 0,
    
    x1 = 1.02 * x_unit,
    y1 = 1.02 * y_unit
  )


label_df <- axis_layout %>%
  
  mutate(
    
    label = recode(
      
      as.character(
        variable
      ),
      
      NES_score = "Absolute NES",
      
      heart_failure =
        "Heart failure"
    ),
    
    lx = 1.15 * x_unit,
    
    ly = 1.15 * y_unit,
    
    hjust = ifelse(
      x_unit >= 0,
      0,
      1
    ),
    
    vjust = ifelse(
      
      abs(y_unit) < 1e-6,
      
      0.5,
      
      ifelse(
        y_unit > 0,
        0,
        1
      )
    )
  )


# ============================================================
# Radar plot
# All pathways combined in one radar
# Agonist and Inhibitor shown for each pathway
# ============================================================

p_radar <- ggplot() +
  
  ggforce::geom_circle(
    data = rings,
    aes(
      x0 = 0,
      y0 = 0,
      r = r
    ),
    linetype = "dashed",
    color = "grey75",
    linewidth = 0.6
  ) +
  
  geom_segment(
    data = spokes,
    aes(
      x = x0,
      y = y0,
      xend = x1,
      yend = y1
    ),
    color = "grey75",
    linewidth = 0.6
  ) +
  
  # Ring labels on right side
  geom_text(
    data = tibble(
      x = ring_vals,
      y = 0,
      lab = percent(
        ring_vals,
        accuracy = 1
      )
    ),
    aes(
      x = x,
      y = y,
      label = lab
    ),
    size = 4.5,
    color = "grey35",
    hjust = -0.15,
    vjust = -0.4,
    fontface = "bold"
  ) +
  
  geom_polygon(
    data = df_xy,
    aes(
      x = x,
      y = y,
      group = interaction(
        pathway,
        context
      ),
      color = pathway,
      linetype = context
    ),
    fill = NA,
    linewidth = 1.3
  ) +
  
  geom_point(
    data = df_xy,
    aes(
      x = x,
      y = y,
      color = pathway,
      shape = context
    ),
    size = 3.5,
    stroke = 1
  ) +
  
  geom_text(
    data = label_df,
    aes(
      x = lx,
      y = ly,
      label = label,
      hjust = hjust,
      vjust = vjust
    ),
    size = 5,
    fontface = "bold"
  ) +
  
  scale_color_manual(
    values = c(
      "MITOCHONDRIAL BIOGENESIS"      = "#E68613",
      "MITOCHONDRIAL PROTEIN IMPORT" = "#C23B3B",
      "AUTOPHAGY"                     = "#7A3E9D"
    ),
    drop = FALSE
  ) +
  
  scale_linetype_manual(
    values = c(
      "Agonist" = "solid",
      "Inhibitor" = "dashed"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Agonist" = 16,
      "Inhibitor" = 17
    )
  ) +
  
  coord_equal(
    xlim = c(-1.48, 1.58),
    ylim = c(-1.48, 1.48),
    expand = FALSE
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    color = "Pathway",
    linetype = "Search type",
    shape = "Search type"
  ) +
  
  theme_minimal(
    base_size = 16
  ) +
  
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    
    plot.title = element_text(
      face = "bold",
      size = 20
    ),
    
    plot.subtitle = element_text(
      size = 15,
      margin = margin(b = 12)
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 14
    ),
    
    legend.text = element_text(
      size = 13
    ),
    
    legend.position = "right"
  )


print(p_radar)



# ============================================================
# 11. Extract terms for combined word clouds
# ============================================================

empty_terms <- tibble(
  
  term = character(),
  
  freq = double(),
  
  category = character(),
  
  pathway = character(),
  
  context = character()
)


extract_terms <- function(
    p,
    context
) {
  
  pi <- p$pathway_info %||% list()
  
  la <- p$literature_analysis %||% list()
  
  pathway_name <-
    (pi$name %||% pi$id) %||%
    "Unknown pathway"
  
  
  # ----------------------------------------------------------
  # Genes
  # ----------------------------------------------------------
  
  genes <- la$gene_mentions %||% list()
  
  
  gene_terms <- if (
    length(genes) > 0
  ) {
    
    tibble(
      
      term = map_chr(
        genes,
        ~ .x$gene %||% NA_character_
      ),
      
      freq = map_dbl(
        
        genes,
        
        ~ as.numeric(
          (.x$total_mentions %||%
             .x$total_mentioned) %||%
            0
        )
      ),
      
      category = "Genes",
      
      pathway = pathway_name,
      
      context = context
    )
    
  } else {
    
    empty_terms[0, ]
  }
  
  
  # ----------------------------------------------------------
  # Processes
  # ----------------------------------------------------------
  
  processes <-
    la$process_keywords %||%
    list()
  
  
  process_terms <- if (
    length(processes) > 0
  ) {
    
    tibble(
      
      term = map_chr(
        processes,
        ~ .x$process %||% NA_character_
      ),
      
      freq = map_dbl(
        processes,
        ~ as.numeric(
          .x$count %||% 0
        )
      ),
      
      category = "Processes",
      
      pathway = pathway_name,
      
      context = context
    )
    
  } else {
    
    empty_terms[0, ]
  }
  
  
  # ----------------------------------------------------------
  # Diseases
  # ----------------------------------------------------------
  
  diseases <-
    la$disease_keywords %||%
    list()
  
  
  disease_terms <- if (
    length(diseases) > 0
  ) {
    
    tibble(
      
      term = map_chr(
        diseases,
        ~ .x$disease %||% NA_character_
      ),
      
      freq = map_dbl(
        diseases,
        ~ as.numeric(
          .x$count %||% 0
        )
      ),
      
      category = "Diseases",
      
      pathway = pathway_name,
      
      context = context
    )
    
  } else {
    
    empty_terms[0, ]
  }
  
  
  # ----------------------------------------------------------
  # Compounds
  # ----------------------------------------------------------
  
  compounds <-
    la$compound_keywords %||%
    list()
  
  
  compound_terms <- if (
    length(compounds) > 0
  ) {
    
    tibble(
      
      term = map_chr(
        compounds,
        ~ .x$compound %||% NA_character_
      ),
      
      freq = map_dbl(
        compounds,
        ~ as.numeric(
          .x$count %||% 0
        )
      ),
      
      category = "Compounds",
      
      pathway = pathway_name,
      
      context = context
    )
    
  } else {
    
    empty_terms[0, ]
  }
  
  
  bind_rows(
    gene_terms,
    process_terms,
    disease_terms,
    compound_terms
  )
}

# ============================================================
# 12. Extract terms from BOTH analyses
# ============================================================

terms_agonist <- map_dfr(
  
  agonist_dat$paths,
  
  ~ extract_terms(
    .x,
    "Agonist"
  )
)


terms_inhibitor <- map_dfr(
  
  inhibitor_dat$paths,
  
  ~ extract_terms(
    .x,
    "Inhibitor"
  )
)


terms_raw <- bind_rows(
  
  terms_agonist,
  
  terms_inhibitor
  
) %>%
  
  filter(
    !is.na(term),
    term != "",
    is.finite(freq),
    freq > 0
  ) %>%
  
  group_by(
    pathway,
    context,
    category,
    term
  ) %>%
  
  summarise(
    freq = sum(freq),
    .groups = "drop"
  ) %>%
  
  mutate(
    
    pathway = factor(
      pathway,
      levels = all_pathways
    ),
    
    context = factor(
      context,
      levels = c(
        "Agonist",
        "Inhibitor"
      )
    ),
    
    category = factor(
      
      category,
      
      levels = c(
        "Genes",
        "Processes",
        "Diseases",
        "Compounds"
      )
    ),
    
    # Reduce dominance of highly frequent terms
    freq_plot = sqrt(
      freq
    )
  )

# ============================================================
# 13. Top terms independently for each pathway + context
# ============================================================

TOP_N_PER_PANEL <- 40


terms_top <- terms_raw %>%
  
  group_by(
    pathway,
    context
  ) %>%
  
  slice_max(
    
    order_by = freq,
    
    n = TOP_N_PER_PANEL,
    
    with_ties = FALSE
  ) %>%
  
  ungroup()

# ============================================================
# Top 9 terms per category, pathway, and context
# ============================================================

TOP_N_PER_CATEGORY <- 10

terms_top <- terms_raw %>%
  
  group_by(
    pathway,
    context,
    category
  ) %>%
  
  slice_max(
    order_by = freq,
    n = TOP_N_PER_CATEGORY,
    with_ties = FALSE
  ) %>%
  
  ungroup()

# ============================================================
# 14. Word-cloud plot
#
# Each ROW = one pathway
# Two COLUMNS = Agonist / Inhibitor
# ============================================================

set.seed(
  123
)


p_cloud <- ggplot(
  
  terms_top,
  
  aes(
    label = term,
    size = freq_plot,
    color = category
  )
  
) +
  
  ggwordcloud::geom_text_wordcloud_area(
    
    rm_outside = TRUE,
    
    show.legend = TRUE
  ) +
  
  scale_size_area(
    
    max_size = 12,
    
    guide = "none"
  ) +
  
  scale_color_brewer(
    palette = "RdBu", type = "div", direction = -1
  ) +
  
  # ----------------------------------------------------------
# KEY CHANGE:
#
# pathway = rows
# Agonist/Inhibitor = columns
# ----------------------------------------------------------
facet_grid(
  
  rows = vars(
    pathway
  ),
  
  cols = vars(
    context
  ),
  
  scales = "free"
) +
  
  labs(
    
    title =
      "Frequent genes, processes, diseases and compounds",
    
    subtitle =
      "Agonist and inhibitor literature searches shown side-by-side for each pathway",
    
    color =
      "Category"
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    
    strip.text =
      element_text(
        face = "bold"
      ),
    
    strip.text.y =
      element_text(
        angle = 0
      ),
    
    panel.spacing =
      grid::unit(
        1,
        "lines"
      ),
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    legend.position =
      "bottom"
  )


print(
  p_cloud
)

# ============================================================
# 15. Export
# ============================================================

ggsave(
  "pathway_radar_agonist_vs_inhibitor.png",
  p_radar,
  width = 15,
  height = 6,
  dpi = 300
)


ggsave(
  "pathway_wordcloud_agonist_vs_inhibitor.png",
  p_cloud,
  width = 14,
  height = 12,
  dpi = 300
)
