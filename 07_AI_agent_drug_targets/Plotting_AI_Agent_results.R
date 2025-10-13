# ---- (optional) install once ----
library(jsonlite)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)
library(forcats)
library(ggforce)
library(ggwordcloud)
library(stringi)

# Input files in JSON format from the AI agent
file1 <- "pubmed_pathway_analysis_cardiomyocytes_DCM_up_agonist_final.json"
file2 <- "pubmed_pathway_analysis_cardiomyocytes_DCM_down_agonist_final.json"

`%||%` <- function(x, y) if (!is.null(x)) x else y
clamp01 <- function(x) pmin(pmax(x, 0), 1)

# -------------------------------
# Read one JSON → list(paths)
# -------------------------------
read_paths <- function(file) {
  dat <- fromJSON(file, simplifyVector = FALSE)
  paths <- dat$pathways %||% list()
  if (is.data.frame(paths)) paths <- purrr::transpose(paths)
  paths
}

# -------------------------------
# Build raw metrics from one source (no scaling yet)
# -------------------------------
metrics_from_paths <- function(paths) {
  map_dfr(paths, function(p) {
    pi <- p$pathway_info %||% list()
    la <- p$literature_analysis %||% list()
    pq <- la$pubmed_query_results %||% list()
    to <- p$total_occurrences %||% list()
    
    tibble(
      pathway       = (pi$name %||% pi$id) %||% "Unknown pathway",
      NES_score     = as.numeric(pi$nes_score %||% NA_real_),
      citations     = as.numeric(pq$total_citations %||% 0),
      disease_count = as.numeric(to$disease %||% 0),
      hf_count      = as.numeric(to$heart_failure %||% 0),
      comp_count    = as.numeric(to$compounds %||% 0),
      gen_count     = as.numeric(to$genetics %||% 0)
    )
  })
}

# -------------------------------
# Build keyword terms (for word cloud) from one source
# -------------------------------
terms_from_paths <- function(paths) {
  empty_df <- tibble(term = character(), freq = double(),
                     category = character(), pathway = character())
  
  extract_terms <- function(p) {
    pi  <- p$pathway_info %||% list()
    la  <- p$literature_analysis %||% list()
    path <- (pi$name %||% pi$id) %||% "Unknown pathway"
    
    # Genes
    gms <- la$gene_mentions %||% list()
    g <- if (length(gms) > 0) {
      tibble(
        term     = map_chr(gms, ~ .x$gene %||% NA_character_),
        freq     = map_dbl(gms, ~ as.numeric((.x$total_mentions %||% .x$total_mentioned) %||% NA_real_)),
        category = "Genes", pathway = path
      )
    } else empty_df[0, ]
    
    # Processes
    pks <- la$process_keywords %||% list()
    proc <- if (length(pks) > 0) {
      tibble(
        term     = map_chr(pks, ~ .x$process %||% NA_character_),
        freq     = map_dbl(pks, ~ as.numeric(.x$count %||% NA_real_)),
        category = "Processes", pathway = path
      )
    } else empty_df[0, ]
    
    # Diseases
    dks <- la$disease_keywords %||% list()
    dis <- if (length(dks) > 0) {
      tibble(
        term     = map_chr(dks, ~ .x$disease %||% NA_character_),
        freq     = map_dbl(dks, ~ as.numeric(.x$count %||% NA_real_)),
        category = "Diseases", pathway = path
      )
    } else empty_df[0, ]
    
    # Compounds
    cks <- la$compound_keywords %||% list()
    cpd <- if (length(cks) > 0) {
      tibble(
        term     = map_chr(cks, ~ .x$compound %||% NA_character_),
        freq     = map_dbl(cks, ~ as.numeric(.x$count %||% NA_real_)),
        category = "Compounds", pathway = path
      )
    } else empty_df[0, ]
    
    bind_rows(g, proc, dis, cpd)
  }
  
  map_dfr(paths, extract_terms)
}

# ===============================
# 1) Build clean data frames from BOTH files
# ===============================
paths1 <- read_paths(file1)
paths2 <- read_paths(file2)

metrics_raw <- bind_rows(
  metrics_from_paths(paths1),
  metrics_from_paths(paths2)
)

# lock pathway levels for consistent colors/order
all_pathways <- sort(unique(metrics_raw$pathway))
metrics_raw <- metrics_raw %>% mutate(pathway = factor(pathway, levels = all_pathways))

# >>> Aggregate across files <<<
# Sum counts/citations; mean NES across appearances (ignore NA)
metrics_agg <- metrics_raw %>%
  group_by(pathway) %>%
  summarise(
    NES_score     = if (all(is.na(NES_score))) NA_real_ else mean(NES_score, na.rm = TRUE),
    citations     = sum(citations, na.rm = TRUE),
    disease_count = sum(disease_count, na.rm = TRUE),
    hf_count      = sum(hf_count, na.rm = TRUE),
    comp_count    = sum(comp_count, na.rm = TRUE),
    gen_count     = sum(gen_count, na.rm = TRUE),
    .groups = "drop"
  )

# Compute per-publication rates from aggregated counts
metrics_rates <- metrics_agg %>%
  mutate(
    disease = if_else(citations > 0, disease_count/1,0), # / citations, 0),
    heart_failure = if_else(citations > 0, hf_count/1,0),#      / citations, 0),
    compounds = if_else(citations > 0, comp_count/1,0),#    / citations, 0),
    genetics = if_else(citations > 0, gen_count/1,0)#     / citations, 0)
  )

# Long form in intended axis order
var_levels <- c("NES_score","citations","disease","heart_failure","compounds","genetics")
metrics_long <- metrics_rates %>%
  transmute(
    pathway,
    NES_score, citations, disease, heart_failure, compounds, genetics
  ) %>%
  pivot_longer(
    cols = all_of(var_levels),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    pathway = factor(pathway, levels = all_pathways),
    variable = factor(variable, levels = var_levels)
  )

# Ensure every pathway has all axes
metrics_long <- metrics_long %>%
  group_by(pathway) %>%
  complete(variable = factor(var_levels, levels = var_levels),
           fill = list(value = NA_real_)) %>%
  ungroup()

# ===============================
# 2) Normalization (global, once)
# ===============================
scale_tbl <- metrics_long %>%
  filter(variable %in% c("citations","disease","heart_failure","compounds","genetics")) %>%
  group_by(variable) %>%
  summarise(minv = min(value, na.rm = TRUE),
            maxv = max(value, na.rm = TRUE),
            .groups = "drop")

metrics_scaled <- metrics_long %>%
  left_join(scale_tbl, by = "variable") %>%
  mutate(
    scaled = case_when(
      variable == "NES_score" ~ {
        s <- (value + 3) / 6               # -3 -> 0, 0 -> 0.5, +3 -> 1
        ifelse(is.na(s), 0.5, clamp01(s))  # neutral for NA
      },
      TRUE ~ {
        v <- ifelse(is.na(value), 0, value)
        ifelse((maxv - minv) > 0, (v - minv) / (maxv - minv), 0.5)
      }
    ),
    scaled = clamp01(scaled)
  ) %>%
  select(pathway, variable, scaled)

# ===============================
# 3) Convert to Cartesian for straight-edged radar
# ===============================
n_axes <- length(var_levels)
axis_layout <- tibble(
  variable = factor(var_levels, levels = var_levels),
  idx      = seq_len(n_axes),
  theta    = pi/2 - 2*pi*(idx-1)/n_axes,
  x_unit   = cos(theta),
  y_unit   = sin(theta)
)

df_xy <- metrics_scaled %>%
  left_join(axis_layout, by = "variable") %>%
  arrange(pathway, idx) %>%
  mutate(x = scaled * x_unit,
         y = scaled * y_unit)

# Helpers
ring_vals <- c(0.25, 0.50, 0.75, 1.00)
rings  <- tibble(r = ring_vals)
spokes <- axis_layout %>% transmute(x0 = 0, y0 = 0, x1 = 1.02*x_unit, y1 = 1.02*y_unit)
label_radius <- 1.10
label_df <- axis_layout %>%
  mutate(
    lx   = label_radius * x_unit,
    ly   = label_radius * y_unit,
    hjust = ifelse(x_unit >= 0, 0, 1),
    vjust = ifelse(abs(y_unit) < 1e-6, 0.5, ifelse(y_unit > 0, 0, 1))
  )

# ===============================
# 4) PLOT — Single combined radar (no faceting)
# ===============================
p_radar <-
  ggplot() +
  ggforce::geom_circle(data = rings, aes(x0 = 0, y0 = 0, r = r),
                       linetype = "dashed", color = "grey80") +
  geom_segment(data = spokes, aes(x = x0, y = y0, xend = x1, yend = y1),
               color = "grey80") +
  geom_text(data = tibble(y = ring_vals, lab = percent(ring_vals)),
            aes(x = 0, y = y, label = lab),
            size = 4, color = "grey40", vjust = -0.4) +
  geom_polygon(
    data = df_xy,
    aes(x = x, y = y, group = pathway, fill = pathway, color = pathway),
    alpha = 0.18, linewidth = 0.6
  ) +
  geom_point(
    data = df_xy,
    aes(x = x, y = y, color = pathway),
    size = 1.8
  ) +
  geom_text(
    data = label_df,
    aes(x = lx, y = ly, label = variable, hjust = hjust, vjust = vjust),
    size = 3.5
  ) +
  scale_color_brewer(palette = "Dark2", direction = -1, drop = FALSE) +
  scale_fill_brewer(palette = "Dark2", direction = -1, drop = FALSE) +
  coord_equal(xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35), expand = FALSE) +
  labs(
    x = NULL, y = NULL, color = "Pathway", fill = "Pathway"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.title = element_blank(),
    legend.position = "right"
  )

# ===============================
# 5) PLOT — Word cloud faceted by pathway (color = category)
# ===============================
terms_raw <- bind_rows(
  terms_from_paths(paths1),
  terms_from_paths(paths2)
) %>%
  filter(!is.na(term), !is.na(freq), freq > 0) %>%
  group_by(pathway, category, term) %>%     # aggregate across files
  summarise(freq = sum(freq), .groups = "drop") %>%
  mutate(
    pathway  = factor(pathway, levels = all_pathways),
    category = factor(category, levels = c("Genes","Processes","Diseases","Compounds")),
    freq_plot = sqrt(freq)
  )

# Keep the top terms per pathway to reduce clutter
top_n_per_pathway <- 40
terms_top <- terms_raw %>%
  group_by(pathway) %>%
  slice_max(order_by = freq, n = top_n_per_pathway, with_ties = FALSE) %>%
  ungroup()

wrap_spaces <- function(x, width = 16) {
  vapply(x, function(s) paste(stri_wrap(s, width = width, whitespace_only = TRUE), collapse = "\n"),
         FUN.VALUE = character(1))
}
terms_top$pathway_lbl <- wrap_spaces(terms_top$pathway)

set.seed(123)
p_cloud <-
  ggplot(terms_top, aes(label = term, size = freq_plot, color = category)) +
  ggwordcloud::geom_text_wordcloud_area(rm_outside = TRUE, show.legend = TRUE) +
  scale_size_area(max_size = 14, guide = "none") +
  scale_color_brewer(palette = "RdBu", type = "div", direction = -1) +
  facet_wrap(~ pathway_lbl, scales = "free", nrow = 2) +
  labs(color = "Category") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines"),
    legend.position = "bottom"
  )


# ===============================
# Render plots for visualization
# ===============================
print(p_radar)
print(p_cloud)

