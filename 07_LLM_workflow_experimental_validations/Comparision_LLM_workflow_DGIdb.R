library(jsonlite)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)
library(readr)

# Compare DGIdb drug-gene interactions against PubMed literature compounds
# ============================================================================
# 1. Input literature JSON files
# ============================================================================

agonist_json <- "./Heart atlas/Claude_MCP/pubmed_pathway_analysis_2026-08-05_agonist.json"

inhibitor_json <- "./Heart atlas/Claude_MCP/pubmed_pathway_analysis_2026-08-05_inhibitor.json"


# ============================================================================
# 2. Helper functions
# ============================================================================

# Drug-name normalization used as the current fallback matching method.
#
# Examples:
#   "N-Acetylcysteine" -> "nacetylcysteine"
#   "N acetyl cysteine" -> "nacetylcysteine"

normalize_drug <- function(x) {
  
  x %>%
    str_to_lower() %>%
    str_trim() %>%
    str_replace_all("[^a-z0-9]", "")
}


# Safely collapse unique non-missing values
collapse_unique <- function(x) {
  
  x <- unique(x[!is.na(x) & x != ""])
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(sort(x), collapse = "; ")
}


# Try to extract CHEMBL identifier from DGIdb concept ID.
# This intentionally returns NA when the DGIdb concept is not represented
# as a ChEMBL identifier.

extract_chembl_id <- function(x) {
  
  x_upper <- toupper(x)
  
  result <- str_extract(
    x_upper,
    "CHEMBL[0-9]+"
  )
  
  result
}


# ============================================================================
# 3. Classify ALL DGIdb interactions
# ============================================================================

dgidb_classified <- drug_gene_pathway %>%
  
  mutate(
    
    interaction_group = case_when(
      
      interaction_type %in% c(
        "activator [ACTIVATING]",
        "modulator"
      ) ~ "Agonist",
      
      interaction_type %in% c(
        "inhibitor [INHIBITORY]",
        "negative modulator [INHIBITORY]"
      ) ~ "Inhibitor",
      
      is.na(interaction_type) |
        interaction_type == "" ~ "Unspecified",
      
      # Keep other DGIdb interaction annotations rather than silently
      # throwing them away.
      TRUE ~ "Other"
    ),
    
    drug_normalized = normalize_drug(drug),
    
    chembl_id = extract_chembl_id(
      drug_concept_id
    )
  )


# Inspect how DGIdb interactions are classified
dgidb_class_counts <- dgidb_classified %>%
  count(
    interaction_group,
    sort = TRUE
  )

print(dgidb_class_counts)

# ============================================================================
# 4. Read literature JSON files
# ============================================================================

agonist_literature <- fromJSON(
  agonist_json,
  simplifyVector = FALSE
)

inhibitor_literature <- fromJSON(
  inhibitor_json,
  simplifyVector = FALSE
)


# ============================================================================
# 5. Extract literature compounds
# ============================================================================

extract_literature_compounds <- function(json_object, category) {
  
  map_dfr(
    json_object$pathways,
    
    function(x) {
      
      compounds <- x$literature_analysis$compound_keywords
      
      if (is.null(compounds) || length(compounds) == 0) {
        return(tibble())
      }
      
      map_dfr(
        compounds,
        
        function(y) {
          
          tibble(
            
            category = category,
            
            pathway = x$pathway_info$name,
            
            literature_drug = y$compound,
            
            drug_normalized = normalize_drug(
              y$compound
            ),
            
            literature_count = as.numeric(
              y$count
            ),
            
            pubchem_cid = if (
              is.null(y$pubchem_cid)
            ) {
              NA_character_
            } else {
              as.character(y$pubchem_cid)
            },
            
            pubchem_validated = if (
              is.null(y$pubchem_validated)
            ) {
              FALSE
            } else {
              as.logical(y$pubchem_validated)
            }
          )
        }
      )
    }
  )
}


literature_all <- bind_rows(
  
  extract_literature_compounds(
    agonist_literature,
    "Agonist"
  ),
  
  extract_literature_compounds(
    inhibitor_literature,
    "Inhibitor"
  )
)


# ============================================================================
# 6. Reduce literature to unique drug + pathway + category
# ============================================================================

literature_unique <- literature_all %>%
  
  group_by(
    category,
    pathway,
    drug_normalized
  ) %>%
  
  summarise(
    
    literature_drug = first(
      literature_drug
    ),
    
    literature_count = sum(
      literature_count,
      na.rm = TRUE
    ),
    
    pubchem_cid = collapse_unique(
      pubchem_cid
    ),
    
    pubchem_validated = any(
      pubchem_validated,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


# ============================================================================
# 7. Reduce DGIdb to unique drug + pathway + interaction group
# ============================================================================

dgidb_unique <- dgidb_classified %>%
  
  # "Other" is retained in dgidb_classified but not used in this specific
  # agonist/inhibitor comparison.
  filter(
    interaction_group %in% c(
      "Agonist",
      "Inhibitor",
      "Unspecified"
    )
  ) %>%
  
  group_by(
    interaction_group,
    pathway,
    drug_normalized
  ) %>%
  
  summarise(
    
    dgidb_drug = first(
      drug
    ),
    
    drug_concept_id = collapse_unique(
      drug_concept_id
    ),
    
    chembl_id = collapse_unique(
      chembl_id
    ),
    
    genes = collapse_unique(
      gene
    ),
    
    n_genes = n_distinct(
      gene
    ),
    
    interaction_types = collapse_unique(
      interaction_type
    ),
    
    max_interaction_score = if (
      all(is.na(interaction_score))
    ) {
      NA_real_
    } else {
      max(
        interaction_score,
        na.rm = TRUE
      )
    },
    
    sources = collapse_unique(
      sources
    ),
    
    pubmed_ids = collapse_unique(
      pubmed_ids
    ),
    
    .groups = "drop"
  )


# ============================================================================
# 8. Build DGIdb datasets for comparison
#
# Directional interactions:
#
#   DGIdb Agonist   -> compare with agonist literature
#   DGIdb Inhibitor -> compare with inhibitor literature
#
# Unspecified interactions:
#
#   Compare with BOTH literature categories because DGIdb does not provide
#   directionality.
# ============================================================================


# --------------------------------------------------------------------------
# 8a. Directional DGIdb records
# --------------------------------------------------------------------------

dgidb_directional <- dgidb_unique %>%
  
  filter(
    interaction_group %in% c(
      "Agonist",
      "Inhibitor"
    )
  ) %>%
  
  mutate(
    
    # Literature category against which it should be tested
    category = interaction_group,
    
    dgidb_direction_status = "Direction specified"
  )


# --------------------------------------------------------------------------
# 8b. DGIdb records with interaction_type = NA
#
# Duplicate each into two comparison contexts:
#
#   Unspecified -> agonist literature
#   Unspecified -> inhibitor literature
#
# This does NOT imply that the compound is an agonist or inhibitor.
# It only asks whether the drug occurs in either literature set.
# --------------------------------------------------------------------------

dgidb_unspecified <- dgidb_unique %>%
  
  filter(
    interaction_group == "Unspecified"
  )


dgidb_unspecified_for_comparison <- bind_rows(
  
  dgidb_unspecified %>%
    mutate(
      category = "Agonist"
    ),
  
  dgidb_unspecified %>%
    mutate(
      category = "Inhibitor"
    )
  
) %>%
  
  mutate(
    dgidb_direction_status = "Direction unknown"
  )


# Combine
dgidb_for_comparison <- bind_rows(
  dgidb_directional,
  dgidb_unspecified_for_comparison
)

# ============================================================================
# 8b. DGIdb interactions with unknown direction
# ============================================================================

dgidb_unspecified <- dgidb_unique %>%
  filter(
    interaction_group == "Unspecified"
  )


# Duplicate temporarily so every unknown-direction DGIdb compound can
# be checked against BOTH literature categories.
dgidb_unspecified_for_comparison <- bind_rows(
  
  dgidb_unspecified %>%
    mutate(category = "Agonist"),
  
  dgidb_unspecified %>%
    mutate(category = "Inhibitor")
  
) %>%
  
  mutate(
    dgidb_direction_status = "Direction unknown"
  )


# ============================================================================
# 8c. Combine directional and unknown-direction DGIdb interactions
# ============================================================================

dgidb_for_comparison <- bind_rows(
  dgidb_directional,
  dgidb_unspecified_for_comparison
)


# ============================================================================
# 9. Compare DGIdb with literature
# ============================================================================

drug_comparison_raw <- full_join(
  
  dgidb_for_comparison,
  
  literature_unique,
  
  by = c(
    "category",
    "pathway",
    "drug_normalized"
  )
  
) %>%
  
  mutate(
    
    source = case_when(
      
      # DGIdb + literature, regardless of whether DGIdb direction is known
      !is.na(dgidb_drug) &
        !is.na(literature_drug) ~
        "Matched",
      
      # Directional DGIdb interaction, absent from literature
      !is.na(dgidb_drug) &
        is.na(literature_drug) &
        dgidb_direction_status == "Direction specified" ~
        "DGIdb only",
      
      # Unknown-direction DGIdb interaction
      !is.na(dgidb_drug) &
        is.na(literature_drug) &
        dgidb_direction_status == "Direction unknown" ~
        "DGIdb candidate",
      
      # Literature / LLM only
      is.na(dgidb_drug) &
        !is.na(literature_drug) ~
        "LLM-based",
      
      TRUE ~ NA_character_
    ),
    
    drug_name = coalesce(
      dgidb_drug,
      literature_drug
    )
  )

# ============================================================================
# 10. Identify ALL unknown-direction DGIdb compounds that matched
#
# If an unknown DGIdb interaction matches either the Agonist OR Inhibitor
# literature, it must NOT subsequently appear in the standalone DGIdb group.
# ============================================================================

unknown_dgidb_matched <- drug_comparison_raw %>%
  
  filter(
    source == "Matched",
    dgidb_direction_status == "Direction unknown"
  ) %>%
  
  distinct(
    pathway,
    drug_normalized
  )


# ============================================================================
# 11. Construct truly unmatched unknown-direction DGIdb set
# ============================================================================

dgidb_unknown_unmatched <- dgidb_unspecified %>%
  
  # Remove anything that matched either literature category
  anti_join(
    unknown_dgidb_matched,
    by = c(
      "pathway",
      "drug_normalized"
    )
  ) %>%
  
  mutate(
    
    category = "DGIdb",
    
    source = "DGIdb",
    
    dgidb_direction_status = "Direction unknown",
    
    literature_drug = NA_character_,
    literature_count = NA_real_,
    pubchem_cid = NA_character_,
    pubchem_validated = NA,
    
    drug_name = dgidb_drug
  )

# ============================================================================
# 12. Final clean comparison
# ============================================================================

drug_comparison_clean <- drug_comparison_raw %>%
  
  # Remove temporary duplicated unknown-direction unmatched records
  filter(
    source != "DGIdb candidate"
  ) %>%
  
  bind_rows(
    dgidb_unknown_unmatched
  )

# ============================================================================
# 13. Inspect matched DGIdb interactions with unknown direction
#
# Although these are now called simply "Matched", it is useful to retain
# the direction-status column so we can see which matches came from NA
# interaction types.
# ============================================================================

matched_unknown_direction <- drug_comparison_clean %>%
  
  filter(
    source == "Matched",
    dgidb_direction_status == "Direction unknown"
  ) %>%
  
  select(
    category,
    pathway,
    drug_name,
    genes,
    drug_concept_id,
    chembl_id,
    literature_count,
    pubchem_cid,
    sources
  ) %>%
  
  arrange(
    category,
    pathway,
    drug_name
  )


cat(
  "\n============================================================\n"
)

cat(
  "MATCHED DRUGS WITH UNKNOWN DGIdb INTERACTION DIRECTION\n"
)

cat(
  "============================================================\n\n"
)

print(
  matched_unknown_direction,
  n = Inf
)


write_csv(
  matched_unknown_direction,
  "matched_unknown_DGIdb_direction.csv"
)


# ============================================================================
# 14. Counts for overall plot
# ============================================================================

source_levels <- c(
  "Matched",
  "DGIdb only",
  "DGIdb",
  "LLM-based"
)


category_levels <- c(
  "Agonist",
  "Inhibitor",
  "DGIdb"
)


comparison_counts <- drug_comparison_clean %>%
  
  count(
    category,
    source,
    name = "n"
  ) %>%
  
  complete(
    
    category = category_levels,
    
    source = source_levels,
    
    fill = list(
      n = 0
    )
  ) %>%
  
  mutate(
    
    category = factor(
      category,
      levels = category_levels
    ),
    
    source = factor(
      source,
      levels = source_levels
    )
  )


print(comparison_counts)


write_csv(
  comparison_counts,
  "DGIdb_vs_LLM_counts.csv"
)


# ============================================================================
# 15. Overall comparison plot
# ============================================================================

# ============================================================================
# 13. Overall plot
#
# facet_grid with space = "free_x" makes panel width proportional to the
# number of bars, so individual bars retain approximately the same width.
# ============================================================================
comparison_counts_plot <- comparison_counts %>%
  filter(n > 0) %>%
  droplevels()

p1 <- ggplot(
  comparison_counts_plot,
  aes(
    x = source,
    y = n,
    fill = source
  )
) +
  
  geom_col(
    width = 0.7,
    show.legend = FALSE
  ) +
  
  geom_text(
    aes(label = n),
    vjust = -0.3,
    size = 3.5
  ) +
  
  facet_grid(
    cols = vars(category),
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_fill_manual(
    values = c(
      "Matched"    = "#8B0000",
      "DGIdb only" = "#D73027",
      "LLM-based"  = "#FC8D8D",
      "DGIdb"      = "#808080"
    )
  ) +
  
  labs(
    x = NULL,
    
    y = "Compound-pathway associations"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    
    plot.title = element_text(
      face = "bold"
    ),
    
    axis.text.x = element_text(
      angle = 35,
      size = 14,
      hjust = 1
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    
    panel.spacing.x = unit(
      1,
      "lines"
    )
  )


print(p1)


ggsave(
  "DGIdb_vs_LLM_overall_faceted.png",
  p1,
  width = 11,
  height = 5,
  dpi = 300
)


ggsave(
  "DGIdb_vs_LLM_overall.png",
  p1,
  width = 10,
  height = 6,
  dpi = 300
)


# ============================================================================
# 16. Counts by pathway
# ============================================================================

pathway_counts <- drug_comparison_clean %>%
  
  count(
    category,
    pathway,
    source,
    name = "n"
  ) %>%
  
  complete(
    
    category = category_levels,
    
    pathway,
    
    source = source_levels,
    
    fill = list(
      n = 0
    )
  ) %>%
  
  mutate(
    
    category = factor(
      category,
      levels = category_levels
    ),
    
    source = factor(
      source,
      levels = source_levels
    )
  )


write_csv(
  pathway_counts,
  "DGIdb_vs_LLM_pathway_counts.csv"
)


# ============================================================================
# 17. Plot by pathway
# ============================================================================

pathway_counts_plot <- pathway_counts %>%
  filter(n > 0) %>%
  droplevels()


p2 <- ggplot(
  pathway_counts_plot,
  aes(
    x = pathway,
    y = n,
    fill = source
  )
) +
  
  geom_col(
    position = position_dodge2(
      width = 0.9,
      preserve = "total",
      padding = 0.05
    ),
    width = 0.85
  ) +
  
  geom_text(
    aes(label = n),
    position = position_dodge2(
      width = 0.9,
      preserve = "total",
      padding = 0.05
    ),
    vjust = -0.25,
    size = 3
  ) +
  
  facet_grid(
    cols = vars(category),
    scales = "free_x",
    space = "free_x"
  ) +
  
  # Remove unused x levels AND outer whitespace
  scale_x_discrete(
    drop = TRUE,
    expand = expansion(
      mult = c(0.02, 0.02)
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Matched"    = "#8B0000",
      "DGIdb only" = "#D73027",
      "LLM-based"  = "#FC8D8D",
      "DGIdb"      = "#808080"
    ),
    drop = TRUE
  ) +
  
  labs(
    x = NULL,
    y = "Compound-pathway associations",
    fill = "Source"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    
    strip.text = element_text(
      face = "bold"
    ),
    
    # Reduce gap BETWEEN facets
    panel.spacing.x = grid::unit(
      0.2,
      "lines"
    ),
    
    legend.position = "top"
  )


print(p2)

ggsave(
  "DGIdb_vs_LLM_by_pathway.png",
  p2,
  width = 13,
  height = 7,
  dpi = 300
)


# ============================================================================
# 18. Literature / LLM recovery by DGIdb
# ============================================================================

llm_recovery <- drug_comparison_clean %>%
  
  filter(
    category %in% c(
      "Agonist",
      "Inhibitor"
    ),
    !is.na(literature_drug)
  ) %>%
  
  mutate(
    
    match_class = case_when(
      
      source == "Matched" ~
        "Matched with DGIdb",
      
      source == "LLM-based" ~
        "LLM-based only",
      
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(
    !is.na(match_class)
  ) %>%
  
  count(
    category,
    match_class,
    name = "n"
  ) %>%
  
  group_by(
    category
  ) %>%
  
  mutate(
    percentage = 100 * n / sum(n)
  ) %>%
  
  ungroup()


print(llm_recovery)


write_csv(
  llm_recovery,
  "LLM_recovery_by_DGIdb.csv"
)