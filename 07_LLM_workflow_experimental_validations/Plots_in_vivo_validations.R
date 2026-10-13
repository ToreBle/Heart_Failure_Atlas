library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(ggpubr)

# Read Excel file
df_wide <- read_excel(
  "./Heart atlas/Experimental_Validation/Raw data echo.xlsx",
  sheet = "EF (%) raw data"
)

df <- df_wide %>%
  pivot_longer(
    cols = c(
      `Naive`,
      `Isoproterenol`,
      `Isoproterenol+ZLN005`
    ),
    names_to = "Group",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

group_order <- c(
  "Naive",
  "Isoproterenol",
  "Isoproterenol+ZLN005"
)

df$Group <- factor(df$Group, levels = group_order)

group_colors <- c(
  "Naive" = "gray65",
  "Isoproterenol" = "#D95F5F",
  "Isoproterenol+ZLN005" = "#8B0000"
)

# Check number of observations
df %>%
  count(Group)

# Summary
summary_df <- df %>%
  group_by(Group) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  )

# T-tests
ttest_results <- df %>%
  t_test(
    Value ~ Group,
    comparisons = list(
      c("Naive", "Isoproterenol"),
      c("Naive", "Isoproterenol+ZLN005"),
      c("Isoproterenol", "Isoproterenol+ZLN005")
    ),
    p.adjust.method = "none"
  ) %>%
  mutate(
    p_label = paste0("p = ", signif(p, 3))
  )

# Calculate suitable heights based on the data
y_max <- max(df$Value, na.rm = TRUE)

# Different height for each comparison
ttest_results$y.position <- c(
  y_max * 1.10,
  y_max * 1.30,
  y_max * 1.20
)

print(ttest_results)


p <- ggplot(
  summary_df,
  aes(x = Group, y = mean, fill = Group)
) +
  
  geom_col(
    width = 0.9
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se
    ),
    width = 0.10,
    linewidth = 0.9
  ) +
  
  # Raw individual data points
  geom_jitter(
    data = df,
    aes(x = Group, y = Value),
    inherit.aes = FALSE,
    width = 0.06,
    size = 3
  ) +
  
  # P-values
  stat_pvalue_manual(
    ttest_results,
    label = "p_label",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.7,
    size = 6,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    values = group_colors
  ) +
  
  scale_x_discrete(
    limits = group_order,
    labels = c(
      "Naive" = "Naive",
      "Isoproterenol" = "MI",
      "Isoproterenol+ZLN005" = "MI +\nZLN005"
    ),
    
    # More space between y-axis and first bar
    expand = expansion(add = c(0.55, 0.15))
  ) +
  
  # Y-axis starts at zero
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.08))
  ) +
  
  coord_cartesian(
    ylim = c(0, y_max * 1.45),
    clip = "off"
  ) +
  
  labs(
    x = NULL,
    y = "Ejection fraction [%]"
  ) +
  
  theme_classic(base_size = 20) +
  
  theme(
    legend.position = "none",
    
    axis.title.y = element_text(
      size = 22,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 19,
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 19,
      color = "black",
      angle = 0,
      hjust = 0.5
    ),
    
    axis.line = element_line(
      linewidth = 0.9,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.9,
      color = "black"
    )
  )

p


troponin_wide <- read_excel(
  "./Heart atlas/Raw data troponin and neonatal CM cell death",
  sheet = "Tabelle1"
)

troponin <- troponin_wide %>%
  pivot_longer(
    cols = c(
      `Isoproterenol`,
      `Isoproterenol+ZLN005`
    ),
    names_to = "Group",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))


###
summary_troponin <- troponin %>%
  group_by(Group) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  )


# T-tests
ttest_results <- troponin %>%
  t_test(
    Value ~ Group,
    comparisons = list(
      c("Isoproterenol", "Isoproterenol+ZLN005")
    ),
    p.adjust.method = "none"
  ) %>%
  mutate(
    p_label = paste0("p = ", signif(p, 3))
  )

# Calculate suitable heights based on the data
y_max <- max(troponin$Value, na.rm = TRUE)

# Different height for each comparison
ttest_results$y.position <- c(
  y_max * 1.10
)

print(ttest_results)


plot_troponin <- ggplot(
  summary_troponin,
  aes(x = Group, y = mean, fill = Group)
) +
  
  geom_col(
    width = 0.9
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se
    ),
    width = 0.10,
    linewidth = 0.9
  ) +
  
  # Raw individual data points
  geom_jitter(
    data = troponin,
    aes(x = Group, y = Value),
    inherit.aes = FALSE,
    width = 0.06,
    size = 3
  ) +
  
  # P-values
  stat_pvalue_manual(
    ttest_results,
    label = "p_label",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.7,
    size = 6,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    values = group_colors
  ) +
  
  scale_x_discrete(
    expand = expansion(add = c(0.55, 0.15))
  ) +
  
  coord_cartesian(
    ylim = c(0, y_max * 1.2),
    clip = "off"
  ) +
  
  labs(
    x = NULL,
    y = "Troponin [ng/L]"
  ) +
  
  theme_classic(base_size = 20) +
  
  theme(
    legend.position = "none",
    
    axis.title.y = element_text(
      size = 22,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 19,
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 19,
      color = "black",
      angle = 0,
      hjust = 0.5
    ),
    
    axis.line = element_line(
      linewidth = 0.9,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.9,
      color = "black"
    )
  )
plot_troponin


neonatal_wide <- read_excel(
  "./Heart atlas/Experimental_Validation/Raw data troponin and neonatal CM cell death.xlsx",
  sheet = "Tabelle2"
)

neonatal <- neonatal_wide %>%
  pivot_longer(
    cols = c(
      `Isoproterenol`,
      `Isoproterenol+ZLN005`
    ),
    names_to = "Group",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))



###
summary_neonatal <- neonatal %>%
  group_by(Group) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  )


# T-tests
ttest_results <- neonatal %>%
  t_test(
    Value ~ Group,
    comparisons = list(
      c("Isoproterenol", "Isoproterenol+ZLN005")
    ),
    p.adjust.method = "none"
  ) %>%
  mutate(
    p_label = paste0("p = ", signif(p, 3))
  )

# Calculate suitable heights based on the data
y_max <- max(neonatal$Value, na.rm = TRUE)

# Different height for each comparison
ttest_results$y.position <- c(
  y_max * 1.10
)

print(ttest_results)


plot_neonatal <- ggplot(
  summary_neonatal,
  aes(x = Group, y = mean, fill = Group)
) +
  
  geom_col(
    width = 0.9
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se
    ),
    width = 0.10,
    linewidth = 0.9
  ) +
  
  # Raw individual data points
  geom_jitter(
    data = neonatal,
    aes(x = Group, y = Value),
    inherit.aes = FALSE,
    width = 0.06,
    size = 3
  ) +
  
  # P-values
  stat_pvalue_manual(
    ttest_results,
    label = "p_label",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.7,
    size = 6,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    values = group_colors
  ) +
  
  scale_x_discrete(
    expand = expansion(add = c(0.55, 0.15))
  ) +
  
  coord_cartesian(
    ylim = c(0, y_max * 1.2),
    clip = "off"
  ) +
  
  labs(
    x = NULL,
    y = "Cell death [%]"
  ) +
  
  theme_classic(base_size = 20) +
  
  theme(
    legend.position = "none",
    
    axis.title.y = element_text(
      size = 22,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 19,
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 19,
      color = "black",
      angle = 0,
      hjust = 0.5
    ),
    
    axis.line = element_line(
      linewidth = 0.9,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.9,
      color = "black"
    )
  )
plot_neonatal
