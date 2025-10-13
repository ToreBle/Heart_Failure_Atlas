library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

Meta_Heart <- read.csv("Heart_atlas_meta_data.csv", row.names = 1)

Meta_Heart$Condition[Meta_Heart$Condition=="myocardial infarction"] <- "MI"
Meta_Heart$Condition[Meta_Heart$Condition=="Long-term ischemic"] <- "MI"

# Define common themes and color palettes
base_theme <- theme(
  panel.background = element_blank(),
  axis.text = element_text(size = 16),
  axis.title = element_text(size = 16),
  legend.position = "top",
  legend.text = element_text(size = 16),
  legend.title = element_blank()
)

Sex_colors <- c("#FDBF6F", "#FF7F00")

# Boxplot function
make_boxplot <- function(data, x, y, fill, fill_palette = NULL, fill_colors = NULL, ylab = "", xlab = "") {
  p <- ggplot(data, aes_string(x = x, y = y, fill = fill)) +
    geom_boxplot(position = "dodge") +
    labs(y = ylab, x = xlab) +
    base_theme
  if (!is.null(fill_palette)) p <- p + scale_fill_brewer(palette = fill_palette)
  if (!is.null(fill_colors)) p <- p + scale_fill_manual(values = fill_colors)
  return(p)
}

# Barplot function
make_barplot <- function(data, x, y, fill, fill_palette = NULL, fill_colors = NULL) {
  p <- ggplot(data, aes_string(x = x, y = y, fill = fill)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(x = "", y = "") +
    base_theme
  if (!is.null(fill_palette)) p <- p + scale_fill_manual(values = brewer.pal(min(9, nlevels(factor(data[[fill]]))), fill_palette))
  if (!is.null(fill_colors)) p <- p + scale_fill_manual(values = fill_colors)
  return(p)
}

make_barplot_dodge <- function(data, x, y, fill, fill_palette = NULL, fill_colors = NULL) {
  p <- ggplot(data, aes_string(x = x, y = y, fill = fill)) +
    geom_bar(position = "dodge", stat = "identity") +
    labs(x = "", y = "") +
    base_theme
  if (!is.null(fill_palette)) p <- p + scale_fill_manual(values = brewer.pal(min(9, nlevels(factor(data[[fill]]))), fill_palette))
  if (!is.null(fill_colors)) p <- p + scale_fill_manual(values = fill_colors)
  return(p)
}

# Age-BMI Plot
Disease <- Meta_Heart %>% count(Age, Sex, BMI, Condition, LVEF, Donor)
Pat <- Disease %>% count(Age, BMI, Condition)
Plot_Age_BMI <- make_boxplot(Pat %>% filter(BMI != "NA"), "Age", "BMI", "Condition", fill_palette = "Dark2", ylab = "BMI", xlab = "Age")

# LVEF Plot
Pat <- Disease %>% filter(Condition %in% c("DCM", "HCM", "Control")) %>% count(LVEF, Sex, Condition)
Plot_LVEF <- make_boxplot(Pat, "Condition", "LVEF", "Sex", fill_colors = Sex_colors, ylab = "LVEF")

ggarrange(Plot_Age_BMI, Plot_LVEF, ncol = 1)

# Donor Composition Barplot by Study
Disease <- Meta_Heart %>% count(Hypertension, LVEF, Diabetes, Sex, Ischemia, Study, Condition, Donor)
Donor_plot <- ggplot(Disease, aes(x = Donor, y = n, fill = Study)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(. ~ Condition, scales = "free", space = "free") +
  labs(x = "", y = "single nuclei") +
  base_theme +
  scale_fill_brewer(palette = "Dark2") +
  theme(strip.text = element_text(size = 14), axis.text.x = element_blank())

print(Donor_plot)

# Sex Distribution Barplot
Pat <- Disease %>% count(Condition, Sex)
Sex_plot <- make_barplot_dodge(Pat, "Condition", "n", "Sex", fill_colors = Sex_colors) +
  labs(y = "Donor")

print(Sex_plot)

# MI Condition Plot
Disease <- Meta_Heart %>% filter(Condition == "MI") %>% count(Ischemia, Condition, Donor, Sex)
Pat <- Disease %>% count(Condition, Ischemia) %>%
  mutate(Ischemia = ifelse(Ischemia == "Ischemia", "Post MI", Ischemia))
MI_plot <- make_barplot(Pat, "n", "Condition", "Ischemia", fill_palette = "Reds")

# DCM Hypertension Plot
Disease <- Meta_Heart %>% filter(Condition == "DCM") %>% count(Hypertension, Condition, Donor, Sex)
Pat <- Disease %>% count(Condition, Hypertension) %>%
  mutate(Hypertension = ifelse(Hypertension == "", "Unknown", Hypertension))
DCM_plot <- make_barplot(Pat, "n", "Condition", "Hypertension", fill_palette = "Blues")

# DCM Diabetes Plot
Disease <- Meta_Heart %>% filter(Condition == "DCM") %>% count(Diabetes, Condition, Donor, Study)
Pat <- Disease %>% count(Condition, Diabetes) %>%
  mutate(Diabetes = ifelse(Diabetes == "", "Unknown", Diabetes))
DCM_diabetes <- make_barplot(Pat, "n", "Condition", "Diabetes", fill_palette = "Blues")

ggarrange(MI_plot, DCM_plot, DCM_diabetes, ncol = 1)

Cell_type_colors <- c("Adipocytes" = '#1f77b4', "Cardiomyocytes" = '#ff7f0e', "Endocardial" = '#279e68', "Endothelial" = '#d62728', 
                      "Epicardium" = '#aa40fc', "Fibroblast" = '#8c564b',"Ischemic cells (MI)" = '#e377c2', "Lymphatic EC" = '#b5bd61', 
                      "Lymphocytes" = '#17becf', "Mast cells" = '#aec7e8', "Myeloid" = '#ffbb78', "Neuronal" = '#98df8a',
                      "Pericytes" = '#ff9896', "VSMC" = '#c5b0d5')


high_cell_types <- Meta_Heart %>%
  count(Annotation_1) %>%
  arrange(n) %>%
  mutate(Annotation_1 = factor(Annotation_1, levels = Annotation_1)) %>%
  ggplot(aes(y = Annotation_1, x = n, fill = Annotation_1)) +
  geom_bar(position = "dodge", stat = "identity")+base_theme+
  scale_fill_manual(values = Cell_type_colors)+theme(legend.position = "", axis.text.x = element_text(size=12, colour = "black"))+
  labs(x="", y="")+xlim(0,610000)
low_cell_types <- Meta_Heart %>%
  count(Annotation_1) %>%
  arrange(n) %>% filter(n<100000)%>%
  mutate(Annotation_1 = factor(Annotation_1, levels = Annotation_1)) %>%
  ggplot(aes(y = Annotation_1, x = n, fill = Annotation_1)) +
  geom_bar(position = "dodge", stat = "identity")+base_theme+
  scale_fill_manual(values = Cell_type_colors)+theme(legend.position = "", axis.text.x = element_text(size=12, colour = "black"))+
  labs(x="", y="")
ggarrange(high_cell_types, low_cell_types, ncol = 2)

Disease <- Meta_Heart %>% count(Donor, Study, Sex, Age, Condition, Death, BMI, LVEF, Hypertension, Diabetes, Genotype, Assay, Technology)
write.csv(Disease, "Supplementary Table 1 - Clinical meta data.csv")
