## Protein-Coding Isoform Boxplot Analysis with Statistical Testing
## Creates boxplots comparing ISG vs random genes with Wilcoxon tests
## Shows distribution of protein-coding isoforms per gene across species
## Usage: Ensure protein_coding_isoforms_isg.csv and protein_coding_isoforms_random.csv (from data_extraction and non_ISGs.R) are in working directory
## Created by B270551
##---------------------------------------------------------##

library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)

isg_data <- read_csv("protein_coding_isoforms_isg.csv")
random_data <- read_csv("protein_coding_isoforms_random.csv")

isg_data$gene_set <- "ISG"
random_data$gene_set <- "Random"

combined_data <- bind_rows(isg_data, random_data)

combined_data$species <- as.factor(combined_data$species)
combined_data$gene_set <- factor(combined_data$gene_set, levels = c("ISG", "Random"))

p <- ggplot(combined_data, aes(x = species, y = value, fill = gene_set)) +
  geom_boxplot(
    position = position_dodge(width = 0.6),
    width = 0.5,
    outlier.shape = NA
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    position = position_dodge(width = 0.6),
    width = 0.4,
    color = "grey50"
  ) +
  stat_compare_means(
    aes(group = gene_set),
    method = "wilcox.test",
    label = "p.signif",
    label.y = 14.5,
    size = 5,
    bracket.size = 0.6,
    hide.ns = TRUE
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5)
  ) +
  labs(
    title = "Protein-Coding Isoforms per Gene by Species",
    x = "Species",
    y = "Number of Protein-Coding Isoforms",
    fill = "Gene Set"
  ) +
  scale_fill_manual(
    values = c(
      "ISG" = "#66C2A5",
      "Random" = "#B2ABD2"
    ),
    labels = c("ISG" = "ISG", "Random" = "Random")
  ) +
  coord_cartesian(ylim = c(0, 15.5), expand = FALSE)

ggsave("protein_coding_isoforms_boxplot.png", p, width = 12, height = 5, dpi = 300)
cat("Plot saved to protein_coding_isoforms_boxplot.png\n")
print(p)