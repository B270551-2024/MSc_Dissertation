## Splicing Percentage Visualization for ISG vs Random Genes
## Creates dot plot comparing percentage of spliced genes across species
## Visualizes both protein-coding and total isoform percentages
## Usage: Ensure splicing_summary_isg_and_random.csv (from data_analysis_isgs.R) is in working directory, run script
## Created by B270551
##---------------------------------------------------------##

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
library(dplyr)
library(ggplot2)
library(tidyr)

data <- read.csv("splicing_summary_isg_and_random.csv", stringsAsFactors = FALSE)

data_long <- data %>%
  pivot_longer(
    cols = c(perc_spliced_protein_coding, perc_spliced_total_isoforms),
    names_to = "metric_type",
    values_to = "percentage"
  ) %>%
  mutate(
    metric = case_when(
      metric_type == "perc_spliced_protein_coding" & gene_set == "ISG" ~ "Protein_Coding_ISG",
      metric_type == "perc_spliced_protein_coding" & gene_set == "Random" ~ "Protein_Coding_Random",
      metric_type == "perc_spliced_total_isoforms" & gene_set == "ISG" ~ "Total_Isoforms_ISG",
      metric_type == "perc_spliced_total_isoforms" & gene_set == "Random" ~ "Total_Isoforms_Random"
    ),
    metric_group = case_when(
      metric_type == "perc_spliced_protein_coding" ~ "Protein_Coding",
      metric_type == "perc_spliced_total_isoforms" ~ "Total_Isoforms"
    )
  )

p <- ggplot(data_long, aes(x = species, y = percentage)) +
  geom_point(aes(color = metric, shape = metric), size = 4, position = position_dodge(width = 0.4)) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5)
  ) +
  labs(
    title = "Percentage of Spliced Genes by Species (Protein Coding and Total Isoforms)",
    x = "Species",
    y = "Percentage of Spliced Genes (%)"
  ) +
  scale_color_manual(
    values = c(
      "Protein_Coding_ISG" = "#1B9E77",  
      "Protein_Coding_Random" = "#7570B3",  
      "Total_Isoforms_ISG" = "#66C2A5",  
      "Total_Isoforms_Random" = "#B2ABD2"  
    ),
    labels = c(
      "Protein_Coding_ISG" = "Protein Coding Isoforms (ISG)",
      "Protein_Coding_Random" = "Protein Coding Isoforms (Random)",
      "Total_Isoforms_ISG" = "Total Isoforms (ISG)",
      "Total_Isoforms_Random" = "Total Isoforms (Random)"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Protein_Coding_ISG" = 16,  
      "Protein_Coding_Random" = 16,  
      "Total_Isoforms_ISG" = 17,  
      "Total_Isoforms_Random" = 17  
    ),
    labels = c(
      "Protein_Coding_ISG" = "Protein Coding Isoforms (ISG)",
      "Protein_Coding_Random" = "Protein Coding Isoforms (Random)",
      "Total_Isoforms_ISG" = "Total Isoforms (ISG)",
      "Total_Isoforms_Random" = "Total Isoforms (Random)"
    )
  ) +
  coord_cartesian(ylim = c(0, 100))

ggsave("splicing_percentage_plot.png", p, width = 12, height = 5, dpi = 300)
cat("Plot saved to splicing_percentage_plot.png\n")