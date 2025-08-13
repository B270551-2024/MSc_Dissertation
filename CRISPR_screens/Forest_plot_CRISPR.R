## CRISPR-validated Gene Forest Plot Analysis
## Visualises enrichment of splicing in antiviral/proviral genes with bootstrap results
## Creates grouped forest plot showing odds ratios and confidence intervals
## Usage: Ensure fisher_tests_CRISPR_bootstrap.tsv is in working directory
## Created by B270551
##---------------------------------------------------------##

library(ggplot2)
library(readr)
library(dplyr)

df <- read_tsv("fisher_tests_CRISPR_bootstrap.tsv")

df <- df %>%
  mutate(group = case_when(
    grepl("Antiviral", clade) ~ "Antiviral",
    grepl("Proviral", clade) ~ "Proviral",
    TRUE ~ "Control"
  ))

df <- df %>%
  mutate(display_label = case_when(
    clade == "Antiviral_Genes" ~ "Antiviral Genes",
    clade == "Proviral_Genes" ~ "Proviral Genes",
    clade == "Antiviral_Genes_Top100_MAIC" ~ "Antiviral Genes (Top 100)",
    clade == "Proviral_Genes_Top100_MAIC" ~ "Proviral Genes (Top 100)",
    clade == "Antiviral_Genes_Positive_Score" ~ "Antiviral Genes (Positively Selected)",
    clade == "Proviral_Genes_Positive_Score" ~ "Proviral Genes (Positively Selected)",
    clade == "Negative_Control" ~ "Random Genes",
    TRUE ~ clade
  ))

df <- df %>%
  mutate(label = paste0(display_label, " (", significant_p_count, "/", total_iterations, ")"))

df$label <- factor(df$label, levels = rev(df$label))
df$group <- factor(df$group, levels = c("Antiviral", "Proviral", "Control"))

group_colours <- c("Antiviral" = "#1b9e77", "Proviral" = "#e6ab02", "Control" = "gray40")

plot <- ggplot(df, aes(x = label, y = mean_odds_ratio, colour = group)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_ci_low, ymax = mean_ci_high), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "gray50") +
  scale_colour_manual(values = group_colours, name = "Gene Set Type") +
  coord_flip() +
  labs(
    title = "Forest Plot: Enrichment of Splicing in CRISPR-validated Genes",
    x = "Gene Set",
    y = "Odds Ratio (95% CI)"
  ) +
  theme_minimal(base_size = 14)

ggsave("forest_plot_bootstrapped_enrichment_grouped.png", plot = plot, width = 11, height = 6, dpi = 300)