## GWAS Categories Forest Plot Analysis
## Compares enrichment of protein-coding genes across GWAS categories using two overlap methods
## Visualises odds ratios with confidence intervals for different disease-associated regions
## Usage: Ensure fisher_tests_gwas_WINDOW_BOOTSTRAP.tsv and fisher_tests_gwas_INTERSECT_BOOTSTRAP.tsv are in working directory
## Created by B270551
##---------------------------------------------------------##

library(ggplot2)
library(dplyr)
library(readr)

results_window <- read_tsv("fisher_tests_gwas_protein_coding_WINDOW_BOOTSTRAP.tsv") %>%
  mutate(Source = "BEDTools Window")
results_bedtools <- read_tsv("fisher_tests_gwas_protein_coding_INTERSECT_BOOTSTRAP.tsv") %>%
  mutate(Source = "BEDTools Intersect")

combined_results <- bind_rows(results_window, results_bedtools)

combined_results <- combined_results %>%
  mutate(
    Category = case_when(
      clade == "all_protein_coding" ~ "GWAS-associated regions (All)",
      clade == "positive_gwas_all_protein_coding" ~ "GWAS-associated regions (All, positively selected)",
      clade == "brain_protein_coding" ~ "GWAS-associated regions (Brain)",
      clade == "cancer_protein_coding" ~ "GWAS-associated regions (Cancer)",
      clade == "immune_protein_coding" ~ "GWAS-associated regions (Immune)",
      clade == "non_gwas_protein_coding" ~ "Non GWAS associated regions",
      clade == "random_genes_protein_coding" ~ "Random Genes"
    ),
    Category = factor(Category, levels = rev(c(
      "GWAS-associated regions (All)",
      "GWAS-associated regions (All, positively selected)",
      "GWAS-associated regions (Brain)",
      "GWAS-associated regions (Cancer)",
      "GWAS-associated regions (Immune)",
      "Non GWAS associated regions",
      "Random Genes"
    )))
  )

plot <- ggplot(combined_results, aes(x = Category, y = mean_odds_ratio, colour = Source)) +
  geom_point(data = combined_results %>% filter(!Category %in% c("Non GWAS associated regions", "Random Genes")),
             position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(data = combined_results %>% filter(!Category %in% c("Non GWAS associated regions", "Random Genes")),
                aes(ymin = mean_ci_low, ymax = mean_ci_high),
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_point(data = combined_results %>% filter(Category %in% c("Non GWAS associated regions", "Random Genes")),
             aes(x = Category, y = mean_odds_ratio),
             colour = "gray50", position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(data = combined_results %>% filter(Category %in% c("Non GWAS associated regions", "Random Genes")),
                aes(x = Category, y = mean_odds_ratio, ymin = mean_ci_low, ymax = mean_ci_high),
                colour = "gray50", position = position_dodge(width = 0.5), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "gray50") +
  scale_colour_manual(
    values = c(
      "BEDTools Window" = "#2E4057",
      "BEDTools Intersect" = "#8A2BE2"
    ),
    name = "Overlap Method"
  ) +
  labs(
    title = "Enrichment of Protein-Coding Genes in GWAS Categories",
    x = "GWAS Category",
    y = "Odds Ratio (95% CI)"
  ) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_line(colour = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(colour = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.box = "horizontal",
    legend.margin = margin(t = 0, b = 5)
  )

ggsave("forest_plot_GWAS_categories.png", plot, width = 11, height = 6, dpi = 300)