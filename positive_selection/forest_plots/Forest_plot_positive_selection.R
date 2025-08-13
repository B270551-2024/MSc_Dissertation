## Forest Plot for Positive Selection Enrichment Analysis
## Compares odds ratios between human and mouse datasets across clades
## Visualises enrichment of splicing in positively selected genes with confidence intervals
## Usage: Ensure fisher_test_results.tsv and fisher_test_results_MOUSE.tsv (generated from positive_selection main script) are in working directory
## Created by B270551
##---------------------------------------------------------##

library(ggplot2)
library(dplyr)
library(readr)

human_results <- read_tsv("fisher_test_results.tsv") %>%
  mutate(source = "Human")
mouse_results <- read_tsv("fisher_test_results_MOUSE.tsv") %>%
  mutate(source = "Mouse")

combined_results <- bind_rows(human_results, mouse_results)

combined_results$clade <- factor(combined_results$clade, levels = unique(human_results$clade))

plot <- ggplot(combined_results, aes(x = clade, y = odds_ratio, color = source)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  coord_flip() +
  scale_color_manual(values = c("Human" = "darkblue", "Mouse" = "orange")) +
  labs(
    title = "Forest Plot: Enrichment of Splicing in Positively Selected Genes (Top 20%)",
    x = "Clade",
    y = "Odds Ratio (95% CI)",
    color = "Isoform Database"
  ) +
  theme_minimal(base_size = 14)

ggsave("forest_plot_enrichment_POSITIVE_HUMAN_MOUSE.png", plot = plot, width = 10, height = 6, dpi = 300)