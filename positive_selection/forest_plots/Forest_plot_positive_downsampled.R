# Script that does the downsampling assessment 
# Usage: need fisher_test_downsampled_summary.tsv from downsampling script 

library(ggplot2)
library(dplyr)
library(readr)

# Load the summary dataset
results <- read_tsv("fisher_test_downsampled_summary.tsv") %>%
  mutate(clade = if_else(clade == "Negative Control", "Random Genes", clade))

# Define the desired clade order
clade_order <- c("ALL", "Bat", "Mouse", "Human", "Pholidota", "Perissodactyla", 
                 "Lagomorpha", "Eulipotyphla", "Carnivora", "Afrotheria", 
                 "Artiodactyla", "Chiroptera", "Rodentia", "Primates", "Random Genes")

# Ensure clade is a factor with the specified order
results$clade <- factor(results$clade, levels = clade_order)

# Create labels with significant p-value count
results <- results %>%
  mutate(label = paste0(clade, " (", significant_p_count, "/", total_iterations, ")"))

# Create a color column based on clade
results <- results %>%
  mutate(plot_color = if_else(clade == "Random Genes", "grey50", "darkblue"))

# Forest plot
plot <- ggplot(results, aes(x = reorder(label, match(clade, rev(clade_order))), y = mean_odds_ratio, color = plot_color)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_ci_low, ymax = mean_ci_high), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  coord_flip() +
  scale_color_identity() +
  labs(
    title = "Forest Plot: Enrichment of Splicing in Positively Selected Genes (Top 10%)",
    x = "Clade",
    y = "Mean Odds Ratio (95% CI)"
  ) +
  theme_minimal(base_size = 14)

# Save the plot
ggsave("forest_plot_enrichment_downsampled.png", plot = plot, width = 10, height = 8, dpi = 300)