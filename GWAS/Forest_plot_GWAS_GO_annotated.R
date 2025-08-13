## Immune GWAS Categories Forest Plot Analysis
## Visualises enrichment of protein-coding genes in immune-related GWAS categories
## Compares antiviral, proviral, and other immune genes with controls
## Usage: Ensure fisher_tests_gwas_IMMUNE.tsv is in working directory
## Created by B270551
##---------------------------------------------------------##

library(ggplot2)
library(readr)
library(dplyr)

df <- read_tsv("fisher_tests_gwas_IMMUNE.tsv")

gene_counts <- tibble(
  Category = c("immune_all", "immune_antiviral", "immune_proviral", "immune_non_anti/proviral", "random_genes"),
  N = c(3670, 617, 423, 2630, 500)
)

df <- df %>%
  left_join(gene_counts, by = "Category") %>%
  mutate(
    group = case_when(
      Category %in% c("immune_non_anti/proviral", "random_genes") ~ "Control",
      Category == "immune_all" ~ "All Immune",
      Category == "immune_proviral" ~ "Proviral",
      Category == "immune_antiviral" ~ "Antiviral",
      TRUE ~ "Other"
    ),
    Category_label = case_when(
      Category == "immune_all" ~ paste0("Immune Genes (All) (N = ", N, ")"),
      Category == "immune_antiviral" ~ paste0("Immune Genes (Antiviral) (N = ", N, ")"),
      Category == "immune_proviral" ~ paste0("Immune Genes (Proviral) (N = ", N, ")"),
      Category == "immune_non_anti/proviral" ~ paste0("Immune Genes (Non-Pro/Antiviral) (N = ", N, ")"),
      Category == "random_genes" ~ paste0("Random Genes (N = ", N, ")"),
      TRUE ~ Category
    )
  )

df$Category_label <- factor(df$Category_label, levels = rev(c(
  paste0("Immune Genes (All) (N = ", df$N[df$Category == "immune_all"][1], ")"),
  paste0("Immune Genes (Antiviral) (N = ", df$N[df$Category == "immune_antiviral"][1], ")"),
  paste0("Immune Genes (Proviral) (N = ", df$N[df$Category == "immune_proviral"][1], ")"),
  paste0("Immune Genes (Non-Pro/Antiviral) (N = ", df$N[df$Category == "immune_non_anti/proviral"][1], ")"),
  paste0("Random Genes (N = ", df$N[df$Category == "random_genes"][1], ")")
)))

df$group <- factor(df$group, levels = c("Control", "All Immune", "Proviral", "Antiviral"))

group_colours <- c(
  "Control" = "gray40",
  "All Immune" = "#1f77b4",
  "Proviral" = "#e6ab02",
  "Antiviral" = "#1b9e77"
)

plot <- ggplot(df, aes(x = Category_label, y = Odds_Ratio, colour = group)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "gray50") +
  scale_colour_manual(values = group_colours, name = "Group") +
  scale_y_continuous(limits = c(0, 20)) +
  coord_flip() +
  labs(
    title = "Enrichment of Protein-Coding Genes in Immune GWAS Categories",
    x = "Gene Set",
    y = "Odds Ratio (95% CI)"
  ) +
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

ggsave("forest_plot_immune_GWAS.png", plot = plot, width = 11, height = 6, dpi = 300)