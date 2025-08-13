## ISG Forest Plot Analysis with Bootstrap Results
## Visualises enrichment of splicing in different ISG gene sets across species
## Creates grouped forest plot comparing ISGs, antiviral genes, interferome, and controls
## Usage: Ensure fisher_tests_ISGs_ALL.tsv is in working directory
## Created by B270551
##---------------------------------------------------------##

library(ggplot2)
library(readr)
library(dplyr)

df <- read_tsv("fisher_tests_ISGs_ALL.tsv")

df <- df %>%
  mutate(group = case_when(
    clade %in% c("Random_Non_ISGs", "Positive_Genes_All") ~ "Negative Control",
    clade == "Mammalian_Interferome" ~ "Interferome",
    clade %in% c("Human_ISGs", "Mouse_ISGs", "FruitBat_ISGs", 
                 "Core_Vert_Mamm_ISGs", "Positively_Selected_ISGs") ~ "ISGs",
    clade %in% c("Antiviral_Genes_All", "Core_Antiviral_ISGs") ~ "Antiviral ISGs",
    TRUE ~ "Other"
  )) %>%
  mutate(display_label = recode(clade,
                                "Random_Non_ISGs" = "Random Non-ISGs",
                                "Positive_Genes_All" = "Positively Selected Genes (All)",
                                "Mammalian_Interferome" = "Mammalian Interferome",
                                "Human_ISGs" = "Human ISGs",
                                "Mouse_ISGs" = "Mouse ISGs",
                                "FruitBat_ISGs" = "Fruit Bat ISGs",
                                "Positively_Selected_ISGs" = "Positively Selected Genes (ISGs)",
                                "Core_Vert_Mamm_ISGs" = "Core ISGs (vert-mamm)",
                                "Antiviral_Genes_All" = "Antiviral Genes (All)",
                                "Core_Antiviral_ISGs" = "Core Antiviral ISGs"
  )) %>%
  mutate(label = paste0(display_label, " (", significant_p_count, "/", total_iterations, ")"))

clade_order <- c(
  "Random_Non_ISGs",
  "Positive_Genes_All",
  "Mammalian_Interferome",
  "Human_ISGs",
  "Mouse_ISGs",
  "FruitBat_ISGs",
  "Positively_Selected_ISGs",
  "Core_Vert_Mamm_ISGs",
  "Antiviral_Genes_All",
  "Core_Antiviral_ISGs"
)

df$clade <- factor(df$clade, levels = rev(clade_order))
label_order <- df$label[match(rev(clade_order), df$clade)]

df$group <- factor(df$group, levels = c(
  "Negative Control",
  "Interferome",
  "ISGs",
  "Antiviral ISGs"
))

group_colours <- c(
  "Negative Control" = "gray40",
  "Interferome" = "#1f77b4",
  "ISGs" = "#ff7f0e",
  "Antiviral ISGs" = "#2ca02c"
)

plot <- ggplot(df, aes(x = clade, y = mean_odds_ratio, colour = group)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_ci_low, ymax = mean_ci_high), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "gray50") +
  scale_colour_manual(values = group_colours, name = "Group") +
  scale_x_discrete(labels = label_order) +
  scale_y_continuous(limits = c(0, 20)) +
  coord_flip() +
  labs(
    title = "Forest Plot: Enrichment of Splicing in ISG Sets",
    x = "Gene Set",
    y = "Odds Ratio (95% CI)"
  ) +
  theme_minimal(base_size = 14)

ggsave("forest_plot_ISGs.png", plot = plot, width = 10, height = 6, dpi = 300)