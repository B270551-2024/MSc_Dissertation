# A script for doing a KEGG analysis 
# Created by Max Falk

packages <- c("gprofiler2", "ggplot2", "dplyr", "readr", "forcats")
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)
if (length(to_install) > 0) install.packages(to_install)

library(gprofiler2)
library(ggplot2)
library(dplyr)
library(readr)
library(forcats)


gene_list <- read_lines("top_10_percent_proviral_genes.tsv") %>% trimws()

background_list <- read_tsv("all_proviral_genes.tsv", col_names = FALSE) %>%
  pull(1) %>%
  trimws() %>%
  unique()

# Runing KEGG enrichment using g:Profiler with custom background
gprof_res <- gost(
  query = gene_list,
  organism = "hsapiens",
  sources = "KEGG",
  correction_method = "fdr",
  user_threshold = 0.05,
  significant = TRUE,
  custom_bg = background_list
)


if (is.null(gprof_res$result) || nrow(gprof_res$result) == 0) {
  stop("No significant KEGG enrichment found.")
}


kegg_df <- gprof_res$result %>%
  filter(term_name != "KEGG root term") %>%
  mutate(
    enrichment = (intersection_size / query_size) / (term_size / effective_domain_size),
    term_name = fct_reorder(term_name, -log10(p_value))
  ) %>%
  arrange(p_value) %>%
  slice_head(n = 20)  # Top 20 terms


#MAIN PLOT - 7C
p3 <- ggplot(kegg_df, aes(x = enrichment, y = term_name, fill = -log10(p_value))) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(
    low = "#abd9e9",
    high = "#d73027",
    name = "-log10(FDR)"
  ) +
  labs(
    title = "KEGG Enrichment Bar Plot\n(proviral spliced vs proviral all)",
    x = "Fold Enrichment",
    y = "Pathway"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

ggsave("KEGG_barplot_enrichment_pvalue_PROVIRAL_SPLICED_VS_PROVIRAL_ALL.png", plot = p3, width = 9, height = 6, dpi = 300)
