## ISG vs Random Gene Splicing Analysis Across Species
## Compares isoform numbers between ISGs and random gene sets for multiple species
## Generates summary statistics and boxplot visualizations
## Usage: Ensure ISG and random CSV files are in working directory, run script
## Created by B270551
##---------------------------------------------------------##

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(dplyr)
library(ggplot2)

analyze_file <- function(species, file, gene_set) {
  if (!file.exists(file)) {
    cat("Error: File", file, "not found for", species, " (", gene_set, ")\n")
    return(NULL)
  }
  
  data <- tryCatch({
    read.csv(file, stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("Error reading file", file, "for", species, " (", gene_set, "):", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(data)) return(NULL)
  
  required_cols <- c("ensembl_gene_id", "external_gene_name", "total_isoforms", 
                     "protein_coding", "non_protein_coding")
  if (!all(required_cols %in% colnames(data))) {
    cat("Error: Missing required columns in", file, "for", species, " (", gene_set, ")\n")
    return(NULL)
  }
  
  total_genes <- nrow(data)
  if (total_genes == 0) {
    cat("Error: No genes found in", file, "for", species, " (", gene_set, ")\n")
    return(NULL)
  }
  
  perc_spliced_protein <- sum(data$protein_coding > 1) / total_genes * 100
  perc_spliced_total <- sum(data$total_isoforms > 1) / total_genes * 100
  
  list(
    summary = data.frame(
      species = species,
      gene_set = gene_set,
      total_genes = total_genes,
      perc_spliced_protein_coding = perc_spliced_protein,
      perc_spliced_total_isoforms = perc_spliced_total
    ),
    protein_coding = data.frame(
      species = species,
      gene_set = gene_set,
      value = data$protein_coding
    ),
    non_protein_coding = data.frame(
      species = species,
      gene_set = gene_set,
      value = data$non_protein_coding
    )
  )
}

species_list <- list(
  list(species = "chicken", isg_file = "chicken_isg_isoforms.csv", random_file = "chicken_random_isoforms.csv"),
  list(species = "cow", isg_file = "cow_isg_isoforms.csv", random_file = "cow_random_isoforms.csv"),
  list(species = "dog", isg_file = "dog_isg_isoforms.csv", random_file = "dog_random_isoforms.csv"),
  list(species = "fruit_bat", isg_file = "fruit_bat_isg_isoforms.csv", random_file = "fruit_bat_random_isoforms.csv"),
  list(species = "horse", isg_file = "horse_isg_isoforms.csv", random_file = "horse_random_isoforms.csv"),
  list(species = "human", isg_file = "human_isg_isoforms.csv", random_file = "human_random_isoforms.csv"),
  list(species = "microbat", isg_file = "microbat_isg_isoforms.csv", random_file = "microbat_random_isoforms.csv"),
  list(species = "mouse", isg_file = "mouse_isg_isoforms.csv", random_file = "mouse_random_isoforms.csv"),
  list(species = "pig", isg_file = "pig_isg_isoforms.csv", random_file = "pig_random_isoforms.csv"),
  list(species = "rat", isg_file = "rat_isg_isoforms.csv", random_file = "rat_random_isoforms.csv"),
  list(species = "sheep", isg_file = "sheep_isg_isoforms.csv", random_file = "sheep_random_isoforms.csv")
)

results <- lapply(species_list, function(s) {
  isg_result <- analyze_file(s$species, s$isg_file, "ISG")
  random_result <- analyze_file(s$species, s$random_file, "Random")
  
  list(isg = isg_result, random = random_result)
})

summary_data <- bind_rows(
  lapply(results, function(r) bind_rows(r$isg$summary, r$random$summary))
)
protein_coding_data <- bind_rows(
  lapply(results, function(r) bind_rows(r$isg$protein_coding, r$random$protein_coding))
)
non_protein_coding_data <- bind_rows(
  lapply(results, function(r) bind_rows(r$isg$non_protein_coding, r$random$non_protein_coding))
)

write.csv(summary_data, "splicing_summary_isg_and_random.csv", row.names = FALSE)
cat("Splicing summary saved to splicing_summary_isg_and_random.csv\n")

write.csv(protein_coding_data[protein_coding_data$gene_set == "ISG", ], 
          "protein_coding_isoforms_isg.csv", row.names = FALSE)
write.csv(protein_coding_data[protein_coding_data$gene_set == "Random", ], 
          "protein_coding_isoforms_random.csv", row.names = FALSE)
write.csv(non_protein_coding_data[non_protein_coding_data$gene_set == "ISG", ], 
          "non_protein_coding_isoforms_isg.csv", row.names = FALSE)
write.csv(non_protein_coding_data[non_protein_coding_data$gene_set == "Random", ], 
          "non_protein_coding_isoforms_random.csv", row.names = FALSE)
cat("Box plot data saved to protein_coding_isoforms_isg.csv, protein_coding_isoforms_random.csv, ",
    "non_protein_coding_isoforms_isg.csv, and non_protein_coding_isoforms_random.csv\n")

p1 <- ggplot(protein_coding_data, aes(x = species, y = value, fill = gene_set)) +
  geom_boxplot(position = position_dodge(width = 0.45), width = 0.4, notch = TRUE, outlier.shape = NA) +
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
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5)
  ) +
  labs(
    title = "Protein-Coding Isoforms per Gene by Species",
    x = "Species",
    y = "Number of Protein-Coding Isoforms",
    fill = "Gene Set"
  ) +
  scale_fill_brewer(palette = "Set2") +
  coord_cartesian(ylim = c(0, 15))

ggsave("protein_coding_boxplot_isg_and_random.png", p1, width = 10, height = 6)
cat("Box plots saved to protein_coding_boxplot_isg_and_random.png")