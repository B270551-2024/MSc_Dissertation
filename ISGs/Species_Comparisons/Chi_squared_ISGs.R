# ISG splicing chi-squared tests
library(dplyr)
library(magrittr)

df <- read.csv("splicing_summary_isg_and_random.csv")

# Function to perform chi-squared tests
perform_chi_sq_tests <- function(df, perc_col, analysis_type) {
  # Calculate spliced and non-spliced counts
  df$spliced <- round(df$total_genes * df[[perc_col]] / 100)
  df$non_spliced <- df$total_genes - df$spliced
  
  # Species-specific tests
  species_results <- df %>%
    group_by(species) %>%
    summarise(
      isg_total = total_genes[gene_set == "ISG"],
      isg_spliced = spliced[gene_set == "ISG"],
      random_total = total_genes[gene_set == "Random"],
      random_spliced = spliced[gene_set == "Random"],
      isg_perc = .data[[perc_col]][gene_set == "ISG"],
      random_perc = .data[[perc_col]][gene_set == "Random"],
      difference = isg_perc - random_perc,
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      p_value = {
        if(isg_total > 0 & random_total > 0 & (isg_spliced + random_spliced) > 0) {
          matrix_data <- matrix(c(isg_spliced, isg_total - isg_spliced,
                                  random_spliced, random_total - random_spliced), 
                                nrow = 2, byrow = TRUE)
          test_result <- tryCatch(chisq.test(matrix_data), error = function(e) NULL)
          if(is.null(test_result)) NA else test_result$p.value
        } else NA
      },
      analysis_type = analysis_type
    ) %>%
    ungroup()
  
  # Overall test
  total_isg <- sum(df$spliced[df$gene_set == "ISG"])
  total_isg_genes <- sum(df$total_genes[df$gene_set == "ISG"])
  total_random <- sum(df$spliced[df$gene_set == "Random"])
  total_random_genes <- sum(df$total_genes[df$gene_set == "Random"])
  
  overall_matrix <- matrix(c(total_isg, total_isg_genes - total_isg,
                             total_random, total_random_genes - total_random), 
                           nrow = 2, byrow = TRUE)
  overall_test <- chisq.test(overall_matrix)
  
  overall_result <- data.frame(
    species = "OVERALL",
    isg_total = total_isg_genes,
    isg_spliced = total_isg,
    random_total = total_random_genes,
    random_spliced = total_random,
    isg_perc = 100 * total_isg / total_isg_genes,
    random_perc = 100 * total_random / total_random_genes,
    difference = 100 * total_isg / total_isg_genes - 100 * total_random / total_random_genes,
    p_value = overall_test$p.value,
    analysis_type = analysis_type
  )
  
  return(rbind(species_results, overall_result))
}

# Run tests for both protein coding and total isoforms
protein_coding_results <- perform_chi_sq_tests(df, "perc_spliced_protein_coding", "protein_coding")
total_isoforms_results <- perform_chi_sq_tests(df, "perc_spliced_total_isoforms", "total_isoforms")

# Combine results and select only required columns
all_results <- rbind(protein_coding_results, total_isoforms_results)
all_results <- all_results[, c("species", "isg_perc", "random_perc", "difference", "p_value", "analysis_type")]

# Save to file
write.csv(all_results, "isg_splicing_chi_squared_results.csv", row.names = FALSE)

# Print summary
cat("PROTEIN CODING GENES:\n")
cat("Overall:", sprintf("ISG %.1f%% vs Random %.1f%%, p = %.2e\n", 
                        protein_coding_results$isg_perc[protein_coding_results$species == "OVERALL"],
                        protein_coding_results$random_perc[protein_coding_results$species == "OVERALL"],
                        protein_coding_results$p_value[protein_coding_results$species == "OVERALL"]))

cat("\nTOTAL ISOFORMS:\n")
cat("Overall:", sprintf("ISG %.1f%% vs Random %.1f%%, p = %.2e\n", 
                        total_isoforms_results$isg_perc[total_isoforms_results$species == "OVERALL"],
                        total_isoforms_results$random_perc[total_isoforms_results$species == "OVERALL"],
                        total_isoforms_results$p_value[total_isoforms_results$species == "OVERALL"]))

cat("\nResults saved to: isg_splicing_chi_squared_results.csv")