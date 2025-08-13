## ISG Isoform Extraction Script for Multiple Species
## Extracts transcript isoform data for interferon-stimulated genes (ISGs) from Ensembl BioMart
## Usage: Place ISG gene lists as .txt files in working directory, run script
## Outputs: Individual CSV files per species + combined results file
## Created by B270551
##---------------------------------------------------------##

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(biomaRt)
library(dplyr)

process_species <- function(species, gene_file, dataset, host = "https://www.ensembl.org") {
  if (!file.exists(gene_file)) {
    cat("Error: File", gene_file, "not found for", species, "\n")
    return(NULL)
  }
  
  genes <- tryCatch({
    readLines(gene_file)
  }, error = function(e) {
    cat("Error reading file", gene_file, "for", species, ":", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(genes)) return(NULL)
  
  genes <- trimws(genes)
  genes <- genes[genes != ""]
  genes <- unique(genes) 
  
  if (length(genes) == 0 || all(is.na(genes))) {
    cat("Error: No valid gene symbols found in", gene_file, "for", species, "\n")
    return(NULL)
  }
  
  cat("Input gene list size for", species, ":", length(genes), "genes\n")
  
  ensembl <- tryCatch({
    mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset, host = host)
    datasets <- listDatasets(mart)
    if (!dataset %in% datasets$dataset) {
      cat("Error: Dataset", dataset, "not valid for", species, ". Available datasets:\n")
      print(datasets)
      return(NULL)
    }
    mart
  }, error = function(e) {
    cat("Error connecting to BioMart for", species, ":", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(ensembl)) return(NULL)
  
  results <- tryCatch({
    getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                         "ensembl_transcript_id", "transcript_biotype"),
          filters = "external_gene_name",
          values = genes,
          mart = ensembl)
  }, error = function(e) {
    cat("Error querying BioMart for", species, ":", e$message, "\n")
    return(data.frame())
  })
  
  if (nrow(results) == 0) {
    cat("No results found for", species, "\n")
    return(data.frame())
  }
  
  raw_output_file <- paste0(species, "_raw_biomart_results.csv")
  write.csv(results, raw_output_file, row.names = FALSE)
  cat("Raw BioMart results saved to", raw_output_file, "\n")
  
  unique_biotypes <- unique(results$transcript_biotype)
  cat("Unique transcript biotypes for", species, ":", paste(unique_biotypes, collapse = ", "), "\n")
  
  summary <- results %>%
    group_by(ensembl_gene_id, external_gene_name) %>%
    summarise(
      total_isoforms = n(),
      protein_coding = sum(transcript_biotype == "protein_coding"),
      non_protein_coding = sum(transcript_biotype != "protein_coding"),
      .groups = "drop"
    ) %>%
    filter(external_gene_name %in% genes)
  
  if (nrow(summary) > length(genes) * 1.5) {
    cat("Warning: Output gene count (", nrow(summary), ") significantly exceeds input gene count (", 
        length(genes), ") for", species, "\n")
  }
  
  output_file <- paste0(species, "_isg_isoforms.csv")
  write.csv(summary, output_file, row.names = FALSE)
  cat("Processed", species, "- Results saved to", output_file, "with", nrow(summary), "genes\n")
  
  return(summary)
}

species_list <- list(
  list(species = "chicken", gene_file = "Chicken_ISGs.txt", dataset = "ggallus_gene_ensembl", host = "https://www.ensembl.org"),
  list(species = "cow", gene_file = "Cow_ISGs.txt", dataset = "btaurus_gene_ensembl", host = "https://www.ensembl.org"),
  list(species = "dog", gene_file = "Dog_ISGs.txt", dataset = "clfamiliaris_gene_ensembl", host = "https://www.ensembl.org"),
  list(species = "fruit_bat", gene_file = "Fruit_bat_ISGs.txt", dataset = "pvampyrus_gene_ensembl", host = "https://www.ensembl.org"),
  list(species = "horse", gene_file = "Horse_ISGs.txt", dataset = "ecaballus_gene_ensembl", host = "https://www.ensembl.org"),
  list(species = "human", gene_file = "Human_ISGs.txt", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org"),
  list(species = "microbat", gene_file = "Microbat_ISGs.txt", dataset = "mlucifugus_gene_ensembl", host = "https://www.ensembl.org"),
  list(species = "mouse", gene_file = "Mouse_ISGs.txt", dataset = "mmusculus_gene_ensembl", host = "https://www.ensembl.org"),
  list(species = "pig", gene_file = "Pig_ISGs.txt", dataset = "sscrofa_gene_ensembl", host = "https://www.ensembl.org"),
  list(species = "rat", gene_file = "Rat_ISGs.txt", dataset = "rnorvegicus_gene_ensembl", host = "https://www.ensembl.org"),
  list(species = "sheep", gene_file = "Sheep_ISGs.txt", dataset = "oaries_gene_ensembl", host = "https://www.ensembl.org")
)

results <- lapply(species_list, function(s) {
  process_species(s$species, s$gene_file, s$dataset, s$host)
})

combined_results <- bind_rows(results, .id = "species")
write.csv(combined_results, "all_species_isg_isoforms.csv", row.names = FALSE)
cat("Combined results saved to all_species_isg_isoforms.csv\n")