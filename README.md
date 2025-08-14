# MSc Dissertation: Alternative Splicing Analysis (Scripts + Data)

**Author:** B270551   
**Programme:** MSc Bioinformatics  
**Year:** 2024-2025

## Overview

This repository contains the complete computational analysis pipeline for investigating alternative splicing patterns across different gene sets, including positively selected genes, interferon-stimulated genes (ISGs), CRISPR-validated genes, and GWAS-associated regions. The analysis employs Fisher's exact tests with bootstrap methodology to assess enrichment of highly spliced genes across various biological contexts. It was created as part of the data analysis required for the dissertation element of the course.

## Repository Structure

### 📁 **master_data/**
Core datasets and reference files used throughout the analysis:
- `ensembl_isoforms.tsv` - Primary Ensembl isoform annotation data
- `ensembl_isoforms_v4.tsv` - Enhanced version with additional biotype classifications
- `ensembl_isoforms_mouse_HUMAN_NAME.tsv` - Mouse-specific isoform data with human gene names
- `all_tested_genes.txt` - Complete list of genes included in positive selection analyses
- `morales_supplementary_table6.txt` - Positive selection reference data from Morales et al. study
- `shaw_ISG_data.txt` - ISG reference data from Shaw et al.

### 🧬 **positive_selection/**
Analysis of alternative splicing in positively selected genes across mammalian clades:

**Data:** `Positive_selection_data/`
- Positive and negative gene sets for 14 mammalian clades
- Control datasets for comparative analysis

**Scripts:**
- `fisher_tests_positive_selection.py` - Fisher's exact tests for enrichment analysis
- `downsampled_fisher_tests.py` - Bootstrap analysis with matched sample sizes (to assess significance in humans)
- `threshold_analysis_polygon.py` - Threshold sweep analysis with polygon plots (to check quantiles from 10-50)
- `extract_genes.py` - Gene extraction and processing utilities

**Visualisation:**
- `forest_plots/Forest_plot_positive_selection.R` - Forest plots creation for visualising odds ratios
- `forest_plots/Forest_plot_positive_downsampled.R` - Downsampled analysis forest plots

**Functional Analysis:**
- `GO_analysis/GO_analysis_BP_updated.py` - Gene Ontology enrichment analysis script (heatmap generation)

### 🛡️ **ISGs/**
Comprehensive analysis of alternative splicing in interferon-stimulated genes:

**Data:** `ISG_Data/`
- Species-specific ISG lists (human, mouse, fruit bat)
- Core vertebrate-mammalian ISGs
- Antiviral gene classifications
- Random control gene sets

**Scripts:**
- `fisher_tests_ISGs.py` - Fishers exact test analysis of ISG splicing patterns
- `Mouse_ISG_extraction.R` - Mouse ISG data extraction (as came from an independent data source)

**Species Comparisons:** `Species_Comparisons/`
- `Extract_Data_ISGs.R` - Multi-species data extraction pipeline (uses shaw et al ISG master data to extract ISG lists)
- `Data_analysis_ISGs.R` - Comparative statistical analysis (analyses splicing, % splicing and number of isoforms)
- `Box_plots_ISGs.R` - Visualisation of species-specific patterns (number of isoforms)
- `Chi_squared_ISGs.R` - Chi-squared tests for distribution differences
- `Splice_graph_ISGs.R` - Splicing pattern visualisations (% splicing)

**Visualisation:**
- `forest_plots/Forest_plot_ISGs.R` - Forest plots for ISG Odds Ratios enrichment results

### 🔬 **CRISPR_screens/**
Analysis of CRISPR-validated antiviral and proviral genes:

**Data Files:**
- Antiviral and proviral gene lists from CRISPR screens (not publically available as obtained from a collaborator)

**Scripts:**
- `fisher_tests_CRISPR.py` - Fisher's exact test analysis of CRISPR-validated genes
- `CRISPR_GO.py` - Gene Ontology analysis for functional validation

**Threshold Analysis:** `threshold_sweep/`
- `create_lists.py` - Generate gene lists across MAIC thresholds (top 100-top 2000)
- `sweep_plots.py` - Visualisation of threshold-dependent effects

**Visualisation:**
- `Forest_plot_CRISPR.R` - Forest plots for CRISPR screen splicing enrichment i.e. odds ratios 

### 🧭 **GWAS/**
Genome-wide association study region analysis:

**Data:** `GWAS_Data/`
- Disease-specific GWAS-associated genes (immune, cancer, brain)
- Intersection and window-based gene overlaps
- `Immune_GWAS/` - Detailed immune-related GWAS data (annotated using GO terms as pro or antiviral)

**Scripts:**
- `bootstrap_fisher_test_INTERSECT.py` - Direct overlap analysis
- `bootstrap_fisher_test_WINDOW.py` - Window-based proximity analysis
- `fisher_tests_GWAS_IMMUNE.py` - Immune-specific subsets GWAS analysis

**Visualisation:**
- `Forest_plot_GWAS_Bedtools.R` - Comparative forest plots for overlap methods
- `Forest_plot_GWAS_GO_annotated.R` - GO-annotated GWAS forest plot visualisations

### 📊 **isoform_distribution/**
Comparative analysis of isoform type distributions:

**Data:** `Isoform_comparison_data/`
- High-splicing gene sets from various analyses
- Categorised gene lists for distribution analysis

**Scripts:**
- `stacked_barplot.py` - Comprehensive isoform distribution analysis with Chi-squared testing - for assessing isoform type changes between groups

### 🔄 **comparing_databases/**
Comparative analysis between different annotation databases:

**Scripts:**
- `splice_analysis_by_database.py` - Cross-database splice pattern comparison
- `splice_analysis_KENDALL_TAU.py` - Correlation analysis using Kendall's tau

### 🎯 **GO_analysis/**
Gene Ontology functional enrichment analysis:

**Scripts:**
- `GO_analysis.py` - Comprehensive GO term enrichment analysis (requires outputs generated from PANTHER GO: https://pantherdb.org/)

### 🛠️ **other_scripts/**
Utility scripts for data processing and annotation:

**Scripts:**
- `annotate_genes.py` - Gene annotation with isoform data
- `extract_top_percent.py` - Extract top percentile genes by splicing complexity

## Methodology


### Visualisation
- **Forest plots** for effect size visualisation across categories
- **Threshold sweep plots** for sensitivity analysis
- **Stacked bar plots** for isoform distribution comparison
- **Box plots** for species-specific pattern comparison


## Dependencies

### Python
- pandas, numpy, scipy
- matplotlib, seaborn
- statsmodels
- gprofiler-official

### R
- ggplot2, dplyr
- ggpubr (for statistical comparisons)
- biomaRt (for gene annotation)

## Usage

Each analysis directory contains self-contained scripts that can be run independently. Ensure all data files are present in the appropriate directories before execution. Scripts are designed to read from relative paths within the repository structure.


## Contact

**Student ID:** B270551  

---

*This repository represents the complete computational component of an MSc dissertation investigating alternative splicing patterns across evolutionary, immunological, and disease contexts.*
