# MSc Dissertation: Alternative Splicing Analysis

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
- `GO_analysis/GO_analysis_BP_updated.py` - Gene Ontology enrichment analysis script (heatmaps)

### 🛡️ **ISGs/**
Comprehensive analysis of alternative splicing in interferon-stimulated genes:

**Data:** `ISG_Data/`
- Species-specific ISG lists (human, mouse, fruit bat)
- Core vertebrate-mammalian ISGs
- Antiviral gene classifications
- Random control gene sets

**Scripts:**
- `fisher_tests_ISGs.py` - Statistical analysis of ISG splicing patterns
- `Mouse_ISG_extraction.R` - Species-specific ISG data extraction

**Species Comparisons:** `Species_Comparisons/`
- `Extract_Data_ISGs.R` - Multi-species data extraction pipeline
- `Data_analysis_ISGs.R` - Comparative statistical analysis
- `Box_plots_ISGs.R` - Visualisation of species-specific patterns
- `Chi_squared_ISGs.R` - Chi-squared tests for distribution differences
- `Splice_graph_ISGs.R` - Splicing pattern visualisations

**Visualisation:**
- `forest_plots/Forest_plot_ISGs.R` - Forest plots for ISG enrichment results

### 🔬 **CRISPR_screens/**
Analysis of CRISPR-validated antiviral and proviral genes:

**Data Files:**
- Antiviral and proviral gene lists from CRISPR screens
- Control datasets for comparative analysis

**Scripts:**
- `fisher_tests_CRISPR.py` - Statistical analysis of CRISPR-validated genes
- `CRISPR_GO.py` - Gene Ontology analysis for functional validation

**Threshold Analysis:** `threshold_sweep/`
- `create_lists.py` - Generate gene lists across splicing thresholds
- `sweep_plots.py` - Visualisation of threshold-dependent effects

**Visualisation:**
- `Forest_plot_CRISPR.R` - Forest plots for CRISPR screen results

### 🧭 **GWAS/**
Genome-wide association study region analysis:

**Data:** `GWAS_Data/`
- Disease-specific GWAS regions (immune, cancer, brain)
- Intersection and window-based gene overlaps
- `Immune_GWAS/` - Detailed immune-related GWAS data

**Scripts:**
- `bootstrap_fisher_test_INTERSECT.py` - Direct overlap analysis
- `bootstrap_fisher_test_WINDOW.py` - Window-based proximity analysis
- `fisher_tests_GWAS_IMMUNE.py` - Immune-specific GWAS analysis

**Visualisation:**
- `Forest_plot_GWAS_Bedtools.R` - Comparative forest plots for overlap methods
- `Forest_plot_GWAS_GO_annotated.R` - GO-annotated GWAS visualisations

### 📊 **isoform_distribution/**
Comparative analysis of isoform type distributions:

**Data:** `Isoform_comparison_data/`
- High-splicing gene sets from various analyses
- Categorised gene lists for distribution analysis

**Scripts:**
- `stacked_barplot.py` - Comprehensive isoform distribution analysis with Chi-squared testing

### 🔄 **comparing_databases/**
Comparative analysis between different annotation databases:

**Scripts:**
- `splice_analysis_by_database.py` - Cross-database splice pattern comparison
- `splice_analysis_KENDALL_TAU.py` - Correlation analysis using Kendall's tau

### 🎯 **GO_analysis/**
Gene Ontology functional enrichment analysis:

**Scripts:**
- `GO_analysis.py` - Comprehensive GO term enrichment analysis

### 🛠️ **other_scripts/**
Utility scripts for data processing and annotation:

**Scripts:**
- `annotate_genes.py` - Gene annotation with isoform data
- `extract_top_percent.py` - Extract top percentile genes by splicing complexity

## Methodology

### Statistical Analysis
- **Fisher's Exact Tests** with bootstrap methodology for robust statistical inference
- **Multiple testing correction** using Benjamini-Hochberg FDR control
- **Effect size estimation** using odds ratios with 95% confidence intervals
- **Threshold sensitivity analysis** to assess robustness of findings

### Data Processing
- **Ensembl annotation** as primary reference for isoform quantification
- **Cross-species orthology mapping** for comparative analysis
- **Quality control filtering** to ensure robust statistical inference
- **Standardised gene naming** (uppercase) for consistent matching

### Visualisation
- **Forest plots** for effect size visualisation across categories
- **Threshold sweep plots** for sensitivity analysis
- **Stacked bar plots** for isoform distribution comparison
- **Box plots** for species-specific pattern comparison

## Key Findings

This analysis provides comprehensive evidence for:
1. **Enrichment of alternative splicing** in positively selected genes across mammalian evolution
2. **Species-specific splicing patterns** in interferon-stimulated genes
3. **Functional validation** of CRISPR screen results through splicing analysis
4. **Disease-relevant splicing patterns** in GWAS-associated regions
5. **Biotype-specific differences** in isoform complexity distributions

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

## Citation

If you use this code or methodology, please cite:
*[Dissertation title and author details to be added upon completion]*

## Contact

**Student ID:** B270551  
**Institution:** University of Edinburgh  
**Programme:** MSc Bioinformatics

---

*This repository represents the complete computational component of an MSc dissertation investigating alternative splicing patterns across evolutionary, immunological, and disease contexts.*
