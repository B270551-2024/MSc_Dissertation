## GWAS Categories Bootstrap Fisher Test Analysis
## Tests enrichment of splicing in GWAS-associated genes across different disease categories
## Analyses both protein-coding and non-protein-coding biotypes with bootstrap downsampling
## Usage: Ensure all GWAS intersection files are in working directory
## Created by B270551
##---------------------------------------------------------##

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2

background_file = "ensembl_isoforms.tsv"
gwas_categories = {
    "all": "all_gwas_regions_intersect_UNIQUE.txt",
    "immune": "immune_gwas_regions_intersect_UNIQUE.txt",
    "cancer": "cancer_gwas_regions_intersect_UNIQUE.txt",
    "brain": "brain_gwas_regions_intersect_UNIQUE.txt",
    "non_gwas": "non_gwas_genes.bed",
    "random_genes": "random_background_genes.bed",
    "positive_gwas_all": "positive_gwas_all_intersect.txt"
}
splicing_threshold_percent = 0.1
bootstrap_iterations = 100
random_state = 42

bg_df = pd.read_csv(background_file, sep="\t")
bg_df["gene_name"] = bg_df["gene_name"].astype(str).str.strip().str.upper()

pseudogene_biotypes = ["pseudogene", "processed_pseudogene", "unprocessed_pseudogene", 
                       "transcribed_pseudogene", "transcribed_processed_pseudogene", 
                       "transcribed_unprocessed_pseudogene"]
biotype_groups = {
    "protein_coding": ["protein_coding"],
    "non_protein_coding": [b for b in bg_df["biotype"].unique() if b != "protein_coding" and b not in pseudogene_biotypes]
}

def bootstrap_fisher_test(gwas_genes, background_genes, label, biotype_group, threshold, n_iter=100, random_state=42):
    np.random.seed(random_state)
    
    bg_filtered = background_genes[background_genes["biotype"].isin(biotype_groups[biotype_group])].copy()
    if bg_filtered.empty:
        print(f"No genes found for biotype group: {biotype_group}")
        return None
    
    bg_filtered = bg_filtered.drop_duplicates(subset="gene_name").copy()
    bg_filtered["is_target"] = bg_filtered["gene_name"].isin(gwas_genes["gene_name"])
    
    q_low = bg_filtered["num_isoforms"].quantile(threshold)
    q_high = bg_filtered["num_isoforms"].quantile(1 - threshold)
    
    bg_filtered["splicing_class"] = bg_filtered["num_isoforms"].apply(
        lambda x: "high" if x >= q_high else ("low" if x <= q_low else "mid")
    )
    
    extremes = bg_filtered[bg_filtered["splicing_class"] != "mid"].copy()
    
    high_group = extremes[extremes["splicing_class"] == "high"]
    low_group = extremes[extremes["splicing_class"] == "low"]
    
    high_target_orig = high_group["is_target"].sum()
    low_target_orig = low_group["is_target"].sum()
    high_bg_orig = len(high_group) - high_target_orig
    low_bg_orig = len(low_group) - low_target_orig
    
    if high_target_orig < 12 or low_target_orig < 3:
        table = [[high_target_orig + 1, low_target_orig + 1], [high_bg_orig + 1, low_bg_orig + 1]]
        odds, p = fisher_exact(table)
        return {
            "clade": f"{label}_{biotype_group}",
            "mean_odds_ratio": odds,
            "se_odds_ratio": 0,
            "mean_ci_low": odds,
            "mean_ci_high": odds,
            "significant_p_count": 1 if p < 0.05 else 0,
            "total_iterations": 1
        }
    
    min_group_size = min(len(high_group), len(low_group))
    odds_ratios, p_values = [], []
    
    for _ in range(n_iter):
        high_sample = high_group.sample(n=min_group_size, replace=False)
        low_sample = low_group.sample(n=min_group_size, replace=False)
        
        high_target = high_sample["is_target"].sum()
        low_target = low_sample["is_target"].sum()
        high_bg = len(high_sample) - high_target
        low_bg = len(low_sample) - low_target
        
        table = [[high_target + 1, low_target + 1], [high_bg + 1, low_bg + 1]]
        
        odds, p = fisher_exact(table)
        odds_ratios.append(odds)
        p_values.append(p)
    
    odds_ratios, p_values = np.array(odds_ratios), np.array(p_values)
    
    return {
        "clade": f"{label}_{biotype_group}",
        "mean_odds_ratio": np.mean(odds_ratios),
        "se_odds_ratio": np.std(odds_ratios, ddof=1) / np.sqrt(len(odds_ratios)),
        "mean_ci_low": np.percentile(odds_ratios, 2.5),
        "mean_ci_high": np.percentile(odds_ratios, 97.5),
        "significant_p_count": np.sum(p_values < 0.05),
        "total_iterations": len(odds_ratios)
    }

results = []

for biotype_group in biotype_groups:
    for category, filepath in gwas_categories.items():
        if filepath.endswith(".bed"):
            gwas_genes = pd.read_csv(filepath, sep="\t", header=None, usecols=[3], names=["gene_name"])
        else:
            gwas_genes = pd.read_csv(filepath, header=None, names=["gene_name"])
        gwas_genes["gene_name"] = gwas_genes["gene_name"].astype(str).str.strip().str.upper()
        
        res = bootstrap_fisher_test(gwas_genes, bg_df[["gene_name", "num_isoforms", "biotype"]], 
                                   category, biotype_group, splicing_threshold_percent, 
                                   bootstrap_iterations, random_state)
        if res is not None:
            results.append(res)

results_df = pd.DataFrame(results)
cols = ["clade", "mean_odds_ratio", "se_odds_ratio", "mean_ci_low", "mean_ci_high", "significant_p_count", "total_iterations"]
results_df = results_df[cols]

results_df[results_df["clade"].str.contains("protein_coding")].to_csv(
    "fisher_tests_gwas_protein_coding_INTERSECT_BOOTSTRAP.tsv", sep="\t", index=False
)
results_df[results_df["clade"].str.contains("non_protein_coding")].to_csv(
    "fisher_tests_gwas_non_protein_coding_INTERSECT_BOOTSTRAP.tsv", sep="\t", index=False
)

print("Analysis complete. Results saved to 'fisher_tests_gwas_protein_coding_INTERSECT_BOOTSTRAP.tsv' and 'fisher_tests_gwas_non_protein_coding_INTERSECT_BOOTSTRAP.tsv'")
print(results_df.to_string(index=False))
