## Immune GWAS Fisher Test Analysis
## Tests enrichment of splicing in immune-related GWAS gene categories
## Compares antiviral, proviral, and other immune genes with controls using Fisher's exact test
## Usage: Ensure all immune GWAS BED files are in working directory
## Created by B270551
##---------------------------------------------------------##

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2

df = pd.read_csv("ensembl_isoforms.tsv", sep="\t")
df = df[df['biotype'] == 'protein_coding']

top_cutoff = np.percentile(df["num_isoforms"], 90)
bottom_cutoff = np.percentile(df["num_isoforms"], 10)

df["splicing_class"] = np.where(df["num_isoforms"] >= top_cutoff, "high",
                                np.where(df["num_isoforms"] <= bottom_cutoff, "low", "mid"))

df["gene_name_clean"] = df["gene_name"].astype(str).str.strip().str.upper()
filtered_df = df[df["splicing_class"].isin(["high", "low"])].copy()

gwas_categories = {
    "immune_all": "immune_gwas_genes_intersect.bed",
    "immune_antiviral": "antiviral_genes.bed",
    "immune_proviral": "proviral_genes.bed",
    "immune_non_anti/proviral": "non_antiviral_genes.bed",
    "random_genes": "random_background_genes.bed"
}

all_results = []

for category, filepath in gwas_categories.items():
    gwas_genes = pd.read_csv(filepath, sep="\t", header=None, usecols=[3], names=["gene_name"])
    gwas_genes_clean = gwas_genes["gene_name"].astype(str).str.strip().str.upper()
    
    df_sub = filtered_df.copy()
    df_sub["gwas_overlap"] = df_sub["gene_name_clean"].isin(gwas_genes_clean)
    
    contingency = pd.crosstab(df_sub["splicing_class"], df_sub["gwas_overlap"])
    contingency = contingency.reindex(index=["high", "low"], columns=[True, False], fill_value=0)
    
    table = [[contingency.loc["high", True] +1, contingency.loc["low", True] +1],
             [contingency.loc["high", False] + 1, contingency.loc["low", False] +1]]
    
    try:
        oddsratio, p = fisher_exact(table)
        ci_low, ci_high = Table2x2(table).oddsratio_confint()
    except Exception:
        oddsratio, p, ci_low, ci_high = None, None, None, None
    
    result = pd.DataFrame({
        "Category": [category] * 2,
        "Splicing_Category": ["highly_spliced", "lowly_spliced"],
        "GWAS_Overlap_True": [contingency.loc["high", True], contingency.loc["low", True]],
        "GWAS_Overlap_False": [contingency.loc["high", False], contingency.loc["low", False]],
        "Odds_Ratio": [oddsratio] * 2,
        "CI_Lower": [ci_low] * 2,
        "CI_Upper": [ci_high] * 2,
        "P_Value": [p] * 2
    })
    all_results.append(result)

final_df = pd.concat(all_results, ignore_index=True)
final_df.to_csv("fisher_tests_gwas_IMMUNE.tsv", sep="\t", index=False)
