## CRISPR Bootstrap Fisher Test Analysis
## Performs bootstrap Fisher tests with downsampling for antiviral/proviral gene sets
## Tests enrichment in highly spliced genes with balanced high/low splicing groups
## Usage: Ensure all input files are in working directory, run script
## Created by B270551
##---------------------------------------------------------##

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2

background_file = "ensembl_isoforms.tsv"
antiviral_file = "ranked_antiviral_annotated.tsv"
proviral_file = "ranked_proviral_annotated.tsv"
top100_antiviral_file = "first_100_antiviral_annotated.tsv"
top100_proviral_file = "first_100_proviral_annotated.tsv"
av_pos_file = "antiviral_positive_annotated.tsv"
pv_pos_file = "proviral_positive_annotated.tsv"
negative_control = "random_100_protein_coding.tsv"
splicing_threshold_percent = 0.1 #CAN BE ALTERED
bootstrap_iterations = 100

bg_df = pd.read_csv(background_file, sep="\t")
bg_df = bg_df[(bg_df["biotype"] == "protein_coding") & (bg_df["num_isoforms"].notna())]

av_df = pd.read_csv(antiviral_file, sep="\t", header=None, names=["gene_name", "num_isoforms"]).dropna(subset=["num_isoforms"])
pv_df = pd.read_csv(proviral_file, sep="\t", header=None, names=["gene_name", "num_isoforms"]).dropna(subset=["num_isoforms"])
first100_av_df = pd.read_csv(top100_antiviral_file, sep="\t", header=None, names=["gene_name", "num_isoforms"]).dropna(subset=["num_isoforms"])
first100_pv_df = pd.read_csv(top100_proviral_file, sep="\t", header=None, names=["gene_name", "num_isoforms"]).dropna(subset=["num_isoforms"])
av_pos_df = pd.read_csv(av_pos_file, sep="\t", header=None, names=["gene_name", "num_isoforms"]).dropna(subset=["num_isoforms"])
pv_pos_df = pd.read_csv(pv_pos_file, sep="\t", header=None, names=["gene_name", "num_isoforms"]).dropna(subset=["num_isoforms"])
control_df = pd.read_csv(negative_control, sep="\t", header=None, names=["gene_name", "num_isoforms"]).dropna(subset=["num_isoforms"])

av_pv_genes = set(av_df["gene_name"]).union(set(pv_df["gene_name"]))
control_df = control_df[~control_df["gene_name"].isin(av_pv_genes)]

def bootstrap_fisher_test(target_genes, background_genes, label, threshold, n_iter=100, random_state=42):
    np.random.seed(random_state)
    
    background_genes = background_genes.drop_duplicates(subset="gene_name").copy()
    background_genes["is_target"] = background_genes["gene_name"].isin(target_genes["gene_name"])
    
    q_low = background_genes["num_isoforms"].quantile(threshold)
    q_high = background_genes["num_isoforms"].quantile(1 - threshold)
    
    background_genes["splicing_class"] = background_genes["num_isoforms"].apply(
        lambda x: "high" if x >= q_high else ("low" if x <= q_low else "mid")
    )
    
    extremes = background_genes[background_genes["splicing_class"] != "mid"].copy()
    
    print(f"\n[{label}] Splicing class counts (full extremes):")
    print(extremes["splicing_class"].value_counts())
    
    high_group = extremes[extremes["splicing_class"] == "high"]
    low_group = extremes[extremes["splicing_class"] == "low"]
    
    n_high, n_low = len(high_group), len(low_group)
    min_group_size = min(n_high, n_low)
    
    odds_ratios, p_values = [], []
    
    for i in range(n_iter):
        if n_high > n_low:
            high_sample = high_group.sample(n=min_group_size, replace=False)
            low_sample = low_group.copy()
        elif n_low > n_high:
            low_sample = low_group.sample(n=min_group_size, replace=False)
            high_sample = high_group.copy()
        else:
            high_sample, low_sample = high_group.copy(), low_group.copy()
        
        high_target, high_bg = high_sample["is_target"].sum(), len(high_sample) - high_sample["is_target"].sum()
        low_target, low_bg = low_sample["is_target"].sum(), len(low_sample) - low_sample["is_target"].sum()
        
        table = [[high_target, low_target], [high_bg, low_bg]]
        
        if any(sum(row) == 0 for row in table) or any(sum(col) == 0 for col in zip(*table)):
            continue
        
        odds, p = fisher_exact(table)
        odds_ratios.append(odds)
        p_values.append(p)
    
    odds_ratios, p_values = np.array(odds_ratios), np.array(p_values)
    
    return {
        "clade": label,
        "mean_odds_ratio": np.mean(odds_ratios),
        "se_odds_ratio": np.std(odds_ratios, ddof=1),
        "mean_ci_low": np.percentile(odds_ratios, 2.5),
        "mean_ci_high": np.percentile(odds_ratios, 97.5),
        "significant_p_count": np.sum(p_values < 0.05),
        "total_iterations": len(odds_ratios)
    }

results = []
for (target, label) in [
    (av_df, "Antiviral_Genes"),
    (pv_df, "Proviral_Genes"),
    (first100_av_df, "Antiviral_Genes_Top100_MAIC"),
    (first100_pv_df, "Proviral_Genes_Top100_MAIC"),
    (av_pos_df, "Antiviral_Genes_Positive_Score"),
    (pv_pos_df, "Proviral_Genes_Positive_Score"),
    (control_df, "Negative_Control")
]:
    res = bootstrap_fisher_test(target, bg_df[["gene_name", "num_isoforms"]], label, 
                               splicing_threshold_percent, bootstrap_iterations, 42)
    results.append(res)

results_df = pd.DataFrame(results)
results_df.to_csv("fisher_tests_CRISPR_bootstrap.tsv", sep="\t", index=False)
print("\nBootstrap downsampling Fisher test results saved to 'fisher_tests_CRISPR_bootstrap.FINAL_2407.tsv'")
print(results_df.to_string(index=False))
