## ISG Splice Enrichment Analysis with Bootstrap Downsampling
## Compares multiple ISG lists for enriched alternative splicing using bootstrap Fisher tests
## Tests enrichment in highly spliced genes with balanced sampling approach
## Usage: Ensure all ISG annotation files are in working directory
## Created by B270551
##---------------------------------------------------------##

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2

background_file = "ensembl_isoforms.tsv"
splicing_threshold_percent = 0.1
bootstrap_iterations = 100

bg_df = pd.read_csv(background_file, sep="\t")
bg_df = bg_df[(bg_df["biotype"] == "protein_coding") & (bg_df["num_isoforms"].notna())]
bg_df["gene_name"] = bg_df["gene_name"].str.upper()

file_list = [
    ("random_200_protein_coding.tsv", "Random_Non_ISGs"),
    ("positive_genes_ALL_annotated.tsv", "Positive_Genes_All"),
    ("merged_unique_ISGs_annotated.tsv", "Mammalian_Interferome"),
    ("human_ISGs_annotated.tsv", "Human_ISGs"),
    ("mouse_ISGs_annotated.tsv", "Mouse_ISGs"),
    ("fruitbat_ISGs_annotated.tsv", "FruitBat_ISGs"),
    ("positively_selected_ISGs_annotated.tsv", "Positively_Selected_ISGs"),
    ("core_ISGs_vert-mamm_annotated.tsv", "Core_Vert_Mamm_ISGs"),
    ("antiviral_genes_all_annotated.tsv", "Antiviral_Genes_All"),
    ("core_antiviral_ISGs_annotated.tsv", "Core_Antiviral_ISGs")
]

gene_dfs = []
for filename, label in file_list:
    df = pd.read_csv(filename, sep="\t", header=None, names=["gene_name", "num_isoforms"])
    df["gene_name"] = df["gene_name"].str.upper()
    df = df.dropna()
    gene_dfs.append((df, label))

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
    
    high_group = extremes[extremes["splicing_class"] == "high"]
    low_group = extremes[extremes["splicing_class"] == "low"]
    
    high_target_orig = high_group["is_target"].sum()
    low_target_orig = low_group["is_target"].sum()
    high_bg_orig = len(high_group) - high_target_orig
    low_bg_orig = len(low_group) - low_target_orig
    
    min_group_size = min(len(high_group), len(low_group))
    
    if high_target_orig < 12 or low_target_orig < 3:
        table = [[high_target_orig + 1, low_target_orig + 1], [high_bg_orig + 1, low_bg_orig + 1]]
        odds, p = fisher_exact(table)
        ci_low, ci_high = Table2x2(table).oddsratio_confint()
        return {
            "clade": label, "mean_odds_ratio": odds, "se_odds_ratio": None,
            "mean_ci_low": ci_low, "mean_ci_high": ci_high,
            "significant_p_count": 1 if p < 0.05 else 0, "total_iterations": 1
        }
    
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
        "clade": label,
        "mean_odds_ratio": np.mean(odds_ratios),
        "se_odds_ratio": np.std(odds_ratios) / np.sqrt(len(odds_ratios)),
        "mean_ci_low": np.percentile(odds_ratios, 2.5),
        "mean_ci_high": np.percentile(odds_ratios, 97.5),
        "significant_p_count": np.sum(p_values < 0.05),
        "total_iterations": len(odds_ratios)
    }

results = []
for df, label in gene_dfs:
    res = bootstrap_fisher_test(df, bg_df[["gene_name", "num_isoforms"]], label, splicing_threshold_percent, bootstrap_iterations)
    results.append(res)

results_df = pd.DataFrame(results)
results_df = results_df[["clade", "mean_odds_ratio", "se_odds_ratio", "mean_ci_low", "mean_ci_high", "significant_p_count", "total_iterations"]]
results_df.to_csv("fisher_tests_ISGs.tsv", sep="\t", index=False)
print("\nSaved results to 'fisher_tests_ISGs.tsv'")
print(results_df.to_string(index=False))
