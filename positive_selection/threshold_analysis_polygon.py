##----------POSITIVE SELECTION ANALYSIS: THRESHOLD SWEEP (ALL GENES)------##
## Testing whether positively selected genes are enriched in highly spliced genes
## across different quantile thresholds (10% to 50%) - Ensembl only
## Created by B270551 - Polygon plots with CI shading + significance stars
##---------------------------------------------------------##

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2

threshold_range = np.arange(0.1, 0.55, 0.05)
all_pos_file, non_pos_file = "positive_genes_ALL.tsv", "not_positive_genes_NONE.txt"
isoform_file = "ensembl_isoforms.tsv"

all_pos_genes = pd.read_csv(all_pos_file, header=None, names=["gene_name"])
pos_set = set(all_pos_genes["gene_name"])
non_pos_genes = pd.read_csv(non_pos_file, header=None, names=["gene_name"])
neg_set = set(non_pos_genes["gene_name"])

df = pd.read_csv(isoform_file, sep="\t")
df = df[(df["biotype"] == "protein_coding") & (df["num_isoforms"].notna()) & (df["gene_name"].isin(pos_set.union(neg_set)))]
df["is_selected"] = df["gene_name"].isin(pos_set)

results = []
for threshold_percent in threshold_range:
    q_low, q_high = df["num_isoforms"].quantile(threshold_percent), df["num_isoforms"].quantile(1 - threshold_percent)
    df_temp = df.copy()
    df_temp["splicing_class"] = df_temp["num_isoforms"].apply(lambda x: "high" if x >= q_high else ("low" if x <= q_low else "mid"))
    df_extremes = df_temp[df_temp["splicing_class"] != "mid"]
    
    ct = pd.crosstab(df_extremes["is_selected"], df_extremes["splicing_class"])
    hs, ls = ct.get("high", {}).get(True, 0), ct.get("low", {}).get(True, 0)
    hns, lns = ct.get("high", {}).get(False, 0), ct.get("low", {}).get(False, 0)
    
    table = [[hs, ls], [hns, lns]]
    try:
        odds, p = fisher_exact(table)
        ci_low, ci_high = Table2x2(table).oddsratio_confint()
    except:
        odds = p = ci_low = ci_high = None
    
    results.append({"threshold_percent": threshold_percent, "odds_ratio": odds, "p_value": p, "ci_low": ci_low, "ci_high": ci_high})

results_df = pd.DataFrame(results)
results_df.to_csv("threshold_sweep_results_Ensembl_ALL.tsv", sep="\t", index=False)

x_vals, y_vals = results_df["threshold_percent"] * 100, results_df["odds_ratio"]
plt.figure(figsize=(8, 5))
plt.errorbar(x_vals, y_vals, yerr=[y_vals - results_df["ci_low"], results_df["ci_high"] - y_vals], 
             fmt='o-', capsize=4, label="Odds Ratio", color='steelblue')
plt.fill_between(x_vals, results_df["ci_low"], results_df["ci_high"], color='lightblue', alpha=0.3, label="95% CI")

for x, y, pval in zip(x_vals, y_vals, results_df["p_value"]):
    if pd.notnull(pval):
        star = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else ''
        if star:
            plt.text(x, y + 0.05, star, ha='center', va='bottom', fontsize=10, color='black')

plt.axhline(1, color='grey', linestyle='--', linewidth=1)
plt.xlabel("Splicing Threshold (%)", fontsize=12, weight='bold')
plt.ylabel("Odds Ratio", fontsize=12, weight='bold')
plt.title("Positive Selection vs Splicing (ALL) — Threshold Sweep", fontsize=14, weight='bold')
plt.grid(True, linestyle='--', alpha=0.4)
plt.legend()
plt.tight_layout()
plt.savefig("fisher_threshold_sweep_POLYGON_ALL_V2_Ensembl.png", dpi=300)
plt.close()

print("Done. Results and plots saved.")
