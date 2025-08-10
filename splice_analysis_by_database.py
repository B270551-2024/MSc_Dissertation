#-------------------------
# A script comparing coding vs non-coding genes for each database
# Outputs: mann-whitney (do coding genes have more isoforms than non-coding) // isoforms number comparisons // are coding more likely to be highly spliced barplots
# Created by B270551
# 22/05/2025 
# Modified to include Fisher's exact tests for highly spliced gene proportions
#-------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr, fisher_exact, kendalltau
import numpy as np
import matplotlib.ticker as mtick

percentile_cutoff = 90
sns.set(style="whitegrid")

ensembl = pd.read_csv("ensembl_isoforms.tsv", sep="\t")
refseq = pd.read_csv("refseq_isoforms.tsv", sep="\t")

ensembl = ensembl[~ensembl["biotype"].str.contains("pseudogene", case=False, na=False)]
refseq = refseq[~refseq["biotype"].str.contains("pseudogene", case=False, na=False)]

ensembl["coding_status"] = ensembl["biotype"].apply(lambda x: "protein_coding" if x == "protein_coding" else "noncoding")
ensembl["is_coding"] = ensembl["coding_status"] == "protein_coding"

# Step 1 - Mann Whitney U & boxplots
def compare_isoforms(df, db_label):
    coding = df[df["is_coding"]]["num_isoforms"]
    noncoding = df[~df["is_coding"]]["num_isoforms"]
    if coding.empty or noncoding.empty:
        print(f"\n[{db_label}] Skipping: one or both groups are empty.")
        return {"database": db_label, "p_value": None, "median_coding": None, "median_noncoding": None, "note": "One or both groups empty"}

    stat, p = mannwhitneyu(coding, noncoding, alternative="two-sided")
    print(f"\n[{db_label}] Mann–Whitney U Test:")
    print(f"Number of genes — Coding: {len(coding)}, Noncoding: {len(noncoding)}")
    print(f"p-value = {p:.2e}")
    print(f"Median isoforms — Coding: {coding.median()}, Noncoding: {noncoding.median()}")

    palette = {"protein_coding": "#1f77b4", "noncoding": "#ff7f0e"}
    df["log_isoforms"] = np.log1p(df["num_isoforms"])
    sns.boxplot(data=df, x="coding_status", y="log_isoforms", palette=palette, width=0.5)
    sns.despine()
    plt.title(f"{db_label}: Isoform Counts by Coding Status", fontsize=13, weight='bold')
    plt.ylabel("log1p(Number of Isoforms)")
    plt.xlabel("")
    plt.tight_layout()
    plt.savefig(f"boxplot_{db_label.lower()}_isoforms_coding_vs_noncoding.png", dpi=300)
    plt.clf()

    return {"database": db_label, "p_value": f"{p:.2e}", "median_coding": coding.median(), "median_noncoding": noncoding.median(), "note": ""}

results = [compare_isoforms(ensembl, "Ensembl")]

merged = pd.merge(ensembl[["gene_name", "num_isoforms", "coding_status"]], refseq[["gene_name", "num_isoforms"]], on="gene_name", suffixes=("_ensembl", "_refseq"))
merged["is_coding"] = merged["coding_status"] == "protein_coding"

results.append(compare_isoforms(merged.rename(columns={"num_isoforms_refseq": "num_isoforms"})[["num_isoforms", "coding_status", "is_coding"]].copy(), "RefSeq"))

pd.DataFrame(results).to_csv("mannwhitney_isoform_test_results.tsv", sep="\t", index=False)

corr, p_corr = spearmanr(merged["num_isoforms_ensembl"], merged["num_isoforms_refseq"])
print(f"\n[Ensembl vs RefSeq] Spearman correlation = {corr:.2f}, p = {p_corr:.4g}")

merged["delta_isoforms"] = merged["num_isoforms_ensembl"] - merged["num_isoforms_refseq"]
merged["abs_diff"] = merged["delta_isoforms"].abs()
merged.sort_values("abs_diff", ascending=False).head(100)[["gene_name", "num_isoforms_ensembl", "num_isoforms_refseq", "delta_isoforms", "coding_status"]].to_csv("genes_with_largest_isoform_discrepancy.tsv", sep="\t", index=False)

# Step 2 - scatterplots
plt.figure(figsize=(12, 5))
sns.scatterplot(data=merged, x="num_isoforms_ensembl", y="num_isoforms_refseq", hue="coding_status", alpha=0.6, palette={"protein_coding": "#1f77b4", "noncoding": "#ff7f0e"})
plt.title("Isoform Counts: Ensembl vs RefSeq", fontsize=16, weight='bold')
plt.xlabel("Number of Ensembl Isoforms", fontsize=13, weight='bold')
plt.ylabel("Number of RefSeq Isoforms", fontsize=13, weight='bold')
plt.legend(title="Coding Status", fontsize=11, title_fontsize=12)
plt.tight_layout()
plt.savefig("scatterplot_ensembl_vs_refseq_isoforms_FINAL.png", dpi=300)
plt.clf()

plt.figure(figsize=(9, 4.5))
sns.histplot(data=merged, x="delta_isoforms", kde=True, color='steelblue')
plt.axvline(0, color='red', linestyle='--')
plt.title("Difference in Isoform Counts (Ensembl - RefSeq)", fontsize=16, weight='bold')
plt.xlabel("Isoform Count Difference", fontsize=13, weight='bold')
plt.ylabel("Gene Count", fontsize=13, weight='bold')
plt.xlim(-50, 50)
plt.ylim(0, 4000)
sns.despine()
plt.tight_layout()
plt.savefig("histogram_ensembl_minus_refseq_isoforms.png", dpi=300)
plt.clf()

def perform_fisher_test(df, db_label, cutoff_value, isoform_col="num_isoforms"):
    df["highly_spliced"] = df[isoform_col] > cutoff_value
    contingency = pd.crosstab(df["is_coding"], df["highly_spliced"])
    print(f"\n[{db_label}] Fisher's Exact Test Contingency Table:")
    print("Rows: is_coding (False=non-coding, True=protein-coding)")
    print("Cols: highly_spliced (False=not highly spliced, True=highly spliced)")
    print(contingency)
    
    try:
        a = contingency.loc[True, True]
        b = contingency.loc[True, False]
        c = contingency.loc[False, True]
        d = contingency.loc[False, False]
        
        table = [[a, b], [c, d]]
        odds_ratio, p_value = fisher_exact(table, alternative='two-sided')
        
        prop_coding_highly_spliced = a / (a + b) if (a + b) > 0 else 0
        prop_noncoding_highly_spliced = c / (c + d) if (c + d) > 0 else 0
        
        log_or = np.log(odds_ratio) if odds_ratio > 0 else np.nan
        se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d) if min(a,b,c,d) > 0 else np.nan
        
        if not np.isnan(se_log_or):
            ci_lower = np.exp(log_or - 1.96 * se_log_or)
            ci_upper = np.exp(log_or + 1.96 * se_log_or)
        else:
            ci_lower = ci_upper = np.nan
        
        print(f"\nFisher's Exact Test Results:")
        print(f"Odds Ratio = {odds_ratio:.3f} (95% CI: {ci_lower:.3f} - {ci_upper:.3f})")
        print(f"p-value = {p_value:.2e}")
        print(f"Proportion highly spliced - Protein-coding: {prop_coding_highly_spliced:.1%}")
        print(f"Proportion highly spliced - Non-coding: {prop_noncoding_highly_spliced:.1%}")
        
        return {"database": db_label, "cutoff_percentile": percentile_cutoff, "cutoff_value": cutoff_value, "protein_coding_highly_spliced": a, "protein_coding_not_highly_spliced": b, "noncoding_highly_spliced": c, "noncoding_not_highly_spliced": d, "odds_ratio": odds_ratio, "odds_ratio_ci_lower": ci_lower, "odds_ratio_ci_upper": ci_upper, "p_value": p_value, "prop_coding_highly_spliced": prop_coding_highly_spliced, "prop_noncoding_highly_spliced": prop_noncoding_highly_spliced}
        
    except Exception as e:
        print(f"Error in Fisher's test for {db_label}: {e}")
        return {"database": db_label, "error": str(e), "cutoff_percentile": percentile_cutoff, "cutoff_value": cutoff_value}

fisher_results = []
cutoff_ensembl = np.percentile(ensembl["num_isoforms"], percentile_cutoff)
fisher_results.append(perform_fisher_test(ensembl, "Ensembl", cutoff_ensembl))

cutoff_refseq = np.percentile(merged["num_isoforms_refseq"], percentile_cutoff)
fisher_results.append(perform_fisher_test(merged.rename(columns={"num_isoforms_refseq": "num_isoforms"}), "RefSeq", cutoff_refseq))

pd.DataFrame(fisher_results).to_csv("fisher_exact_test_results.tsv", sep="\t", index=False)
print(f"\nFisher's exact test results saved to: fisher_exact_test_results.tsv")

ensembl["highly_spliced_plot"] = ensembl["num_isoforms"] > cutoff_ensembl
bar_ensembl = ensembl.groupby("coding_status")["highly_spliced_plot"].value_counts(normalize=True).unstack()
bar_ensembl.rename(columns={True: "Highly Spliced", False: "Not or Lowly Spliced"}, inplace=True)
bar_ensembl = bar_ensembl[["Not or Lowly Spliced", "Highly Spliced"]]

colors = ["#cccccc", "#1f77b4"]
ax1 = bar_ensembl.plot(kind="bar", stacked=True, figsize=(7, 5), edgecolor="black", color=colors)

if len(fisher_results) > 0 and 'p_value' in fisher_results[0]:
    p_val = fisher_results[0]['p_value']
    plt.text(0.02, 0.98, f"Fisher's p = {p_val:.2e}", transform=ax1.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left', bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

plt.title(f"Ensembl: Proportion of Highly Spliced Genes\n(> {percentile_cutoff}th Percentile)", fontsize=15, weight='bold')
plt.ylabel("Proportion", fontsize=13, weight='bold')
plt.xlabel("Gene Type", fontsize=13, weight='bold')
plt.xticks(rotation=0, fontsize=11)
plt.yticks(fontsize=11)
plt.legend(title="", fontsize=12, loc='lower right')
ax1.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
sns.despine()
plt.tight_layout()
plt.savefig("barplot_ensembl_highly_spliced.png", dpi=300)
plt.clf()

merged["refseq_spliced"] = merged["num_isoforms_refseq"] > cutoff_refseq
bar_refseq = merged.groupby("coding_status")["refseq_spliced"].value_counts(normalize=True).unstack()
bar_refseq.rename(columns={True: "Highly Spliced", False: "Not or Lowly Spliced"}, inplace=True)
bar_refseq = bar_refseq[["Not or Lowly Spliced", "Highly Spliced"]]

ax2 = bar_refseq.plot(kind="bar", stacked=True, figsize=(7, 5), edgecolor="black", color=colors)

if len(fisher_results) > 1 and 'p_value' in fisher_results[1]:
    p_val = fisher_results[1]['p_value']
    plt.text(0.02, 0.98, f"Fisher's p = {p_val:.2e}", transform=ax2.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left', bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

plt.title(f"RefSeq: Proportion of Highly Spliced Genes\n(> {percentile_cutoff}th Percentile)", fontsize=15, weight='bold')
plt.ylabel("Proportion", fontsize=13, weight='bold')
plt.xlabel("Gene Type", fontsize=13, weight='bold')
plt.xticks(rotation=0, fontsize=11)
plt.yticks(fontsize=11)
plt.legend(title="", fontsize=12, loc='lower right')
ax2.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
sns.despine()
plt.tight_layout()
plt.savefig("barplot_refseq_highly_spliced.png", dpi=300)
plt.clf()

print("\n" + "="*60)
print("ANALYSIS COMPLETE")
print("="*60)
print("Files generated:")
print("- mannwhitney_isoform_test_results.tsv (Mann-Whitney U test results)")
print("- fisher_exact_test_results.tsv (Fisher's exact test results)")
print("- genes_with_largest_isoform_discrepancy.tsv (top 100 discrepant genes)")
print("- Various plots (.png files)")
print("="*60)
