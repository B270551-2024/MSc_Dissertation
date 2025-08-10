#-------------------------
# Kendall tau scatter plots comparing RefSeq vs Ensembl isoform counts
# Modified to show smoothed scatter plots with Ensembl on y-axis, RefSeq on x-axis
# Created by B270551
# Modified: 29/07/2025 
#-------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import kendalltau

ensembl = pd.read_csv("ensembl_isoforms.tsv", sep="\t")
refseq = pd.read_csv("refseq_isoforms.tsv", sep="\t")

ensembl = ensembl[~ensembl["biotype"].str.contains("pseudogene", case=False, na=False)]
refseq = refseq[~refseq["biotype"].str.contains("pseudogene", case=False, na=False)]

ensembl["coding_status"] = ensembl["biotype"].apply(lambda x: "protein_coding" if x == "protein_coding" else "noncoding")

merged = pd.merge(ensembl[["gene_name", "num_isoforms", "coding_status"]], refseq[["gene_name", "num_isoforms"]], on="gene_name", suffixes=("_ensembl", "_refseq"))

sns.set(style="whitegrid")

axis_min = min(np.log(merged["num_isoforms_refseq"].min()), np.log(merged["num_isoforms_ensembl"].min()))
axis_max = max(np.log(merged["num_isoforms_refseq"].max()), np.log(merged["num_isoforms_ensembl"].max()))

def plot_kendall_smoothscatter(df, gene_type, filename):
    df_sub = df[df["coding_status"] == "protein_coding"] if gene_type == "protein_coding" else df[df["coding_status"] != "protein_coding"]
    title = "Isoform Counts: Protein Coding Genes" if gene_type == "protein_coding" else "Isoform Counts: Non-Protein Coding Genes"
    
    tau, p_value = kendalltau(df_sub["num_isoforms_refseq"], df_sub["num_isoforms_ensembl"])
    x_data = np.log(df_sub["num_isoforms_refseq"])
    y_data = np.log(df_sub["num_isoforms_ensembl"])
    
    plt.figure(figsize=(8, 8))
    plt.scatter(x_data, y_data, s=8, c='black', alpha=0.5, edgecolors='none', rasterized=True)
    sns.regplot(x=x_data, y=y_data, scatter=False, color='red', line_kws={'linewidth': 2}, ci=None, ax=plt.gca())
    
    plt.xlim(axis_min, axis_max)
    plt.ylim(axis_min, axis_max)
    plt.plot([axis_min, axis_max], [axis_min, axis_max], 'k--', alpha=0.5, linewidth=1)
    
    plt.title(title, fontsize=16, weight='bold')
    plt.xlabel("log(RefSeq Isoforms)", fontsize=13, weight='bold')
    plt.ylabel("log(Ensembl Isoforms)", fontsize=13, weight='bold')
    
    p_text = "p < 0.001" if p_value < 0.001 else f"p = {p_value:.3f}"
    plt.text(0.05, 0.95, f"Kendall's τ = {tau:.3f}\n{p_text}", horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes, fontsize=12, weight='bold', bbox=dict(facecolor='white', alpha=0.8, edgecolor='black', boxstyle='round,pad=0.5'))
    plt.text(0.05, 0.05, f"n = {len(df_sub):,}", horizontalalignment='left', verticalalignment='bottom', transform=plt.gca().transAxes, fontsize=11, bbox=dict(facecolor='white', alpha=0.8, edgecolor='black', boxstyle='round,pad=0.3'))
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    plt.clf()
    
    return tau, p_value

print("Generating Kendall tau scatter plots...")

tau_coding, p_coding = plot_kendall_smoothscatter(merged, "protein_coding", "scatter_kendall_protein_coding_smoothed.png")
print(f"Protein coding genes: τ = {tau_coding:.3f}, p = {p_coding:.2e}")

tau_noncoding, p_noncoding = plot_kendall_smoothscatter(merged, "non_protein_coding", "scatter_kendall_non_protein_coding_smoothed.png")
print(f"Non-protein coding genes: τ = {tau_noncoding:.3f}, p = {p_noncoding:.2e}")

print("Plots saved successfully!")
