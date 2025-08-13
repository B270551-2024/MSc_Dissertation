## GO Enrichment Bubble Plot Analysis
## Creates bubble plot visualising top GO terms from PANTHER enrichment analysis
## Bubble size represents gene count, colour represents significance (-log10 FDR)
## Usage: Ensure PANTHER output file is in working directory, update filename and parameters (see inline comments)
## Need to user panther overrepresentation test + export table & use this as input file
## Created by B270551
##---------------------------------------------------------##

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

with open("immune_gwas_antiviral_spliced_vs_immune_gwas_antiviral.txt") as f: # INPUT FILE
    for i, line in enumerate(f):
        if line.startswith("PANTHER GO-Slim Biological Process"):
            skiprows = i
            break

df = pd.read_csv("immune_gwas_antiviral_spliced_vs_immune_gwas_antiviral.txt", sep="\t", skiprows=skiprows) # INPUT FILE

df.rename(columns={
    "PANTHER GO-Slim Biological Process": "GO_term",
    "Client Text Box Input (fold Enrichment)": "FoldEnrichment",
    "Client Text Box Input (raw P-value)": "Pvalue",
    "Client Text Box Input (FDR)": "FDR",
    "Client Text Box Input (64)": "GeneCount" # Need to alter this so that it reflects number of genes input to Panther
}, inplace=True)

df["FoldEnrichment"] = pd.to_numeric(df["FoldEnrichment"], errors='coerce')
df["FDR"] = pd.to_numeric(df["FDR"], errors='coerce')
df["GeneCount"] = pd.to_numeric(df["GeneCount"], errors='coerce')

df = df.dropna(subset=["FoldEnrichment", "FDR", "GeneCount"])
df = df[df["GeneCount"] != 1]

top_df = df.sort_values("FoldEnrichment", ascending=False).head(20)
top_df = top_df.sort_values("FoldEnrichment", ascending=True).reset_index(drop=True)
top_df["logFDR"] = -np.log10(top_df["FDR"])

plt.figure(figsize=(8,6), constrained_layout=True)
scatter = plt.scatter(
    x=top_df["FoldEnrichment"],
    y=top_df.index + 1,
    s=top_df["GeneCount"] * 10,
    c=top_df["logFDR"],
    cmap="Reds",
    edgecolor='black',
    alpha=0.8
)

plt.yticks(ticks=top_df.index + 1, labels=top_df["GO_term"])
plt.xlabel("Fold Enrichment")
plt.ylabel("")
plt.title("Top 20 GO Terms ( GWAS - Antiviral (Spliced) vs GWAS - Antiviral)") # TITLE

cbar = plt.colorbar(scatter)
cbar.set_label("-log10(FDR)")

for size in [min(top_df["GeneCount"]), np.median(top_df["GeneCount"]), max(top_df["GeneCount"])]:
    plt.scatter([], [], s=size * 10, c='gray', alpha=0.6, edgecolors='black', label=f"{int(size)} genes")
plt.legend(scatterpoints=1, frameon=True, labelspacing=1, title="Gene Count", loc='lower right')

plt.savefig("GO_enrichment_bubbleplot.png", dpi=300, bbox_inches='tight')
plt.show()