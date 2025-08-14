## GO Enrichment Validation for Pro/Antiviral Gene Lists
## Validates antiviral and proviral gene lists using GO:BP enrichment analysis
## Creates comparative bar plot showing top enriched terms for each list
## Usage: Ensure first_100_proviral_annotated.tsv and first_100_antiviral_annotated.tsv are in working directory
## Created by B270551
##---------------------------------------------------------##

import pandas as pd
from gprofiler import GProfiler
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

proviral_df = pd.read_csv("first_100_proviral_annotated.tsv", sep="\t", usecols=[0], names=["Gene"], header=0)
antiviral_df = pd.read_csv("first_100_antiviral_annotated.tsv", sep="\t", usecols=[0], names=["Gene"], header=0)
proviral_genes = proviral_df["Gene"].dropna().tolist()
antiviral_genes = antiviral_df["Gene"].dropna().tolist()

gp = GProfiler(return_dataframe=True)

enrich_proviral = gp.profile(organism='hsapiens', query=proviral_genes, sources=["GO:BP"])
enrich_antiviral = gp.profile(organism='hsapiens', query=antiviral_genes, sources=["GO:BP"])

top10_proviral = enrich_proviral[enrich_proviral["p_value"] < 0.05].sort_values("p_value").head(10)
top10_antiviral = enrich_antiviral[enrich_antiviral["p_value"] < 0.05].sort_values("p_value").head(10)

cols_to_report = ["native", "name", "p_value", "term_size", "query_size", "intersection_size"]
top10_proviral[cols_to_report].to_csv("top10_GO_BP_proviral.tsv", sep="\t", index=False)
top10_antiviral[cols_to_report].to_csv("top10_GO_BP_antiviral.tsv", sep="\t", index=False)

top10_proviral["List"] = "Proviral"
top10_antiviral["List"] = "Antiviral"

combined = pd.concat([top10_proviral, top10_antiviral], ignore_index=True)
combined["-log10(padj)"] = -np.log10(combined["p_value"])

plt.figure(figsize=(12, 7))
sns.barplot(data=combined, x="-log10(padj)", y="name", hue="List")
plt.xlabel("-log10(adjusted p-value)")
plt.ylabel("GO:BP Term")
plt.title("Top 10 Enriched GO:BP Terms per Gene List")
plt.tight_layout()
plt.gca().invert_yaxis()

plt.savefig("go_enrichment_barplot.png", dpi=300)
plt.show()
