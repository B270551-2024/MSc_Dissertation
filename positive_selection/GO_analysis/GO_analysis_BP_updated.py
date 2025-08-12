## GO Biological Process enrichment analysis for positively selected genes across mammalian clades
## Uses g:Profiler with custom background and creates heatmap visualization
## Created by B270551
##---------------------------------------------------------##

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from gprofiler import GProfiler
import numpy as np

clades = {
    "Afrotheria": "positive_genes_Afrotheria.tsv",
    "Primates": "positive_genes_Primates.tsv",
    "Rodents": "positive_genes_Rodentia.tsv",
    "Lagomorpha": "positive_genes_Lagomorpha.tsv",
    "Eulipotyphla": "positive_genes_Eulipotyphla.tsv",
    "Perissodactyla": "positive_genes_Perissodactyla.tsv",
    "Artiodactyla": "positive_genes_Artiodactyla.tsv",
    "Carnivora": "positive_genes_Carnivora.tsv",
    "Pholidota": "positive_genes_Pholidota.tsv",
    "Chiroptera": "positive_genes_Chiroptera.tsv"
}

go_terms_to_display = [
    "GO:0000003", "GO:0002376", "GO:0009987", "GO:0016032",
    "GO:0023052", "GO:0032501", "GO:0032502", "GO:0040011",
    "GO:0044419", "GO:0050896", "GO:0051179", "GO:0065007"
]

try:
    background_genes = pd.read_csv("all_tested_genes.txt", header=None)[0].tolist()
    print(f"Loaded {len(background_genes)} background genes. First 10: {background_genes[:10]}")
except FileNotFoundError:
    print("Error: all_tested_genes.txt not found. Please check the file path.")
    exit()

gp = GProfiler(return_dataframe=True)
results_df = pd.DataFrame(index=go_terms_to_display)

for clade, file in clades.items():
    try:
        genes = pd.read_csv(file, header=None)[0].tolist()
        print(f"\nLoaded {len(genes)} genes for {clade}. First 5: {genes[:5]}")
    except FileNotFoundError:
        print(f"Error: {file} not found. Skipping {clade}.")
        continue
    
    try:
        enrichment = gp.profile(organism='hsapiens', query=genes, sources=['GO:BP'], background=background_genes)
    except Exception as e:
        print(f"Error in g:Profiler for {clade}: {e}")
        continue
    
    print(f"\nEnrichment results for {clade} (top 5):")
    print(enrichment[['native', 'p_value', 'name']].head())
    
    enrichment = enrichment[enrichment['native'].isin(go_terms_to_display)]
    if enrichment.empty:
        print(f"Warning: No enrichment results for {clade} with specified GO terms.")
    enrichment = enrichment.set_index('native')
    results_df[clade] = enrichment['p_value']

print("\nRaw p-value DataFrame:")
print(results_df)

neg_log10_p = results_df.apply(lambda x: -np.log10(x) if pd.notna(x) else np.nan)
print("\n-log10(p) DataFrame:")
print(neg_log10_p)

mask_nonsig = results_df > 0.05

plt.figure(figsize=(10, 6))
sns_heatmap = sns.heatmap(neg_log10_p, annot=True, fmt=".1f", linewidths=0.5, linecolor='gray',
                         cbar_kws={'label': '-log10(p-value)', 'ticks': [0, 10, 20, 30, 40]},
                         vmin=0, vmax=40, square=False, cmap='YlOrRd', mask=neg_log10_p.isna(), cbar=True)

ax = plt.gca()
for (y, x), val in np.ndenumerate(mask_nonsig.values):
    if val and not neg_log10_p.isna().iat[y, x]:
        ax.add_patch(plt.Rectangle((x, y), 1, 1, fill=True, color='lightgrey', lw=0))

plt.title("GO Biological Process Enrichment (Selected Terms)")
plt.tight_layout()
plt.savefig("GO_selected_terms_heatmap_FINAL.png")
plt.show()
