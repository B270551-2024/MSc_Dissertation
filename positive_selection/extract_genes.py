###---------------------
# File for extracting + and - selected genes for each clade (from supp table 6)
# Created by B270551
# 22/05/2025 
###--------------------




import pandas as pd

df = pd.read_csv("supplementary_table6.txt", sep="\t", encoding="latin1")

# Define clade columns of interest
clade_map = {
    "Primates": "Primates",
    "Rodentia": "Rodentia",
    "Chiroptera": "Chiroptera",
    "Artiodactyla": "Artiodactyla",
    "Afrotheria": "Afrotheria",
    "Carnivora": "Carnivora",
    "Eulipotyphla": "Eulipotyphla",
    "Lagomorpha": "Lagomorpha",
    "Perissodactyla": "Perissodactyla",
    "Pholidota": "Pholidota",
    "Humans": "hg38",
    "Mice": "mm10",
    "Bat": "HLrhiSin1" # Bat = chinese horsehoe bat 
}

# Loop through each clade and extract gene list
for clade, col_name in clade_map.items():
    if col_name not in df.columns:
        raise ValueError(f"Column '{col_name}' not found in file. Please check the header row.")

    selected = df[df[col_name] == 1]
    gene_list_pos = selected["geneSymbol"].dropna().drop_duplicates().sort_values()

    # Save to file
    outname_pos = f"positive_genes_{clade}.tsv"
    gene_list_pos.to_csv(outname_pos, index=False, header=False)
    print(f"Saved {len(gene_list_pos)} genes to {outname_pos}")



# Now master gene list for comparison
gene_list = df.iloc[:, 1].dropna().drop_duplicates().sort_values()
output_file = "positive_genes_ALL.tsv"
gene_list.to_csv(output_file, index=False, header=False)
