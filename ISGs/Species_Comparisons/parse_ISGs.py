import pandas as pd

# Load your input data (assuming tab-separated)
df = pd.read_csv("shaw_ISG_data.txt", sep="\t")

# Filter genes with Log2FC > 2.0
df_filtered = df[df['Log2FC'] > 2.0]

# Group by species and save each group to a file
for species, group in df_filtered.groupby("Species"):
    gene_names = group["Gene Name"].dropna().unique()
    filename = f"{species}_ISGs.txt".replace(" ", "_")
    with open(filename, "w") as f:
        for gene in gene_names:
            f.write(f"{gene}\n")
    print(f"✅ Saved {filename} with {len(gene_names)} genes.")
