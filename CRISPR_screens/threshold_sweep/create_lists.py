# A script which creates lists of top 100 to top 2500 genes within each MAIC list 
# Done so that I can asses if there's a drop off in significance (i.e. to see whether 100 is a reasonable cut off)
# Created by Max Falk
# 19/07/2025 


import pandas as pd

# --- Parameters ---
proviral_file = "ranked_proviral_annotated.tsv"
antiviral_file = "ranked_antiviral_annotated.tsv"
sizes = [50, 100, 500, 1000, 2000]
output_prefixes = {
    "proviral": "first_{}_proviral_annotated.tsv",
    "antiviral": "first_{}_antiviral_annotated.tsv"
}

# --- Function to split and save gene lists ---
def split_gene_list(input_file, output_prefix, sizes):
    # Read the input file (no headers, tab-separated)
    df = pd.read_csv(input_file, sep="\t", header=None, names=["gene_name", "num_isoforms"])
    
    # Ensure all rows are included, even with NA values
    for size in sizes:
        if size <= len(df):
            output_file = output_prefix.format(size)
            df.head(size).to_csv(output_file, sep="\t", index=False, header=False)
            print(f"Saved {output_file} with {size} genes")

# --- Process both files ---
split_gene_list(proviral_file, output_prefixes["proviral"], sizes)
split_gene_list(antiviral_file, output_prefixes["antiviral"], sizes)
