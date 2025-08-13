#-------------------------------------------#
# A general quick script that extracts the  top X% of splice isoforms from list of annotated genes
# Created by B270551
# Example use top_percent(input file, output file, %)
# 04/06/2025 
#-----------------------------#

import pandas as pd

def top_percent(file_in, file_out, pct=10):
    df = pd.read_csv(file_in, sep='\t')
    df['num_isoforms'] = pd.to_numeric(df['num_isoforms'], errors='coerce')
    df.dropna(subset=['num_isoforms'], inplace=True)
    n = max(int(len(df) * pct / 100), 1)
    top_genes = df.sort_values('num_isoforms', ascending=False).head(n)
    top_genes[['gene_name']].to_csv(file_out, sep='\t', index=False, header=False)

# Usage example:
# Need to alter input file, output file and % threshold
top_percent('genes_overlapping_all_gwas_annotated.tsv', 'top_10_percent_gwas_genes.tsv', 10)

