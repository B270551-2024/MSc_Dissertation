#--------------------
# Script for Parsing gtf files 
# Output = tsv file: gene name, num_isoforms, biotype, source 
# Made by Max Falk
# 20/05/2025 
#---------------------


import pandas as pd
from collections import defaultdict
import gzip
import argparse

def parse_gtf(gtf_file):
    gene_transcripts = defaultdict(set)
    gene_biotype = {}

    open_func = gzip.open if gtf_file.endswith(".gz") else open

    with open_func(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if parts[2] not in ('transcript', 'exon'):
                continue
            info = parts[8]
            attrs = {k.strip(): v.strip('"') for k, v in 
                     [item.strip().split(' ') for item in info.strip(';').split(';') if item]}
            
            gene_id = attrs.get('gene_id')
            transcript_id = attrs.get('transcript_id')
            gene_name = attrs.get('gene_name', gene_id)
            biotype = attrs.get('gene_biotype') or attrs.get('gene_type') or 'unknown' #Might need to be altered to include biotype equivalent in refseq? Currently just = unknown

            if gene_id and transcript_id:
                gene_transcripts[gene_name].add(transcript_id)
                gene_biotype[gene_name] = biotype

    df = pd.DataFrame([
        {"gene_name": gene, "num_isoforms": len(tx), "biotype": gene_biotype.get(gene)}
        for gene, tx in gene_transcripts.items()
    ])
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = parse_gtf(args.gtf)
    df["source"] = "ensembl" if "ensembl" in args.gtf.lower() else "refseq"
    df.to_csv(args.output, sep='\t', index=False)
