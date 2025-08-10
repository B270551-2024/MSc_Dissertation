# Final script for parsing data from Ensembl/RefSeq 
# Created by B270551 22/07/2025
#---------------------------

import pandas as pd
from collections import defaultdict
import gzip
import argparse

def parse_gtf(gtf_file):
    gene_types = [
        'protein_coding', 'retained_intron', 'nonsense_mediated_decay', 'protein_coding_CDS_not_defined'
    ]

    gene_info = defaultdict(lambda: {
        "biotype": None,
        "transcript_ids": set(),
        "transcript_types": defaultdict(set)
    })

    open_func = gzip.open if gtf_file.endswith(".gz") else open

    with open_func(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature != 'transcript':
                continue  # only process transcript lines

            # Parse attributes field into dict
            attrs_str = parts[8]
            attrs = {}
            for attr in attrs_str.strip().split(';'):
                attr = attr.strip()
                if not attr:
                    continue
                if ' ' not in attr:
                    continue
                key, value = attr.split(' ', 1)
                attrs[key] = value.strip('"')

            transcript_id = attrs.get('transcript_id')
            if not transcript_id:
                continue  # skip if no transcript_id

            gene_id = attrs.get('gene_id')
            gene_name = attrs.get('gene_name', gene_id)
            transcript_type = attrs.get('transcript_type') or attrs.get('transcript_biotype')
            gene_biotype = attrs.get('gene_biotype') or attrs.get('gene_type')

            # Skip pseudogenes
            if gene_biotype and 'pseudogene' in gene_biotype.lower():
                continue

            gene_info[gene_name]["biotype"] = gene_biotype
            gene_info[gene_name]["transcript_ids"].add(transcript_id)

            # Add minimal change here to put unknown types under 'non_coding'
            if transcript_type in gene_types:
                gene_info[gene_name]["transcript_types"][transcript_type].add(transcript_id)
            else:
                gene_info[gene_name]["transcript_types"]["non_coding"].add(transcript_id)

    return gene_info, gene_types

def build_summary_df(gene_info, gene_types):
    summary_data = []
    for gene_name, info in gene_info.items():
        row = {
            "gene_name": gene_name,
            "gene_biotype": info["biotype"],
            "num_isoforms_total": len(info["transcript_ids"]),
        }
        # Count for each known gene type
        for gene_type in gene_types:
            row[f"num_{gene_type}_isoforms"] = len(info["transcript_types"].get(gene_type, set()))
        # Add count for non_coding isoforms
        row["num_non_coding_isoforms"] = len(info["transcript_types"].get("non_coding", set()))
        summary_data.append(row)

    summary_df = pd.DataFrame(summary_data)
    return summary_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse GTF file to count isoforms by gene type")
    parser.add_argument("--gtf", required=True, help="Path to GTF file (gzipped or plain)")
    parser.add_argument("--output", required=True, help="Path to output TSV file")
    args = parser.parse_args()

    gene_info, gene_types = parse_gtf(args.gtf)
    df = build_summary_df(gene_info, gene_types)
    df["source"] = "ensembl" if "ensembl" in args.gtf.lower() else "refseq"
    df.to_csv(args.output, sep='\t', index=False)
    print(df.head())
