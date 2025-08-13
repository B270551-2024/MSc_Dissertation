#---------------------------#
# A general script used to annotate gene lists with their splice isoform numbers
# Outputs = e.g. inputfile_annotated.tsv
# Created by B270551
# 10/06/2025
#-----------------------


import csv
import os

# List of gene list files to annotate
input_files = [
    "genes_overlapping_all_gwas.txt" #INPUT FILE
]

# Load isoform data from ensembl_isoforms.tsv
isoform_dict = {}
with open("ensembl_isoforms.tsv", newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        isoform_dict[row["gene_name"]] = row["num_isoforms"]

# Annotate each file
for input_file in input_files:
    base, ext = os.path.splitext(input_file)
    output_file = f"{base}_annotated.tsv"

    with open(input_file) as infile, open(output_file, "w") as out:
        out.write("gene_name\tnum_isoforms\n")
        for line in infile:
            gene = line.strip().split('\t')[0]
            isoforms = isoform_dict.get(gene, "NA")
            out.write(f"{gene}\t{isoforms}\n")

    print(f"Annotated: {output_file}")
