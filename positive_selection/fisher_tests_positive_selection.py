##----------POSITIVE SELECTION ANALYSIS: UPDATED------##
## Testing whether various mammalian clades have significant enrichment in highly spliced variants (for positively selected genes)
## Assessed with a Fisher's test - ENSEMBL only - protein coding
## UPDATED FOR ADDITIONAL CLADE COVERAGE AND CORRECT NEGATIVE SUBSETS
## Created by B270551 - 12/06/2025
##-----------------------------------------------------------##

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2

isoform_file = "ensembl_isoforms.tsv"
clade_files = {
    "Negative Control": {"pos": "positive_genes_CONTROL.tsv", "neg": "not_positive_genes_CONTROL.tsv"},
    "Primates": {"pos": "positive_genes_Primates.tsv", "neg": "not_positive_genes_Primates.tsv"},
    "Rodentia": {"pos": "positive_genes_Rodentia.tsv", "neg": "not_positive_genes_Rodentia.tsv"},
    "Chiroptera": {"pos": "positive_genes_Chiroptera.tsv", "neg": "not_positive_genes_Chiroptera.tsv"},
    "Artiodactyla": {"pos": "positive_genes_Artiodactyla.tsv", "neg": "not_positive_genes_Artiodactyla.tsv"},
    "Afrotheria": {"pos": "positive_genes_Afrotheria.tsv", "neg": "not_positive_genes_Afrotheria.tsv"},
    "Carnivora": {"pos": "positive_genes_Carnivora.tsv", "neg": "not_positive_genes_Carnivora.tsv"},
    "Eulipotyphla": {"pos": "positive_genes_Eulipotyphla.tsv", "neg": "not_positive_genes_Eulipotyphla.tsv"},
    "Lagomorpha": {"pos": "positive_genes_Lagomorpha.tsv", "neg": "not_positive_genes_Lagomorpha.tsv"},
    "Perissodactyla": {"pos": "positive_genes_Perissodactyla.tsv", "neg": "not_positive_genes_Perissodactyla.tsv"},
    "Pholidota": {"pos": "positive_genes_Pholidota.tsv", "neg": "not_positive_genes_Pholidota.tsv"},
    "Human": {"pos": "positive_genes_Humans.tsv", "neg": "not_positive_genes_Humans.tsv"},
    "Mouse": {"pos": "positive_genes_Mice.tsv", "neg": "not_positive_genes_Mice.tsv"},
    "Bat": {"pos": "positive_genes_Bat.tsv", "neg": "not_positive_genes_Bat.tsv"}
}
all_pos_file = "positive_genes_ALL.tsv"
neg_all_file = "not_positive_genes_NONE.txt"
threshold = 0.1

iso_df = pd.read_csv(isoform_file, sep="\t")
iso_df = iso_df[(iso_df["biotype"] == "protein_coding") & (iso_df["num_isoforms"].notna())]

def fisher_test(df_base, pos_genes, neg_genes, label, thresh):
    pos_set, neg_set = set(pos_genes["gene_name"]), set(neg_genes["gene_name"])
    df = df_base[df_base["gene_name"].isin(pos_set.union(neg_set))].copy()
    df["is_selected"] = df["gene_name"].isin(pos_set)
    
    q_low, q_high = df["num_isoforms"].quantile(thresh), df["num_isoforms"].quantile(1 - thresh)
    df["splicing_class"] = df["num_isoforms"].apply(lambda x: "high" if x >= q_high else ("low" if x <= q_low else "mid"))
    df_extremes = df[df["splicing_class"] != "mid"]
    
    df[(df["splicing_class"] == "high") & (df["is_selected"])][["gene_name"]].to_csv(f"high_selected_genes_{label}.tsv", sep="\t", index=False)
    
    ct = pd.crosstab(df_extremes["is_selected"], df_extremes["splicing_class"])
    hs, ls = ct.get("high", {}).get(True, 0), ct.get("low", {}).get(True, 0)
    hns, lns = ct.get("high", {}).get(False, 0), ct.get("low", {}).get(False, 0)
    
    table = [[hs + 1, ls + 1], [hns + 1, lns + 1]]
    try:
        odds, p = fisher_exact(table)
        ci_low, ci_high = Table2x2(table).oddsratio_confint()
    except:
        odds = p = ci_low = ci_high = None
    
    return {"clade": label, "high_selected": hs, "low_selected": ls, "high_not_selected": hns, 
            "low_not_selected": lns, "odds_ratio": odds, "p_value": p, "ci_low": ci_low, "ci_high": ci_high}

results = []
for clade, files in clade_files.items():
    pos_genes = pd.read_csv(files["pos"], header=None, names=["gene_name"])
    neg_genes = pd.read_csv(files["neg"], header=None, names=["gene_name"])
    results.append(fisher_test(iso_df, pos_genes, neg_genes, clade, threshold))

all_pos = pd.read_csv(all_pos_file, header=None, names=["gene_name"])
neg_all = pd.read_csv(neg_all_file, header=None, names=["gene_name"])
results.append(fisher_test(iso_df, all_pos, neg_all, "ALL", threshold))

pd.DataFrame(results).to_csv("fisher_test_results_10_NEW.tsv", sep="\t", index=False)
print("Fisher test results saved to 'fisher_test_results_10_NEW.tsv'")
