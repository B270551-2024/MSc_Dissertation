## Updated Fisher test script with downsampling to match smallest dataset size
## Runs 100 iterations with random subsets and computes averages
## Tests significance loss when sample sizes are reduced
## Created by B270551 - 18/07/2025
##-----------------------------------------------------------##

import pandas as pd
import numpy as np
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
all_pos_file, neg_all_file = "positive_genes_ALL.tsv", "not_positive_genes_NONE.txt"
threshold, n_iterations = 0.1, 100 #ALTER THRESHOLD HERE 

iso_df = pd.read_csv(isoform_file, sep="\t")
iso_df = iso_df[(iso_df["biotype"] == "protein_coding") & (iso_df["num_isoforms"].notna())]

min_size = len(pd.read_csv(clade_files["Human"]["pos"], header=None, names=["gene_name"]))
print(f"Smallest dataset size (Human positive genes): {min_size}")

def fisher_test(df_base, pos_genes, neg_genes, label, thresh):
    pos_set, neg_set = set(pos_genes["gene_name"]), set(neg_genes["gene_name"])
    df = df_base[df_base["gene_name"].isin(pos_set.union(neg_set))].copy()
    df["is_selected"] = df["gene_name"].isin(pos_set)
    
    q_low, q_high = df["num_isoforms"].quantile(thresh), df["num_isoforms"].quantile(1 - thresh)
    df["splicing_class"] = df["num_isoforms"].apply(lambda x: "high" if x >= q_high else ("low" if x <= q_low else "mid"))
    df_extremes = df[df["splicing_class"] != "mid"]
    
    ct = pd.crosstab(df_extremes["is_selected"], df_extremes["splicing_class"])
    hs, ls = ct.get("high", {}).get(True, 0), ct.get("low", {}).get(True, 0)
    hns, lns = ct.get("high", {}).get(False, 0), ct.get("low", {}).get(False, 0)
    
    table = [[hs, ls], [hns, lns]]
    try:
        odds, p = fisher_exact(table)
        ci_low, ci_high = Table2x2(table).oddsratio_confint()
    except:
        odds = p = ci_low = ci_high = None
    
    return {"clade": label, "high_selected": hs, "low_selected": ls, "high_not_selected": hns, 
            "low_not_selected": lns, "odds_ratio": odds, "p_value": p, "ci_low": ci_low, "ci_high": ci_high}

results_all = []
for i in range(n_iterations):
    for clade, files in clade_files.items():
        pos_genes = pd.read_csv(files["pos"], header=None, names=["gene_name"])
        neg_genes = pd.read_csv(files["neg"], header=None, names=["gene_name"])
        
        if len(pos_genes) > min_size:
            pos_genes = pos_genes.sample(n=min_size, random_state=i)
        
        result = fisher_test(iso_df, pos_genes, neg_genes, clade, threshold)
        result["iteration"] = i + 1
        results_all.append(result)

    all_pos = pd.read_csv(all_pos_file, header=None, names=["gene_name"])
    neg_all = pd.read_csv(neg_all_file, header=None, names=["gene_name"])
    
    if len(all_pos) > min_size:
        all_pos = all_pos.sample(n=min_size, random_state=i)
    
    result_all = fisher_test(iso_df, all_pos, neg_all, "ALL", threshold)
    result_all["iteration"] = i + 1
    results_all.append(result_all)

results_df = pd.DataFrame(results_all)
summary_results = []

for clade in results_df["clade"].unique():
    cd = results_df[results_df["clade"] == clade]
    vo, vcl, vch = cd["odds_ratio"].dropna(), cd["ci_low"].dropna(), cd["ci_high"].dropna()
    pv = cd["p_value"].dropna()
    
    summary_results.append({
        "clade": clade, "mean_odds_ratio": vo.mean() if not vo.empty else None,
        "se_odds_ratio": vo.std() / np.sqrt(len(vo)) if not vo.empty else None,
        "mean_ci_low": vcl.mean() if not vcl.empty else None,
        "mean_ci_high": vch.mean() if not vch.empty else None,
        "significant_p_count": len(pv[pv < 0.05]) if not pv.empty else 0,
        "total_iterations": n_iterations
    })

summary_df = pd.DataFrame(summary_results)
summary_df.to_csv("fisher_test_downsampled_summary_20.tsv", sep="\t", index=False)
print("Summary results saved to 'fisher_test_downsampled_summary_20.tsv'")
print("\nSummary Results:")
print(summary_df.to_string(index=False))
