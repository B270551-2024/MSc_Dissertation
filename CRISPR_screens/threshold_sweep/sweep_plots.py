import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import fisher_exact

# --- PARAMETERS ---
splicing_threshold = 0.1  # Fixed 10% threshold
bootstrap_iterations = 100  # Number of bootstrap iterations
gene_set_sizes = [100, 500, 1000, 2000]
background_file = "ensembl_isoforms.tsv"
datasets = {
    "Proviral": [
        "first_100_proviral_annotated.tsv",
        "first_500_proviral_annotated.tsv",
        "first_1000_proviral_annotated.tsv",
        "first_2000_proviral_annotated.tsv"
    ],
    "Antiviral": [
        "first_100_antiviral_annotated.tsv",
        "first_500_antiviral_annotated.tsv",
        "first_1000_antiviral_annotated.tsv",
        "first_2000_antiviral_annotated.tsv"
    ]
}

# --- Bootstrap Fisher test function with downsampling ---
def bootstrap_fisher_test(target_genes, background_genes, label, threshold, n_iter=100, random_state=42):
    np.random.seed(random_state)
    
    # Prepare background with is_target flag
    background_genes = background_genes.drop_duplicates(subset="gene_name").copy()
    background_genes["is_target"] = background_genes["gene_name"].isin(target_genes["gene_name"])
    
    # Determine thresholds on full dataset
    q_low = background_genes["num_isoforms"].quantile(threshold)
    q_high = background_genes["num_isoforms"].quantile(1 - threshold)
    
    # Assign splicing classes
    background_genes["splicing_class"] = background_genes["num_isoforms"].apply(
        lambda x: "high" if x >= q_high else ("low" if x <= q_low else "mid")
    )
    
    # Filter only extremes (high and low)
    extremes = background_genes[background_genes["splicing_class"] != "mid"].copy()
    
    print(f"\n[{label}] Splicing class counts (full extremes):")
    print(extremes["splicing_class"].value_counts())
    
    # Split extremes into groups
    high_group = extremes[extremes["splicing_class"] == "high"]
    low_group = extremes[extremes["splicing_class"] == "low"]
    
    # Number in smaller group to downsample larger group
    n_high = len(high_group)
    n_low = len(low_group)
    
    min_group_size = min(n_high, n_low)
    
    # Containers for bootstrap results
    odds_ratios = []
    p_values = []
    
    for i in range(n_iter):
        # Downsample larger group to match smaller group size
        if n_high > n_low:
            high_sample = high_group.sample(n=min_group_size, replace=False)
            low_sample = low_group.copy()
        elif n_low > n_high:
            low_sample = low_group.sample(n=min_group_size, replace=False)
            high_sample = high_group.copy()
        else:
            # Equal sizes, no downsampling needed
            high_sample = high_group.copy()
            low_sample = low_group.copy()
        
        # Build contingency table: rows = is_target (True/False), cols = splicing_class (high/low)
        high_target = high_sample["is_target"].sum()
        high_bg = len(high_sample) - high_target
        low_target = low_sample["is_target"].sum()
        low_bg = len(low_sample) - low_target
        
        # Log contingency table for debugging
        print(f"[{label}] Iteration {i+1}: Contingency table: [[{high_target}, {low_target}], [{high_bg}, {low_bg}]]")
        
        # Add pseudocount to avoid degenerate tables
        pseudocount = 1.0
        table = [[high_target + pseudocount, low_target + pseudocount],
                 [high_bg + pseudocount, low_bg + pseudocount]]
        
        # Run Fisher's exact test
        odds, p = fisher_exact(table)
        odds_ratios.append(odds)
        p_values.append(p)
    
    # Convert results to arrays
    odds_ratios = np.array(odds_ratios)
    p_values = np.array(p_values)
    
    # Calculate bootstrap mean and 95% CI for odds ratio
    mean_odds = np.mean(odds_ratios) if len(odds_ratios) > 0 else np.nan
    ci_low_odds = np.percentile(odds_ratios, 2.5) if len(odds_ratios) > 0 else np.nan
    ci_high_odds = np.percentile(odds_ratios, 97.5) if len(odds_ratios) > 0 else np.nan
    
    # Calculate standard error for odds ratio
    se_odds = np.std(odds_ratios, ddof=1) if len(odds_ratios) > 0 else np.nan
    
    # Count significant p-values (p < 0.05)
    significant_p_count = np.sum(p_values < 0.05) if len(p_values) > 0 else 0
    
    return {
        "clade": label,
        "mean_odds_ratio": mean_odds,
        "se_odds_ratio": se_odds,
        "mean_ci_low": ci_low_odds,
        "mean_ci_high": ci_high_odds,
        "significant_p_count": significant_p_count,
        "total_iterations": len(odds_ratios)
    }

def run_gene_set_sweep(label, gene_files):
    print(f"Running analysis for {label}...")
    
    # Load and filter background
    bg_df = pd.read_csv(background_file, sep="\t")
    bg_df = bg_df[bg_df["biotype"] == "protein_coding"]
    bg_df = bg_df.dropna(subset=["num_isoforms"])
    
    results = []
    
    for size, gene_file in zip(gene_set_sizes, gene_files):
        # Load target gene set
        try:
            target_df = pd.read_csv(gene_file, sep="\t", header=None, names=["gene_name", "num_isoforms"])
            target_df = target_df.dropna(subset=["num_isoforms"])
        except Exception as e:
            print(f"Error loading {gene_file}: {e}")
            continue
        
        # Run bootstrap Fisher test
        res = bootstrap_fisher_test(
            target_genes=target_df,
            background_genes=bg_df[["gene_name", "num_isoforms"]],
            label=f"{label}_Top{size}",
            threshold=splicing_threshold,
            n_iter=bootstrap_iterations,
            random_state=42
        )
        res["gene_set_size"] = size
        results.append(res)
    
    # Save results
    results_df = pd.DataFrame(results)
    results_df = results_df[["gene_set_size", "clade", "mean_odds_ratio", "se_odds_ratio", 
                             "mean_ci_low", "mean_ci_high", "significant_p_count", "total_iterations"]]
    results_df.to_csv(f"gene_set_sweep_results_{label}.tsv", sep="\t", index=False)
    
    # Skip plotting if no valid results
    if results_df["total_iterations"].eq(0).all() or results_df["mean_odds_ratio"].isna().all():
        print(f"Skipping plot for {label}: no valid bootstrap iterations")
        return
    
    # --- Plotting ---
    # Filter valid results for plotting
    valid_results = results_df[
        (results_df["total_iterations"] > 0) & 
        (results_df["mean_odds_ratio"].notna()) & 
        (results_df["mean_ci_low"].notna()) & 
        (results_df["mean_ci_high"].notna()) &
        (results_df["mean_ci_low"] <= results_df["mean_odds_ratio"]) & 
        (results_df["mean_ci_high"] >= results_df["mean_odds_ratio"])
    ]
    
    if valid_results.empty:
        print(f"Skipping plot for {label}: no valid data for plotting")
        return
    
    x_vals = valid_results["gene_set_size"]
    y_vals = valid_results["mean_odds_ratio"]
    
    # Clip yerr to ensure non-negative values
    yerr_lower = np.maximum(y_vals - valid_results["mean_ci_low"], 0)
    yerr_upper = np.maximum(valid_results["mean_ci_high"] - y_vals, 0)
    
    plt.figure(figsize=(8, 5))
    
    # Plot odds ratio with error bars
    plt.errorbar(
        x_vals,
        y_vals,
        yerr=[yerr_lower, yerr_upper],
        fmt='o-', capsize=4, label="Odds Ratio", color='steelblue'
    )
    
    # Confidence interval shading
    plt.fill_between(
        x_vals,
        valid_results["mean_ci_low"],
        valid_results["mean_ci_high"],
        color='lightblue',
        alpha=0.3,
        label="95% CI"
    )
    
    # Add significance stars
    for x, y, p_count, total in zip(x_vals, y_vals, valid_results["significant_p_count"], valid_results["total_iterations"]):
        if total > 0:  # Ensure there are valid iterations
            p_prop = p_count / total
            if p_prop >= 0.95:  # Equivalent to p < 0.05 in most iterations
                star = '***'
            elif p_prop >= 0.90:
                star = '**'
            elif p_prop >= 0.80:
                star = '*'
            else:
                star = ''
            if star:
                plt.text(x, y + 0.05, star, ha='center', va='bottom', fontsize=10, color='black')
    
    # Horizontal reference line
    plt.axhline(1, color='grey', linestyle='--', linewidth=1)
    
    # Labels and styling
    plt.xlabel("Gene Set Size", fontsize=12, weight='bold')
    plt.ylabel("Odds Ratio", fontsize=12, weight='bold')
    plt.title(f"Splicing Enrichment (10% Threshold) — {label}", fontsize=14, weight='bold')
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"gene_set_sweep_plot_{label}.png", dpi=300)
    plt.close()

# --- Run for both datasets ---
for label, gene_files in datasets.items():
    run_gene_set_sweep(label, gene_files)

print("Done. Results and plots saved.")
