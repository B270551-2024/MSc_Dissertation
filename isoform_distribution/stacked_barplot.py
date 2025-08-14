## Isoform Distribution Analysis Across Gene Sets
## Compares isoform type distributions between different gene sets using Chi-squared tests
## Creates stacked bar plots showing relative proportions of isoform types
## Usage: Ensure all gene set files are in working directory, run script
## Created by B270551
##---------------------------------------------------------##

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import numpy as np

np.random.seed(42)

available_fonts = [f.name for f in fm.fontManager.ttflist]
font_family = 'Arial' if 'Arial' in available_fonts else 'DejaVu Sans'
if font_family == 'DejaVu Sans':
    print("Arial font not found, using DejaVu Sans as fallback.")

sns.set_style("whitegrid")
plt.rcParams.update({
    'font.family': font_family, 'font.size': 14, 'axes.titlesize': 22, 'axes.titleweight': 'bold',
    'axes.labelsize': 18, 'axes.labelweight': 'bold', 'legend.fontsize': 14, 'legend.title_fontsize': 15,
    'xtick.labelsize': 13, 'ytick.labelsize': 13, 'axes.grid': True, 'grid.color': '#E5E5E5',
    'grid.alpha': 0.7, 'grid.linestyle': '-', 'grid.linewidth': 0.5, 'axes.grid.axis': 'y',
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.spines.left': True,
    'axes.spines.bottom': True, 'axes.edgecolor': '#333333', 'axes.linewidth': 1.2
})

file_paths = {
    'All Protein Coding': 'ensembl_isoforms_v4.tsv',
    'Positively Selected, Spliced': 'high_selected_genes_ALL.tsv',
    'CRISPR Antiviral, Spliced': 'top_10_percent_proviral_genes.tsv',
    'CRISPR Proviral, Spliced': 'top_10_percent_antiviral_genes.tsv',
    'ISGs: Core Vert-Mamm': 'core_ISGs_vert-mamm.tsv',
    'ISGs: Antiviral': 'antiviral_genes_all.tsv',
    'GWAS: Immune, Spliced': 'genes_overlapping_immune_gwas_annotated_top_10_percent.txt',
    'GWAS: Immune, Antiviral, Spliced': 'top_10_percent_antiviral_gwas_genes.tsv'
}

isoform_cols = [
    'num_protein_coding_isoforms', 'num_retained_intron_isoforms', 'num_nonsense_mediated_decay_isoforms',
    'num_protein_coding_CDS_not_defined_isoforms', 'num_non_coding_isoforms'
]

def load_and_validate_data(file_paths, isoform_cols):
    isoform_summary_file = file_paths['All Protein Coding']
    if not os.path.exists(isoform_summary_file):
        raise FileNotFoundError(f"Isoform summary file {isoform_summary_file} not found.")
    
    try:
        isoform_df = pd.read_csv(isoform_summary_file, sep='\t')
        if not all(col in isoform_df.columns for col in ['gene_name', 'gene_biotype'] + isoform_cols):
            raise ValueError("Required columns missing in isoform summary file.")
        isoform_df['gene_name'] = isoform_df['gene_name'].str.strip().str.upper()
        print(f"Loaded {len(isoform_df)} records from {isoform_summary_file}")
    except Exception as e:
        raise ValueError(f"Error reading isoform summary file: {e}")
    
    isoform_df = isoform_df[isoform_df['gene_biotype'] == 'protein_coding']
    
    gene_sets = {}
    for name, file in file_paths.items():
        if name == 'All Protein Coding':
            gene_sets[name] = set(isoform_df['gene_name'])
            continue
        if not os.path.exists(file):
            raise FileNotFoundError(f"Gene list file {file} not found.")
        try:
            gene_df = pd.read_csv(file, sep='\t', header=None, names=['gene_name'])
            gene_df['gene_name'] = gene_df['gene_name'].str.strip().str.upper()
            if gene_df['gene_name'].isna().any():
                raise ValueError(f"Missing gene names detected in {file}.")
            gene_sets[name] = set(gene_df['gene_name'])
            print(f"Loaded {len(gene_sets[name])} genes from {file}")
        except Exception as e:
            raise ValueError(f"Error reading {file}: {e}")
    
    for name, genes in gene_sets.items():
        missing = genes - set(isoform_df['gene_name'])
        if missing:
            print(f"Warning: {len(missing)} genes from {name} not found in isoform data.")
    
    return isoform_df, gene_sets

def compute_isoform_distribution(isoform_df, gene_set, isoform_cols):
    filtered_df = isoform_df[isoform_df['gene_name'].isin(gene_set)]
    if filtered_df.empty:
        print(f"No matching genes found for the given set.")
        return pd.Series(0, index=isoform_cols)
    isoform_sums = filtered_df[isoform_cols].sum()
    return isoform_sums[isoform_sums > 0]

def downsample_background(isoform_df, background_genes, target_count):
    if target_count <= 0:
        return pd.DataFrame(columns=isoform_cols)
    
    isoform_df['total_isoforms'] = isoform_df[isoform_cols].sum(axis=1)
    background_df = isoform_df[isoform_df['gene_name'].isin(background_genes)].copy()
    background_df = background_df.sample(frac=1, random_state=42).reset_index(drop=True)
    
    selected_genes, current_count = [], 0
    for _, row in background_df.iterrows():
        if current_count + row['total_isoforms'] <= target_count or len(selected_genes) == 0:
            selected_genes.append(row['gene_name'])
            current_count += row['total_isoforms']
        if current_count >= target_count:
            break
    
    downsampled_df = isoform_df[isoform_df['gene_name'].isin(selected_genes)]
    isoform_sums = downsampled_df[isoform_cols].sum()
    return isoform_sums[isoform_sums > 0]

def perform_chi_squared_test(isoform_df, observed_sums, background_genes):
    common_indices = observed_sums.index
    target_count = observed_sums.sum()
    
    background_sums = downsample_background(isoform_df, background_genes, target_count)
    if background_sums.empty:
        return None, None, None
    
    common_indices = observed_sums.index.intersection(background_sums.index)
    observed = observed_sums[common_indices].values
    expected = background_sums[common_indices].values * (observed.sum() / background_sums[common_indices].sum())
    
    if expected.sum() == 0 or observed.sum() == 0:
        return None, None, None
    chi2, p_value = stats.chisquare(observed, expected)
    n = observed.sum()
    k = len(common_indices)
    cramers_v = np.sqrt(chi2 / (n * (k - 1))) if (n * (k - 1)) > 0 else None
    return chi2, p_value, cramers_v

def plot_stacked_bar(distributions, output_file, chi2_results):
    dist_df = pd.DataFrame(distributions).T
    dist_df = dist_df.div(dist_df.sum(axis=1), axis=0) * 100
    dist_df = dist_df.fillna(0)
    
    labels = dist_df.columns.str.replace("num_", "").str.replace("_isoforms", "").str.replace("_", " ").str.title()
    dist_df.columns = labels
    
    fig, ax = plt.subplots(figsize=(20, 10))
    colours = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#6A994E']
    
    bars = dist_df.plot(kind='bar', stacked=True, ax=ax, color=colours, width=0.7, edgecolor='white', linewidth=0.5)
    
    ax.set_title("Isoform Type Distribution Across Gene Sets", fontsize=24, weight='bold', pad=35, loc='centre')
    ax.set_xlabel("Gene Set", fontsize=20, weight='bold', labelpad=20)
    ax.set_ylabel("Percentage (%)", fontsize=20, weight='bold', labelpad=20)
    
    ax.set_xticklabels(dist_df.index, rotation=45, ha='right', fontsize=14)
    ax.tick_params(axis='y', labelsize=14)
    
    ax.legend(title="Isoform Types", loc='centre left', bbox_to_anchor=(1.02, 0.5), 
             fontsize=14, title_fontsize=16, frameon=True, fancybox=True, shadow=True, framealpha=0.9)
    
    for i, (gene_set, (_, p_value, _)) in enumerate(chi2_results.items()):
        if p_value is not None and gene_set != 'All Protein Coding':
            if p_value < 0.001:
                p_text = 'p < 0.001'
            else:
                p_text = f'p = {p_value:.3f}'
            
            if p_value < 0.001:
                p_text += ' ***'
            elif p_value < 0.01:
                p_text += ' **'
            elif p_value < 0.05:
                p_text += ' *'
            
            ax.text(i, -8, p_text, ha='centre', va='top', fontsize=12, weight='bold', rotation=0, colour='#333333')
    
    ax.set_ylim(-15, 105)
    ax.yaxis.grid(True, linestyle='-', alpha=0.3)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.20)
    
    plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight', facecolour='white', edgecolour='none')
    plt.savefig(f"{output_file}.pdf", bbox_inches='tight', facecolour='white', edgecolour='none')
    plt.show()

def save_summary_report(distributions, chi2_results, output_file):
    with open(f"{output_file}.txt", "w") as f:
        f.write("Isoform Distribution Analysis Summary\n")
        f.write("="*50 + "\n\n")
        
        for gene_set, dist in distributions.items():
            f.write(f"\nGene Set: {gene_set}\n")
            f.write("-" * (len(gene_set) + 11) + "\n")
            f.write(f"Total Isoforms: {int(dist.sum()):,}\n")
            f.write("Distribution:\n")
            for isoform_type, count in dist.items():
                percentage = (count / dist.sum()) * 100
                f.write(f"  {isoform_type}: {int(count):,} ({percentage:.1f}%)\n")
            
            chi2, p_value, cramers_v = chi2_results[gene_set]
            if chi2 is not None:
                significance = ""
                if p_value < 0.001:
                    significance = " (***)"
                elif p_value < 0.01:
                    significance = " (**)"
                elif p_value < 0.05:
                    significance = " (*)"
                
                f.write(f"Statistical Test Results:\n")
                f.write(f"  Chi-squared: {chi2:.2f}\n")
                f.write(f"  p-value: {p_value:.4f}{significance}\n")
                f.write(f"  Cramér's V: {cramers_v:.4f}\n")
            else:
                f.write("Statistical Test: Not applicable (reference set)\n")
            f.write("\n")
    
    dist_df = pd.DataFrame(distributions).T
    dist_df.index.name = 'Gene_Set'
    dist_df.to_csv(f"{output_file}_distributions.csv")
    
    dist_pct_df = dist_df.div(dist_df.sum(axis=1), axis=0) * 100
    dist_pct_df.to_csv(f"{output_file}_distributions_percentage.csv")

try:
    isoform_df, gene_sets = load_and_validate_data(file_paths, isoform_cols)
    
    distributions = {}
    for name, genes in gene_sets.items():
        distributions[name] = compute_isoform_distribution(isoform_df, genes, isoform_cols)
    
    background_genes = gene_sets['All Protein Coding']
    chi2_results = {}
    for name in distributions:
        if name == 'All Protein Coding':
            chi2_results[name] = (None, None, None)
            continue
        chi2, p_value, cramers_v = perform_chi_squared_test(isoform_df, distributions[name], background_genes)
        chi2_results[name] = (chi2, p_value, cramers_v)
    
    p_values = [p for _, p, _ in chi2_results.values() if p is not None]
    if p_values:
        adjusted_p = multipletests(p_values, method='fdr_bh')[1]
        i = 0
        for name in chi2_results:
            if chi2_results[name][1] is not None:
                chi2_results[name] = (chi2_results[name][0], adjusted_p[i], chi2_results[name][2])
                i += 1
    
    plot_stacked_bar(distributions, "isoform_distribution_comparison", chi2_results)
    save_summary_report(distributions, chi2_results, "isoform_analysis_summary")
    
    print("\nAnalysis completed successfully!")
    print("Files generated:")
    print("  - isoform_distribution_comparison.png (high-resolution plot)")
    print("  - isoform_distribution_comparison.pdf (vector format)")
    print("  - isoform_analysis_summary.txt (detailed report)")
    print("  - isoform_analysis_summary_distributions.csv (raw counts)")
    print("  - isoform_analysis_summary_distributions_percentage.csv (percentages)")
    
except Exception as e:
    print(f"Error: {e}")
