import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
import sys

def compute_continuous_odds_ratios(adata, continuous_var='bmi_lv', celltype_col='Class'):
    celltypes = adata.obs[celltype_col].unique()

    odds_ratios = []
    conf_intervals = []

    # Loop through each cell type and fit logistic regression
    for celltype in celltypes:
        # Binary outcome: Does this cell type exist for each observation?
        adata.obs['is_celltype'] = (adata.obs[celltype_col] == celltype).astype(int)

        # Logistic regression with continuous variable
        X = sm.add_constant(adata.obs[continuous_var])  # Add constant term
        y = adata.obs['is_celltype']
        model = sm.Logit(y, X)
        result = model.fit(disp=False)

        # Get odds ratio (exponentiate the coefficients)
        odds_ratio = np.exp(result.params[continuous_var])
        ci_lower, ci_upper = np.exp(result.conf_int().loc[continuous_var])

        odds_ratios.append((celltype, odds_ratio))
        conf_intervals.append((celltype, ci_lower, ci_upper))

    odds_ratios_df = pd.DataFrame(odds_ratios, columns=['Cell Type', 'Odds Ratio'])
    conf_intervals_df = pd.DataFrame(conf_intervals, columns=['Cell Type', 'CI Lower', 'CI Upper'])

    return odds_ratios_df, conf_intervals_df


def plot_continuous_odds_ratios(odds_ratios_df, conf_intervals_df, continuous_var='bmi_lv', celltype_col='Class', save=None):
    # Merge odds ratios and confidence intervals
    plot_df = pd.merge(odds_ratios_df, conf_intervals_df, on='Cell Type')

    # Plot
    width = max(6, len(plot_df) * 0.3)
    # plt.figure(figsize=(12, 6))
    plt.figure(figsize=(width, 6))
    sns.barplot(x='Cell Type', y='Odds Ratio', data=plot_df, color='skyblue', ci=None)

    # Add error bars for confidence intervals
    plt.errorbar(x=plot_df['Cell Type'],
                 y=plot_df['Odds Ratio'],
                 yerr=[plot_df['Odds Ratio'] - plot_df['CI Lower'],
                       plot_df['CI Upper'] - plot_df['Odds Ratio']],
                 fmt='none', c='black', capsize=5)

    plt.axhline(1, color='red', linestyle='--')
    y_min = odds_ratios_df['Odds Ratio'].min() * 0.98
    y_max = odds_ratios_df['Odds Ratio'].max() * 1.02
    plt.ylim(top=y_max, bottom=y_min)
    plt.title(f'Odds Ratios {continuous_var} {celltype_col} with 95% confidence intervals')
    plt.xticks(rotation=90)
    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=300)
        print(f"Odds ratio plot saved to {save}")

    plt.show()


if __name__ == "__main__":
    # in_path = "/home/anna_y/data/write/all/AD427MR_50k/AD427MR_50k.h5ad" # For testing
    in_path = "/home/anna_y/data/write/all/AD427MR/AD427MR.h5ad"
    print("in_path:", in_path)

    adata = sc.read_h5ad(in_path)
    print(adata)

    continuous_var = 'bmi_lv'
    celltype_col = sys.argv[1] if len(sys.argv) > 1 else 'Class'

    # Filter out missing BMI values
    adata = adata[~adata.obs[continuous_var].isna()]

    # Compute continuous odds ratios
    odds_ratios_df, conf_intervals_df = compute_continuous_odds_ratios(adata, continuous_var=continuous_var, celltype_col=celltype_col)

    # Plot the odds ratios
    out_path = in_path.replace('.h5ad', f'_continuous_odds_ratios.png')
    out_path = in_path.replace('write', 'results/figures/cell_fraction').replace('.h5ad', f'_cell_fraction_{celltype_col}_odds_ratios.png')

    plot_continuous_odds_ratios(odds_ratios_df, conf_intervals_df, continuous_var, celltype_col, save=out_path)
