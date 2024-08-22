import os
import numpy as np
import pandas as pd
import sys

def filter_deg_results(results_file, p_cutoff=0.05, fdr_cutoff=0.05, log2fc_cutoff=0.2):
    # Read the DEG results file
    results = pd.read_csv(results_file, sep='\t')

    # Filter the results based on the desired log2FC and p-value cutoff
    results_filtered = results[(abs(results['log2FC']) > log2fc_cutoff) &
                           (results['FDR'] < fdr_cutoff) &
                           (results['p_bmi_normalized'] < p_cutoff)]
    out_path = results_file.replace('.tsv', '.Filtered.tsv')
    results_filtered.to_csv(out_path, sep='\t', index=False)
    print(f'Filtered DEG results saved to: {out_path}')
    return out_path

def rank_deg_results(results_file, sort_by='log2FC'):
    results = pd.read_csv(results_file, sep='\t')
    # Sort the results by log2FC in descending order
    results_sorted = results.sort_values(by=sort_by, key=lambda col: abs(col), ascending=False)

    out_path = results_file.replace('.Clean.tsv', '.Ranked.tsv')
    results_sorted.to_csv(out_path, sep='\t', index=False)
    print(f'Sorted DEG results saved to: {out_path}')
    return out_path

def get_top_degs(results_file, n_top=20, positive=True):
    """
    Return the top n genes based on log2FC.
    positive (bool) indicates whether to return positive or negative log2FC genes.
    """
    results = pd.read_csv(results_file, sep='\t')
    if positive == True:
        results = results[results['log2FC'] > 0]
        top_genes = results.head(n_top)['gene'].tolist() if n_top else results['gene'].tolist()
    elif positive == False:
        results = results[results['log2FC'] < 0]
        top_genes = results.tail(n_top)['gene'].tolist() if n_top else results['gene'].tolist()
    else:
        print(f"'positive' parameter not set. Returning top {n_top} DEGs.")
        top_genes = results.head(n_top)['gene'].tolist() if n_top else results['gene'].tolist()
    return top_genes

if __name__ == "__main__":
    results_file = sys.argv[1]
    results_file_ranked = rank_deg_results(results_file)
    results_file_filtered = filter_deg_results(results_file_ranked, p_cutoff=0.05, fdr_cutoff=0.05, log2fc_cutoff=0)

    pos_degs = get_top_degs(results_file_filtered, n_top=None, positive=True)
    neg_degs = get_top_degs(results_file_filtered, n_top=None, positive=False)
    print(f'Top positive DEGs: {pos_degs}')
    print(f'Top negative DEGs: {neg_degs}')

    # save the top genes to a txt file
    pos_txt_path = os.path.dirname(results_file_filtered) + '/positive.txt'
    neg_txt_path = os.path.dirname(results_file_filtered) + '/negative.txt'
    with open(pos_txt_path, 'w') as f:
        f.write('\n'.join(pos_degs))
        print(f'Top positive genes saved to {pos_txt_path}')
    with open(neg_txt_path, 'w') as f:
        f.write('\n'.join(neg_degs))
        print(f'Top negative genes saved to {neg_txt_path}')
