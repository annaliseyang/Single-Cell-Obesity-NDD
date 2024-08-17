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

def rank_deg_results(results_file, sort_by='log2FC'):
    results = pd.read_csv(results_file, sep='\t')
    # Sort the results by log2FC in descending order
    results_sorted = results.sort_values(by=sort_by, key=lambda col: abs(col), ascending=False)

    out_path = results_file.replace('.Clean.tsv', '.Ranked.tsv')
    results_sorted.to_csv(out_path, sep='\t', index=False)
    print(f'Sorted DEG results saved to: {out_path}')

if __name__ == "__main__":
    results_file = sys.argv[1]
    rank_deg_results(results_file)
    results_file_ranked = results_file.replace('.Clean.tsv', '.Ranked.tsv')
    filter_deg_results(results_file_ranked)
