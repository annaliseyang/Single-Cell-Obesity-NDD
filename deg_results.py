import os
import numpy as np
import pandas as pd

def filter_deg_results(results_file, p_cutoff=0.05, fdr_cutoff=0.05, log2fc_cutoff=1.0):
    # Read the DEG results file
    results = pd.read_csv(results_file, sep='\t')

    # # Calculate the negative log10 of the FDR
    # results['-log10_FDR'] = -np.log10(results['FDR'])

    # Filter the results based on the desired log2FC cutoff
    results_filtered = results[abs(results['log2FC']) > log2fc_cutoff & results['FDR'] < fdr_cutoff & results['p_bmi_lv'] < p_cutoff]

    # Sort the results by log2FC in descending order
    results_sorted = results_filtered.sort_values(by='log2FC', ascending=False)

    # Save the filtered and sorted results to a new file
    # out_path = os.path.join(out_dir, 'deg_results_filtered.tsv')
    # out_path = results_file.replace('.tsv', '_filtered.tsv')
    # results_sorted.to_csv(os.path.join(out_dir, 'deg_results_filtered'), sep='\t', index=False)
