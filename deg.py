# from init import *
import numpy as np
import pandas as pd
import scanpy as sc
import json

def rank_genes_groups(adata, groupby='AD_states', type=None):
    """
    Rank genes within each group and save the results to a CSV file.

    Parameters:
    adata (AnnData): The input AnnData object.
    groupby (str): The column name in adata.obs that defines the groups for ranking.
    """
    adata = adata.copy()
    adata.X = np.log1p(adata.X)

    sc.tl.rank_genes_groups(adata, groupby=groupby)
    print(adata.uns['rank_genes_groups'], flush=True)
    result = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    parameters = adata.uns['rank_genes_groups']['params']

    # save the ranked genes and parameters to CSV files and JSON file respectively
    result.to_csv(f"rank_genes_groups_{type}_by_{groupby}.csv", index=False)
    json.dump(parameters, open(f"rank_genes_groups_{type}_by_{groupby}_parameters.json", 'w'))
    # sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, groupby=groupby, standard_scale='var', save=f'_{type}_{groupby}.png')

    return result

def get_top_n_genes(group: str, csv_file, n_genes=5):
    """
    Get the top n genes for a given group from the provided CSV file.

    Parameters:
    group (str): The group to get the top genes for.
    csv_file (str): The CSV file containing the DEG results.
    n_genes (int): The number of top genes to return.

    Returns:
    list: A list of the top n genes for the given group.
    """
    deg_df = pd.read_csv(csv_file)
    top_genes = deg_df[group].tolist()
    top_genes = top_genes[:n_genes]
    return top_genes

if __name__ == "__main__":
    obesity_genes = get_top_n_genes('bmi_30+', 'rank_genes_groups_bmi.csv', n_genes=20)
    anti_obesity_genes = get_top_n_genes('bmi_<20', 'rank_genes_groups_bmi.csv', n_genes=20)
    early_AD_genes = get_top_n_genes('earlyAD', 'rank_genes_groups_AD.csv', n_genes=20)
    late_AD_genes = get_top_n_genes('lateAD', 'rank_genes_groups_AD.csv', n_genes=20)

    print(f"{obesity_genes=}")
    print(f"{anti_obesity_genes=}")
    print(f"{early_AD_genes=}")
    print(f"{late_AD_genes=}")
