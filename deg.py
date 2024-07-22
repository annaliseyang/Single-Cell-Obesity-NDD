# from init import *
import pandas as pd

def get_top_n_genes(group, csv_file, n_genes=5):
    deg_df = pd.read_csv(csv_file)
    top_genes = deg_df[group].tolist()
    top_genes = top_genes[:n_genes]
    return top_genes

if __name__ == "__main__":
    obesity_genes = get_top_n_genes('bmi_30+', 'rank_genes_groups_bmi.csv', n_genes=20)
    anti_obesity_genes = get_top_n_genes('bmi_<20', 'rank_genes_groups_bmi.csv', n_genes=20)
    early_AD_genes = get_top_n_genes('earlyAD', 'rank_genes_groups_AD.csv', n_genes=20)
    late_AD_genes = get_top_n_genes('lateAD', 'rank_genes_groups_AD.csv', n_genes=20)

    print(obesity_genes)
    print(anti_obesity_genes)
    print(early_AD_genes)
    print(late_AD_genes)
