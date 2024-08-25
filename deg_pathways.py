import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

def plot_heatmap(celltype, pos_path, neg_path, out_path, var='BestLogPInGroup', all=True):
    pos_results = pd.read_csv(pos_path)
    neg_results = pd.read_csv(neg_path)

    pos_go, pos_description = pos_results['GO'].tolist(), pos_results['Description'].tolist()
    neg_go, neg_description = neg_results['GO'].tolist(), neg_results['Description'].tolist()

    all_go = list(set(pos_go) | set(neg_go))
    common_go = list(set(pos_go) & set(neg_go))

    print(f"Positive GO terms: {len(pos_go)}")
    print(f"Negative GO terms: {len(neg_go)}")
    print(f"Total GO terms: {len(all_go)}")
    print(f"Common GO terms: {len(common_go)} {common_go}")

    go_terms = all_go if all else common_go
    descriptions = [pos_description[pos_go.index(go)] if go in pos_go else neg_description[neg_go.index(go)] for go in go_terms]
    labels = [f"{go}: {description}" for go, description in zip(go_terms, descriptions)]

    df = pd.DataFrame(0, index=labels, columns=['Positive', 'Negative'])

    pos_col = [-pos_results.at[pos_go.index(term), var] if term in pos_go else 0 for term in go_terms]
    neg_col = [-neg_results.at[neg_go.index(term), var] if term in neg_go else 0 for term in go_terms]

    df['Positive'] = pos_col
    df['Negative'] = neg_col
    # df['sum'] = df['Positive'] + df['Negative']
    df.sort_values(by=['Positive', 'Negative'], ascending=False, inplace=True)
    print(df)

    # create heatmap
    plt.figure(figsize=(12, 8))
    plt.title(f"{celltype} Pathway Enrichment Heatmap ({var})")
    sns.heatmap(df, cmap='Blues', annot=False, vmin=0, vmax=df.max().max())

    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.ylabel("Functional Ontology Term")
    # plt.legend(loc='right', title=var)
    plt.savefig(out_path, bbox_inches='tight')

    print(f"Heatmap saved to {out_path}")


if __name__ == "__main__":
    celltypes = ['Exc_50k', 'Ast', 'Mic_Immune', 'OPC', 'Oli', 'Vasc_Epithelia', 'Inh']
    # for dir in os.listdir("/home/anna_y/data/results/deg_bmi_normalized_v1_metascape/Class/"):
    for celltype in celltypes:
        print(celltype)
        # continue
        # celltype = dir.split('_')[0]
        # sign = dir.split('_')[1]
        pos_path = f"/home/anna_y/data/results/deg_bmi_normalized_v1_metascape/Class/{celltype}_positive/Enrichment_GO/_FINAL_GO.csv"
        neg_path = f"/home/anna_y/data/results/deg_bmi_normalized_v1_metascape/Class/{celltype}_negative/Enrichment_GO/_FINAL_GO.csv"
        out_path_common = f"/home/anna_y/data/results/figures/deg_bmi_normalized_v1_GO/Class/{celltype}_heatmap_common.png"
        out_path_all = f"/home/anna_y/data/results/figures/deg_bmi_normalized_v1_GO/Class/{celltype}_heatmap_all.png"

        plot_heatmap(celltype, pos_path, neg_path, out_path_common, all=False)
        plot_heatmap(celltype, pos_path, neg_path, out_path_all, all=True)
