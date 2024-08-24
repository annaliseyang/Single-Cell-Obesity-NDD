import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_heatmap(celltype, pos_path, neg_path, out_path, all=True):
    pos_results = pd.read_csv(pos_path)
    neg_results = pd.read_csv(neg_path)

    pos_go, pos_logp = pos_results['GO'].tolist(), pos_results['LogP'].tolist()
    neg_go, neg_logp = neg_results['GO'].tolist(), neg_results['LogP'].tolist()
    all_go = list(set(pos_go) | set(neg_go))
    common_go = list(set(pos_go) & set(neg_go))

    print(f"Positive GO terms: {len(pos_go)}")
    print(f"Negative GO terms: {len(neg_go)}")
    print(f"Total GO terms: {len(all_go)}")
    print(f"Common GO terms: {len(common_go)} {common_go}")

    go_terms = all_go if all else common_go
    df = pd.DataFrame(0, index=go_terms, columns=['Positive', 'Negative'])
    # df.loc[common_go, 'Positive'] = pos_results[pos_results['GO'].isin(common_go)]['LogP'].tolist()
    # df.loc[common_go, 'Negative'] = neg_results[neg_results['GO'].isin(common_go)]['LogP'].tolist()
    print(pos_results[pos_results['GO'].isin(go_terms)]['LogP'].tolist())
    print(pos_results[pos_results['GO'].isin(go_terms)]['BestLogPInGroup'].tolist())

    # pos_col = -pd.Series(pos_results[pos_results['GO'].isin(go_terms)]['LogP'].tolist(), index=go_terms)
    # neg_col = -pd.Series(neg_results[neg_results['GO'].isin(go_terms)]['LogP'].tolist(), index=go_terms)
    pos_col = [-pos_results.at[pos_go.index(term), 'BestLogPInGroup'] if term in pos_go else 0 for term in go_terms]
    neg_col = [-neg_results.at[neg_go.index(term), 'BestLogPInGroup'] if term in neg_go else 0 for term in go_terms]

    df['Positive'] = pos_col
    df['Negative'] = neg_col
    df.sort_values(by=['Positive', 'Negative'], ascending=False, inplace=True)
    print(df)

    # create heatmap
    plt.figure(figsize=(12, 8))
    plt.title(f"{celltype} GO Enrichment Heatmap")
    sns.heatmap(df, annot=False, vmin=0, vmax=df.max().max())

    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.xlabel("Condition")
    plt.ylabel("GO Terms")
    plt.savefig(out_path, bbox_inches='tight')

    print(f"Heatmap saved to {out_path}")


if __name__ == "__main__":
    pos_path = "/home/anna_y/data/results/deg_bmi_normalized_v1_metascape/Class/Exc_50k_positive/Enrichment_GO/_FINAL_GO.csv"
    neg_path = "/home/anna_y/data/results/deg_bmi_normalized_v1_metascape/Class/Exc_50k_negative/Enrichment_GO/_FINAL_GO.csv"
    out_path_common = "/home/anna_y/data/results/figures/deg_bmi_normalized_v1_GO/Class/Exc_50k_heatmap_common.png"
    out_path_all = "/home/anna_y/data/results/figures/deg_bmi_normalized_v1_GO/Class/Exc_50k_heatmap_all.png"

    plot_heatmap("Exc", pos_path, neg_path, out_path_common, all=False)
    plot_heatmap("Exc", pos_path, neg_path, out_path_all, all=True)
