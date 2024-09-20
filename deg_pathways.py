import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

def combine_go_results(results_dir):
    """Save a single CSV file with all GO enrichment results."""
    sub_dirs = os.listdir(results_dir)
    final_go_files = [f"{results_dir}/{subdir}/Enrichment_GO/_FINAL_GO.csv" for subdir in sub_dirs]

    out = pd.DataFrame()
    # for file in final_go_files:
    for subdir in sub_dirs:
        results_file = f"{results_dir}/{subdir}/Enrichment_GO/_FINAL_GO.csv"
        if not os.path.exists(results_file):
            print(f"File not found: {results_file}")
            continue
        df = pd.read_csv(os.path.join(results_dir, results_file))
        df['gene_list'] = subdir
        df['file'] = results_file
        out = pd.concat([out, df])

    print(out.head())
    out.to_csv(results_dir + '/FINAL_GO_ALL.csv', index=False)
    return out


def get_description(go_term, results_path):
    df = pd.read_csv(results_path)
    description = df.loc[df['GO'] == go_term, 'Description'].values[0]
    return description


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
    # print(df)

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


def degs_in_pathway(go_term, results_path, out_path):
    """Return a list of genes in the given pathway."""
    # df = pd.read_excel(io=results_path, sheet_name='Annotation')
    print(go_term)
    df = pd.read_csv(results_path)
    # print(df.head(5))

    df = df[df['GO'] == go_term]
    print(df['Hits'].values[0])
    genes = df['Hits'].values[0].split('|')
    print(f"Genes in {go_term}: {genes}")
    return genes


def genes_by_celltype_heatmap(genes:list, label:str, results_dir, out_path, var='log2FC'):
    # genes = degs_in_pathway(go_term, results_path, out_path)
    print(f"input genes: {genes}\nlabel: {label}")
    celltypes = os.listdir(results_dir)
    print(celltypes)

    # dataframe to store logFC values for each gene (row) and celltype (column)
    df = pd.DataFrame(0, index=genes, columns=celltypes)
    for celltype in celltypes:
        deg_results_path = os.path.join(results_dir, celltype, f"{celltype}.bmi_normalized.Clean.tsv")
        # print(f"Genes in {go_term}: {genes}")
        deg_results = pd.read_csv(deg_results_path, sep='\t')
        # log_fc_df = df[df['gene'].isin(genes)]
        deg_results = deg_results[deg_results['gene'].isin(genes)]
        print(deg_results.head(5)) # print first 5 rows of the dataframe
        # log_fc = deg_results.loc['log2FC'].tolist()
        genes_in_celltype = deg_results['gene'].tolist()
        df.loc[genes_in_celltype, celltype] = deg_results[var].tolist()

    # create heatmap
    fig, ax = plt.subplots(figsize=(12, 8))
    title = label if label else f"{go_term} Genes Heatmap ({var})"
    plt.title(title)
    sns.heatmap(df, cmap='coolwarm', annot=None, ax=ax)

    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.ylabel("Gene")
    plt.tight_layout()
    plt.savefig(out_path)

def genes_in_pathway_heatmap(go_term, results_path, out_path, var='log2FC'):
    description = get_description(go_term, "/home/anna_y/data/results/deg_bmi_normalized_v1_metascape/Class/FINAL_GO_ALL.csv")
    print(f"Description: {description}")
    # results_path = f"/home/anna_y/data/results/deg_bmi_normalized_v1_metascape/Class/{celltype}_positive/Enrichment_GO/_FINAL_GO.csv"
    # var = 'log2FC'
    heatmap_title = f"{go_term} - {description} ({var})"
    # out_path = f"/home/anna_y/data/results/figures/deg_bmi_normalized_v1_metascape/Class/{celltype}_positive_{go_term}.png"
    genes = degs_in_pathway(go_term, results_path, out_path)
    genes_by_celltype_heatmap(genes, heatmap_title, "/home/anna_y/data/results/deg_bmi_normalized_v1/Class/", out_path)
    print(f"Genes heatmap saved to {out_path}")

def get_top_go_terms(results_path):
    df = pd.read_csv(results_path)
    df.sort_values(by='BestLogPInGroup', ascending=False, inplace=True)
    # print(df.head(10))
    go_list = df['GO'].tolist()
    return go_list

if __name__ == "__main__":
    celltypes = ['Exc_50k', 'Ast', 'Mic_Immune', 'OPC', 'Oli', 'Vasc_Epithelia', 'Inh']
    # for dir in os.listdir("/home/anna_y/data/results/deg_bmi_normalized_v1_metascape/Class/"):
    for celltype in celltypes:
        # continue
        print(f"\nProcessing {celltype}...")
        # continue
        # celltype = dir.split('_')[0]
        # sign = dir.split('_')[1]
        pos_path = f"/home/anna_y/data/results/deg_bmi_normalized_v1_metascape/Class/{celltype}_positive/Enrichment_GO/_FINAL_GO.csv"
        neg_path = f"/home/anna_y/data/results/deg_bmi_normalized_v1_metascape/Class/{celltype}_negative/Enrichment_GO/_FINAL_GO.csv"

        # get top GO terms
        pos_go_terms = get_top_go_terms(pos_path)
        neg_go_terms = get_top_go_terms(neg_path)
        common_go_terms = set(pos_go_terms) & set(neg_go_terms)
        print(f"Common GO terms: {common_go_terms}...")

        # plot heatmaps of top pathways
        out_path_common = f"/home/anna_y/data/results/figures/deg_bmi_normalized_v1_metascape/Class/{celltype}_heatmap_common.png"
        out_path_all = f"/home/anna_y/data/results/figures/deg_bmi_normalized_v1_metascape/Class/{celltype}_heatmap_all.png"
        plot_heatmap(celltype, pos_path, neg_path, out_path_common, all=False)
        plot_heatmap(celltype, pos_path, neg_path, out_path_all, all=True)

        # plot heatmap of top genes in each pathway
        n_top = 3
        for go_term in pos_go_terms[:n_top]:
            out_path_genes = f"/home/anna_y/data/results/figures/deg_bmi_normalized_v1_metascape/Class/{celltype}_positive_{go_term}.png"
            genes_in_pathway_heatmap(go_term, pos_path, out_path_genes)
            # genes_in_pathway_heatmap(go_term, neg_path, out_path_common)
        for go_term in neg_go_terms[:n_top]:
            out_path_genes = f"/home/anna_y/data/results/figures/deg_bmi_normalized_v1_metascape/Class/{celltype}_negative_{go_term}.png"
            # genes_in_pathway_heatmap(go_term, pos_path, out_path_common)
            genes_in_pathway_heatmap(go_term, neg_path, out_path_genes)
        for go_term in common_go_terms:
            out_path_genes = f"/home/anna_y/data/results/figures/deg_bmi_normalized_v1_metascape/Class/{celltype}_common_{go_term}.png"
            genes_in_pathway_heatmap(go_term, pos_path, out_path_genes)

        # break


    # celltype = 'Exc_50k'
    # pathways = ['GO:0043525', 'GO:0042110', 'GO:0046649', 'hsa05224'] # top pathways for this cell type

    # for go_term in pathways:
    #     description = get_description(go_term, "/home/anna_y/data/results/deg_bmi_normalized_v1_metascape/Class/FINAL_GO_ALL.csv")
    #     print(f"Description: {description}")
    #     results_path = f"/home/anna_y/data/results/deg_bmi_normalized_v1_metascape/Class/{celltype}_positive/Enrichment_GO/_FINAL_GO.csv"
    #     var = 'log2FC'
    #     heatmap_title = f"{go_term}: {description} ({var})"
    #     out_path = f"/home/anna_y/data/results/figures/deg_bmi_normalized_v1_metascape/Class/{celltype}_positive_{go_term}.png"
    #     genes = degs_in_pathway(go_term, results_path, out_path)
    #     genes_by_celltype_heatmap(genes, heatmap_title, "/home/anna_y/data/results/deg_bmi_normalized_v1/Class/", out_path)
    #     print(f"Genes heatmap saved to {out_path}")

    # all_go_terms = combine_go_results("/home/anna_y/data/results/deg_bmi_normalized_v1_metascape/Class")
