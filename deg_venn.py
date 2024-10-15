import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from deg_results import get_top_degs
import math
import os

def get_gene_list_from_txt(txt_file):
    """return a list of top DEGs from a txt file."""
    with open(txt_file, 'r') as f:
        degs = f.read().strip().split("\n")
        return degs

def get_gene_lists_dict(in_dir, classes):
    """return a dictionary of gene lists, where the keys are cell types and the values are sets of top DEGs."""
    gene_lists_dict_all = {}
    gene_lists_dict_pos = {}
    gene_lists_dict_neg = {}
    for class_name in classes:
        pos_txt = in_dir + "/" + class_name + "/positive.txt"
        neg_txt = in_dir + "/" + class_name + "/negative.txt"
        pos_degs = get_gene_list_from_txt(pos_txt)
        neg_degs = get_gene_list_from_txt(neg_txt)

        gene_lists_dict_all[class_name] = set(pos_degs + neg_degs)
        gene_lists_dict_pos[class_name] = set(pos_degs)
        gene_lists_dict_neg[class_name] = set(neg_degs)
    return gene_lists_dict_all, gene_lists_dict_pos, gene_lists_dict_neg

def venn_diagram(gene_lists_dict, save=None, sign:int = 0, ax=None):
    # Calculate total number of genes
    total = len(set.union(*gene_lists_dict.values()))
    print(f"Total number of genes: {total}")

    # Create Venn diagram
    plt.figure(figsize=(8, 8))
    classes = tuple(gene_lists_dict.keys())
    gene_sets = tuple(gene_lists_dict.values())
    # print(f"Gene sets: {gene_sets}, classes: {classes} ")
    venn3(gene_sets, classes, subset_label_formatter=lambda count: f"{count}\n({count/total*100:.2f}%)", ax=ax)
    sign_str = ["All", "Positive", "Negative"][sign]
    plt.title(f"Venn Diagram of {sign_str} DEGs in {classes} cells")
    if save:
        plt.savefig(save)

def venn_diagram(gene_lists_dict, save=None, sign: int = 0, ax=None):
    # Calculate total number of genes
    total = len(set.union(*gene_lists_dict.values()))
    print(f"Total number of genes: {total}")

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    else:
        fig = ax.get_figure()

    classes = tuple(gene_lists_dict.keys())
    gene_sets = tuple(gene_lists_dict.values())

    # plot the Venn diagram
    venn3(gene_sets, classes, subset_label_formatter=lambda count: f"{count}\n({count/total*100:.2f}%)", ax=ax)

    sign_str = ["All", "Positive", "Negative"][sign]
    ax.set_title(f"Venn Diagram of {sign_str} DEGs in {classes} cells")

    if save:
        fig.savefig(save)
        print(f"Venn diagram saved to {save}")

    plt.tight_layout()
    return fig, ax

def venn_celltypes(classes, save=None, sign:int = 0, ax=None):
    """
    Generate a Venn diagram showing the overlap of DEGs between different cell types.
    Selectively plot DEGs by sign: 0 for all, 1 for positive, -1 for negative.
    """
    gene_lists_dict_class = get_gene_lists_dict("/home/anna_y/data/results/deg_bmi_normalized_v1/Class", classes)[sign]
    common_genes = set.intersection(*gene_lists_dict_class.values())
    # classes = tuple(gene_lists_dict_class.keys())
    sign_str = ["All", "Positive", "Negative"][sign]
    print(f"Number of common {sign_str} DEGs among {classes}: {len(common_genes)}\n{common_genes}")

    # out_path = save if save else f"/home/anna_y/data/results/figures/deg_bmi_normalized_v1/venn_{sign_str}_degs_{'_'.join(classes)}.png"
    venn_diagram(gene_lists_dict=gene_lists_dict_class, save=save, sign=sign, ax=ax)
    print(f"Venn diagram saved to {save}")

def generate_combinations(classes, n):
    """Generate all combinations of n classes chosen from the given list."""
    if n == 0:
        return [[]]
    if n == len(classes):
        return [classes]

    combos = []
    current_class = classes[0]
    rec_result = generate_combinations(classes[1:], n-1)
    for combo in rec_result:
        combos.append([current_class] + combo)
    rec_result_without = generate_combinations(classes[1:], n)
    combos += rec_result_without
    return combos

def venn_celltypes_combinations(classes, n=3, n_cols=4, save=None, sign:int = 0):
    """Generate a Venn diagram showing the overlap of DEGs between different combinations of cell types."""
    combos = generate_combinations(classes, n)
    n_rows = math.ceil(len(combos) / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 6*n_rows))
    axes = axes.flatten() if len(combos) > 1 else [axes]

    for i, combo in enumerate(combos):
        print(f"Venn diagram for {combo}")
        venn_celltypes(combo, sign=sign, ax=axes[i])

    plt.tight_layout()
    if save:
        plt.savefig(save)
        print(f"Venn diagrams saved to {save}")
    else:
        out_path = f"/home/anna_y/data/results/figures/deg_bmi_normalized_v1/venn_{sign}_degs_{n}_celltypes.png"
        plt.savefig(out_path)
        print(f"Venn diagram saved to", out_path)

if __name__ == "__main__":
    classes = ["Exc_50k", "Inh", "Oli"]
    all_classes = ["Exc_50k", "Inh", "Oli", "Ast", "Mic_Immune", "OPC"]
    # venn_celltypes(classes)
    triples = generate_combinations(all_classes, 3)
    print(f"Number of combinations: {len(triples)}")
    print(f"Triples: {triples}")

    signs = {"All": 0, "Positive": 1, "Negative": -1}
    for sign_str, sign in signs.items():
        venn_celltypes_combinations(all_classes, n=3, n_cols=5, save=f"/home/anna_y/data/results/figures/deg_bmi_normalized_v1/venn_{sign_str}_degs_{'_'.join(all_classes)}.png", sign=sign)

    # glial_cell_types = ["Ast", "Mic_Immune", "OPC"]
    # venn_celltypes(glial_cell_types, sign=0)
    # venn_celltypes(glial_cell_types, sign=1)
    # venn_celltypes(glial_cell_types, sign=-1)
