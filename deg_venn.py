import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from deg_results import get_top_degs
import os

def get_gene_list_from_txt(txt_file):
    with open(txt_file, 'r') as f:
        degs = f.read().strip().split("\n")
        return degs

def get_gene_lists_dict(in_dir, classes):
    """return a dictionary of gene lists, where the keys are cell types and the values are sets of genes."""
    gene_lists_dict = {}
    for class_name in classes:
        pos_txt = in_dir + "/" + class_name + "/positive.txt"
        neg_txt = in_dir + "/" + class_name + "/negative.txt"
        degs = get_gene_list_from_txt(pos_txt) + get_gene_list_from_txt(neg_txt)
        gene_lists_dict[class_name] = set(degs)
    return gene_lists_dict

def venn_diagram(gene_lists_dict, save=None):
    # Calculate total number of genes
    total = len(set.union(*gene_lists_dict.values()))
    print(f"Total number of genes: {total}")

    # Create Venn diagram
    plt.figure(figsize=(8, 8))
    classes = tuple(gene_lists_dict.keys())
    gene_sets = tuple(gene_lists_dict.values())
    # print(f"Gene sets: {gene_sets}, classes: {classes} ")
    venn3(gene_sets, classes, subset_label_formatter=lambda count: f"{count}\n({count/total*100:.2f}%)")
    plt.title(f"Venn Diagram of DEGs in {classes} cells")
    if save:
        plt.savefig(save)

def venn_celltypes(classes, save=None):
    """Generate a Venn diagram showing the overlap of DEGs between different cell types."""
    gene_lists_dict_class = get_gene_lists_dict("/home/anna_y/data/results/deg_bmi_normalized_v1/Class", classes)
    common_genes = set.intersection(*gene_lists_dict_class.values())
    print(f"Number of common genes among {classes}: {len(common_genes)}\n{common_genes}")
    # classes = tuple(gene_lists_dict_class.keys())
    out_path = save if save else f"/home/anna_y/data/results/figures/deg_bmi_normalized_v1/venn_{'_'.join(classes)}.png"
    venn_diagram(gene_lists_dict=gene_lists_dict_class, save=out_path)
    print(f"Venn diagram saved to {out_path}")


if __name__ == "__main__":
    classes = ["Exc_50k", "Inh", "Oli"]
    venn_celltypes(classes)

    glial_cell_types = ["Ast", "Mic_Immune", "OPC"]
    venn_celltypes(glial_cell_types)
