import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_volcano(class_deg_results, log2fc_cutoff=0.2, ax=None):
    # Basic volcano plot implementation for a single class
    ax.scatter(class_deg_results['log2FC'], -np.log10(class_deg_results['FDR']),
               c=(class_deg_results['FDR'] < 0.05) & (abs(class_deg_results['log2FC']) > log2fc_cutoff),
               cmap='coolwarm', edgecolor='k', s=20)
    ax.axhline(y=-np.log10(0.05), color='black', linestyle='--')
    ax.axvline(x=log2fc_cutoff, color='black', linestyle='--')
    ax.axvline(x=-log2fc_cutoff, color='black', linestyle='--')
    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-Log10 FDR')

def plot_volcano_grid(in_dir, classes, nrows, ncols, log2fc_cutoff=0.2, save=None):
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, figsize=(3*ncols, 3*nrows))
    axs = axs.flatten()  # Flatten to make indexing easier

    for i, class_name in enumerate(classes):
        filepath = os.path.join(in_dir, class_name, f'{class_name}.bmi_normalized.Clean.tsv')
        if not os.path.isfile(filepath):
            print(f'File not found: {filepath}')
            continue

        print(f'Reading file: {filepath}')
        class_deg_results = pd.read_csv(filepath, sep='\t')

        ax = axs[i]
        plot_volcano(class_deg_results, log2fc_cutoff, ax=ax)
        ax.set_title(class_name)

    plt.tight_layout()
    plt.savefig(save)
    # plt.show()

def plot_volcano_all_classes(log2fc_cutoff=0.2, ncols=5):
    in_dir = '/home/anna_y/results/deg_bmi_normalized/Class/'
    # in_paths = os.listdir('/home/anna_y/results/deg_bmi_lv/Class/*/*.Clean.tsv')
    # print(f'Files found: {in_paths}')
    classes = os.listdir(in_dir)
    print(f'Classes found: {classes}')

    nrows = len(classes) // ncols + 1
    print(f'Grid layout: {nrows} rows x {ncols} columns')
    plot_volcano_grid(in_dir, classes, nrows, ncols, log2fc_cutoff, save='figures/deg_bmi_normalized/volcano_all_classes.png')

def plot_volcano_subclasses(parent_class=None, log2fc_cutoff=1.0, ncols=5):
    in_dir = '/home/anna_y/results/deg_bmi_normalized/Subclass/'
    classes = os.listdir(in_dir)
    if parent_class:
        # only include subclasses of the parent class
        classes = [c for c in classes if parent_class in c]
    print(f'Classes found: {classes}')

    nrows = len(classes) // ncols + 1
    print(f'Grid layout: {nrows} rows x {ncols} columns')
    out_path = f'figures/deg_bmi_normalized/volcano_{parent_class}_subclasses.png'
    plot_volcano_grid(in_dir, classes, nrows, ncols, log2fc_cutoff, save=out_path)
    print(f'Volcano plots saved to {out_path}.')

def plot_volcano_subtypes(parent_class=None, log2fc_cutoff=0.2, ncols=5):
    in_dir = '/home/anna_y/results/deg_bmi_normalized/Subtype/'
    classes = os.listdir(in_dir)
    if parent_class:
        classes = [c for c in classes if parent_class in c]
    print(f'Classes found: {classes}')

    nrows = len(classes) // ncols + 1
    print(f'Grid layout: {nrows} rows x {ncols} columns')
    out_path = f'figures/deg_bmi_normalized/volcano_{parent_class}_subtypes.png'
    plot_volcano_grid(in_dir, classes, nrows, ncols, log2fc_cutoff, save=out_path)
    print(f'Volcano plots saved to {out_path}.')


if __name__ == "__main__":
    plot_volcano_all_classes(log2fc_cutoff=0.2, ncols=4)
    plot_volcano_subclasses(log2fc_cutoff=0.2, ncols=6)
    plot_volcano_subtypes(log2fc_cutoff=0.2, ncols=6)

    for parent_class in ['Exc', 'Inh']:
        plot_volcano_subclasses(parent_class, log2fc_cutoff=0.2)
        plot_volcano_subtypes(parent_class, log2fc_cutoff=0.2)

    for parent_class in ['Ast', 'Mic', 'Oli', 'OPC']:
        # plot_volcano_subclasses(parent_class, log2fc_cutoff=0.2)
        plot_volcano_subtypes(parent_class, log2fc_cutoff=0.2)
