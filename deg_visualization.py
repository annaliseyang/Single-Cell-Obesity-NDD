import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_volcano(class_deg_results, log2fc_cutoff=1.0, ax=None):
    # Basic volcano plot implementation for a single class
    ax.scatter(class_deg_results['log2FC'], -np.log10(class_deg_results['FDR']),
               c=(class_deg_results['FDR'] < 0.05) & (abs(class_deg_results['log2FC']) > log2fc_cutoff),
               cmap='coolwarm', edgecolor='k', s=20)
    ax.axhline(y=-np.log10(0.05), color='black', linestyle='--')
    ax.axvline(x=log2fc_cutoff, color='black', linestyle='--')
    ax.axvline(x=-log2fc_cutoff, color='black', linestyle='--')
    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-Log10 FDR')

def plot_volcano_all_classes(log2fc_cutoff=1.0):
    in_dir = '/home/anna_y/results/deg_bmi_lv/Class/'
    # in_paths = os.listdir('/home/anna_y/results/deg_bmi_lv/Class/*/*.Clean.tsv')
    # print(f'Files found: {in_paths}')
    classes = os.listdir(in_dir)
    print(f'Classes found: {classes}')

    nrows = 3
    ncols = (len(classes) + nrows - 1) // nrows  # This ensures enough columns
    print(f'Grid layout: {nrows} rows x {ncols} columns')

    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, figsize=(15, 10))
    axs = axs.flatten()  # Flatten to make indexing easier

    for i, class_name in enumerate(classes):
        filepath = os.path.join(in_dir, class_name, f'{class_name}.bmi_lv.Clean.tsv')
        if not os.path.isfile(filepath):
            print(f'File not found: {filepath}')
            # i -= 1 # Skip the plot if the file does not exist
            continue

        print(f'Reading file: {filepath}')
        class_deg_results = pd.read_csv(filepath, sep='\t')

        ax = axs[i]
        plot_volcano(class_deg_results, log2fc_cutoff, ax=ax)
        ax.set_title(class_name)

    plt.tight_layout()
    plt.savefig('figures/deg/volcano_all_classes.png')
    # plt.show()


if __name__ == "__main__":
    plot_volcano_all_classes(log2fc_cutoff=1.0)
