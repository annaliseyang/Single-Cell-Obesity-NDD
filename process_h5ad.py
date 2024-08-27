import sys
import scanpy as sc
import pandas as pd

def define_bmi_groups(adata):
    if 'bmi_groups_copy' in adata.obs:
        print('bmi groups already defined, skipping...')
        return adata

    bmi_groups_range = {
        'bmi_0-20': 0,
        'bmi_20-25': 20,
        'bmi_25-30': 25,
        'bmi_30+': 30,
    }
    bmi_groups_range_1 = {
        'bmi_0-20': 0,
        'bmi_20-23': 20,
        'bmi_23-26': 23,
        'bmi_26-29': 26,
        'bmi_29-32': 29,
        'bmi_32+': 32,
    }

    for i, bmi_dict in enumerate([bmi_groups_range, bmi_groups_range_1]):
        bmi_groups = pd.Series(index=adata.obs_names, dtype='category')
        bmi_groups = pd.Categorical(bmi_groups, categories=list(bmi_dict.keys()))
        for group, min_bmi in bmi_dict.items():
            # adata.obs['bmi_groups'][adata.obs['bmi_lv'] >= min_bmi] = group
            bmi_groups[adata.obs['bmi_lv'] >= min_bmi] = group
        adata.obs[f'bmi_groups_{i}'] = bmi_groups

    # save a copy of the original bmi_groups column
    adata.obs['bmi_groups_copy'] = adata.obs['bmi_groups'].copy()
    adata.obs['bmi_groups'] = adata.obs['bmi_groups_0'].copy()
    return adata


def define_obesity_groups(adata):
    obesity_groups_dict = {
        'underweight': 0,
        'normal': 18.5,
        'overweight': 25,
        'obesity_1': 30,
        'obesity_2': 35,
        'obesity_3': 40,
    }
    obesity_groups = pd.Series(index=adata.obs_names, dtype='category')
    obesity_groups = pd.Categorical(obesity_groups, categories=list(obesity_groups_dict.keys()))
    for group, min_bmi in obesity_groups_dict.items():
        # adata.obs['obesity_groups'][adata.obs['bmi_lv'] >= min_bmi] = group
        obesity_groups[adata.obs['bmi_lv'] >= min_bmi] = group
    adata.obs['obesity_groups'] = obesity_groups
    return adata


def normalize_bmi(adata):
    avg = adata.obs['bmi_lv'].mean()
    std = adata.obs['bmi_lv'].std()
    bmi_normalized = (adata.obs['bmi_lv'] - avg) / std
    adata.obs['bmi_normalized'] = bmi_normalized
    return adata


if __name__ == '__main__':
    in_path = sys.argv[1] # path to AnnData file
    adata = sc.read_h5ad(in_path)
    print(adata)

    # adata = define_bmi_groups(adata)
    adata = define_obesity_groups(adata)
    print(adata.obs['bmi_groups'].value_counts())
    # print(adata.obs['bmi_groups_0'].value_counts())
    # print(adata.obs['bmi_groups_1'].value_counts())
    print(adata.obs['obesity_groups'].value_counts())
    print(adata.obs['obesity_groups'].head())

    adata = normalize_bmi(adata)
    print(adata.obs['bmi_normalized'].head())

    # save updated AnnData file
    adata.write_h5ad(in_path)
    print(f"adata saved.")
