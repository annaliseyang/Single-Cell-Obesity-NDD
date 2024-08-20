import sys
import scanpy as sc
import pandas as pd

in_path = sys.argv[1] # path to AnnData file
adata = sc.read_h5ad(in_path)

# bmi_groups = ['bmi_0-20', 'bmi_20-25', 'bmi_25-30', 'bmi_30+']
bmi_groups_range = {
    'bmi_0-20': 0,
    'bmi_20-25': 20,
    'bmi_25-30': 25,
    'bmi_30+': 30,
}
# bmi_groups_1 = ['bmi_0-20', 'bmi_20-23', 'bmi_23-26', 'bmi_26-29', 'bmi_29-32', 'bmi_32+']
bmi_groups_range_1 = {
    'bmi_0-20': 0,
    'bmi_20-23': 20,
    'bmi_23-26': 23,
    'bmi_26-29': 26,
    'bmi_29-32': 29,
    'bmi_32+': 32,
}

# bmi_lv = adata.obs['bmi_lv']
# adata.obs['bmi_groups'] = []
for i, bmi_dict in enumerate([bmi_groups_range, bmi_groups_range_1]):
    bmi_groups = pd.Series(index=adata.obs_names, dtype='category')
    bmi_groups = pd.Categorical(bmi_groups, categories=list(bmi_dict.keys()))
    for group, min_bmi in bmi_dict.items():
        # adata.obs['bmi_groups'][adata.obs['bmi_lv'] >= min_bmi] = group
        bmi_groups[adata.obs['bmi_lv'] >= min_bmi] = group
    adata.obs[f'bmi_groups_{i}'] = bmi_groups

# # adata.obs['bmi_groups_1'] = []
# for group, min_bmi in bmi_groups_range_1.items():
#     adata.obs['bmi_groups_1'][adata.obs['bmi_lv'] >= min_bmi] = group

print(adata)
print(adata.obs['bmi_groups_0'].value_counts())
# print(adata.obs[adata.obs['bmi_groups_0'] == 'bmi_30+']['bmi_lv'])
print(adata.obs['bmi_groups_1'].value_counts())
# print(adata.obs[adata.obs['bmi_groups_1'] == 'bmi_32+']['bmi_lv'])
# print(adata.obs[adata.obs['bmi_groups_1'] == 'bmi_20-23']['bmi_lv'])

# save updated AnnData file
adata.write_h5ad(in_path)
