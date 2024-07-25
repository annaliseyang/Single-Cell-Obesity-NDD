from init import *
import csv

AD427_META_PATH = '/net/bmc-lab5/data/kellis/group/Zunpeng_Sharing/MR_ATAC_RNA_2024/metadata/PFC427_sample.meta.txt'
ADMR_META_PATH = '/net/bmc-lab5/data/kellis/group/Zunpeng_Sharing/MR_ATAC_RNA_2024/metadata/All_snATAC_snRNA_Multiome.nSamples_799.with_pathology.Mar10_2024.tsv'

def get_AD427_metadata(AD427_META_PATH):
    with open(AD427_META_PATH, 'r') as file:
        lines = file.readlines()
        columns = lines[0].strip().split('\t')
        column_names = ['apoe_genotype', 'cogdx', 'age_death', 'educ', 'msex', 'race', 'braaksc', 'gpath', 'pmi', 'amyloid', 'plaq_d', 'plaq_n', 'nft', 'tangles', 'arteriol_scler', 'ADdiag3types', 'ADdiag2types']
        column_index = [columns.index(name) for name in column_names]
        # print(lines)
        AD_count = 0
        AD_types = set()
        out = {} # dict mapping projid to info
        for line in lines[1:]:
            # loop through each row/individual
            seq = line.strip().split('\t')[0]
            row = line.strip().split('\t')[1:]
            assert len(columns) == len(row), f"Unequal number of columns: {len(columns)} vs {len(row)}"
            ADdiag3types, ADdiag2types = line.strip().split('\t')[-2:]
            AD_types.update(ADdiag3types, ADdiag2types)
            info = {
                k: v for k, v in zip(columns, row)
            }
            diag = {}
            for name, index in zip(column_names, column_index):
                val = row[index]
                # print(f"{name}: {val}")
                diag[name] = val
            if 'AD' in ADdiag3types or 'AD' in ADdiag2types:
                AD_count += 1
            # print(projid, ADdiag3types, ADdiag2types)
            projid = info['projid']
            # print(seq, diag)
            out[projid] = info
    return out

def get_ADMR_metadata(ADMR_META_PATH):
    with open(ADMR_META_PATH, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        data = list(reader)
    out = {} # dict mapping projid to info
    for dict in data:
        projid = dict['projid']
        out[projid] = dict
    return out

@log
def integrate_metadata(adata, AD427_metadata, ADMR_metadata):
    AD427_keys = list(AD427_metadata.values())[0].keys()
    ADMR_keys = list(ADMR_metadata.values())[0].keys()
    all_keys = set(AD427_keys) | set(ADMR_keys)
    features = {k: [] for k in all_keys}
    # diag = {'ADdiag3types': [], 'ADdiag2types': [], 'Pathology': []}
    for projid in adata.obs['projid']:
        projid = str(int(projid))
        info = AD427_metadata.get(projid) if projid in AD427_metadata else ADMR_metadata.get(projid)
        for key, val in features.items():
            val.append(info.get(key))
    for key, val in features.items():
        adata.obs[key] = val
    return adata


def get_unique_projids_bmi_pathology(adata):
    '''return a dataframe with bmi and pathology for unique projids in adata.obs'''
    unique_projids = adata.obs['projid'].unique()
    unique_bmi_pathology = adata.obs[['projid', 'bmi_lv', 'ADdiag3types', 'ADdiag2types']].loc[adata.obs['projid'].isin(unique_projids)].drop_duplicates()
    return unique_bmi_pathology.dropna()


def plot_bmi_AD(adata):
    df = get_unique_projids_bmi_pathology(adata)
    # boxplot for unique projids
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    sns.boxplot(data=df, x='ADdiag3types', y='bmi_lv', ax=axs[0])
    sns.boxplot(data=df, x='ADdiag2types', y='bmi_lv', ax=axs[1])
    plt.savefig('figures/boxplot_bmi_AD_unique.png')


def plot_bmi_pathology(adata, column):
    '''
    plot bmi vs pathology for unique projids, return dataframe with bmi and pathology
    column is a categorical column like 'ADdiag3types' or 'ADdiag2types
    '''
    unique_projids = adata.obs['projid'].unique()
    unique_bmi_pathology = adata.obs[['projid', 'bmi_lv', column]].loc[adata.obs['projid'].isin(unique_projids)].drop_duplicates()
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.boxplot(data=unique_bmi_pathology, x=column, y='bmi_lv', ax=ax)
    plt.savefig(f'figures/boxplot_bmi_{column}_unique.png')
    return unique_bmi_pathology.dropna()

AD427_metadata = get_AD427_metadata(AD427_META_PATH)
ADMR_metadata = get_ADMR_metadata(ADMR_META_PATH)

AD427_keys = list(AD427_metadata.values())[0].keys()
ADMR_keys = list(ADMR_metadata.values())[0].keys()

if __name__ == "__main__":
    # common_keys = set(AD427_keys).intersection(set(ADMR_keys))
    # print('Number of AD427 keys:', len(AD427_keys))
    # print('Number of ADMR keys:', len(ADMR_keys))
    # print('Common keys:', common_keys, 'Number:', len(common_keys))

    # integrate_metadata(adata, AD427_metadata, ADMR_metadata)
    # print("Metadata integrated!", flush=True)
    # print(adata, flush=True)

    # merge ADdiag3types and Pathology into a single column AD_states
    AD_states = []
    for i in range(len(adata.obs['ADdiag3types'])):
        if adata.obs['ADdiag3types'][i] is not None and pd.notna(adata.obs['ADdiag3types'][i]):
            AD_states.append(adata.obs['ADdiag3types'][i])
        else:
            AD_states.append(adata.obs['Pathology'][i])
    adata.obs['AD_states'] = AD_states

    adata.write('/home/anna_y/data/write/AD427_ADMR_meta_Jul22_2024.h5ad')

    pass

# # dataset = ['AD427' if str(int(projid)) in AD427_metadata else 'ADMR' for projid in adata.obs['projid']]
# # adata.obs['dataset'] = ['AD427' if str(int(projid)) in AD427_metadata else 'ADMR' for projid in adata.obs['projid']]
# AD427 = adata[adata.obs['dataset'] == 'AD427']
# ADMR = adata[adata.obs['dataset'] == 'ADMR']

# pathology_cols = ['apoe_genotype', 'cogdx', 'age_death', 'educ', 'msex', 'race', 'braaksc', 'gpath', 'pmi', 'amyloid', 'plaq_d', 'plaq_n', 'nft', 'tangles', 'arteriol_scler', 'ADdiag3types', 'ADdiag2types']
# for col in pathology_cols:
#     plot_bmi_pathology(adata, col)

# # all data
# fig, axs = plt.subplots(1, 2, figsize=(12, 6))
# sns.boxplot(data=AD427.obs, x='ADdiag3types', y='bmi_lv', ax=axs[0])
# sns.boxplot(data=AD427.obs, x='ADdiag2types', y='bmi_lv', ax=axs[1])
# plt.savefig('figures/boxplot_bmi_AD.png')

# plt.figure(figsize=(8, 6))
# # sns.boxplot(data=AD427.obs, x='ADdiag3types', y='bmi_lv')
# sns.boxplot(data=AD427.obs, x='ADdiag2types', y='bmi_lv')
# plt.savefig('figures/boxplot_bmi_ADdiag2types.png')

# bmi_lv = [info['bmi_lv'] for info in AD427_metadata.values()]
# ADdiag3types = [info['ADdiag3types'] for info in AD427_metadata.values()]
# df = pd.DataFrame({'bmi_lv': bmi_lv, 'ADdiag3types': ADdiag3types})
# plt.figure(figsize=(8, 6))
# sns.boxplot(data=df, x='ADdiag3types', y='bmi_lv')
# plt.savefig('figures/boxplot_bmi_lv_ADdiag3types.png')

# plt.figure(figsize=(8, 6))
# sns.boxplot(data=ADMR.obs, x='Pathology', y='bmi_lv')
# plt.show()
# plt.savefig('figures/ADMR_bmi_boxplot.png')
