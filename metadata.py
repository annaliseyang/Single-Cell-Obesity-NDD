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

        print(len(columns))
        print(columns)
        print(lines[1].strip().split('\t'))

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

        # print(AD_types)
        # print(f"Number of AD cases: {AD_count}")

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


if __name__ == "__main__":
    AD427_metadata = get_AD427_metadata(AD427_META_PATH)
    AD427_keys = list(AD427_metadata.values())[0].keys()
    print(AD427_keys)
    assert len(list(AD427_metadata.keys())) == 427
    assert len(set(list(AD427_metadata.keys()))) == 427, "Duplicate sequences found in metadata"

    ADMR_metadata = get_ADMR_metadata(ADMR_META_PATH)
    ADMR_keys = list(ADMR_metadata.values())[0].keys()
    # print(ADMR_metadata.keys())
    print(list(ADMR_metadata.values())[0])
    # print(len(ADMR_metadata))

    common_keys = set(AD427_keys).intersection(set(ADMR_keys))
    print('Number of AD427 keys:', len(AD427_keys))
    print('Number of ADMR keys:', len(ADMR_keys))
    print('Common keys:', common_keys, 'Number:', len(common_keys))
