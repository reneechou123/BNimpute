#!/usr/bin/env python3

import pandas as pd
import os
import numpy as np

metadata = pd.read_csv('rnaseq_metadata_copy.tsv', sep='\t', header=0)
data = pd.read_csv('rnaseq_data_tpm_from_metadata.tsv', sep='\t', header=0, index_col=0)

'''
split metadata based on tissue and cell lines
'''
sample_list = {}
dict_of_tissue = {k: v for k, v in metadata.astype(str).groupby('tissue')}
temp_dataset = pd.DataFrame(columns=metadata.columns)
for k in dict_of_tissue.keys():
    # split embryo by developmental stage
    if k == 'cell line-embryo' or k == 'embryo':
        dict_of_tissue_and_dev_stage = {k: v for k, v in dict_of_tissue[k].astype(str).groupby('developmental stage')}
        for j in dict_of_tissue_and_dev_stage.keys():
            # discard if sample number < 10
            sample_size = dict_of_tissue_and_dev_stage[j].shape[0]
            if sample_size < 10:
                continue
            sample_list[k + '/' + j] = [sample_size, 0]
        continue
    # collect samples with tissue not specified or sample number < 10
    sample_size = dict_of_tissue[k].shape[0]
    if k == 'nan' or dict_of_tissue[k].shape[0] < 10:
        temp_dataset = temp_dataset.append(dict_of_tissue[k])
        continue
    sample_list[k] = [sample_size, 1]

# split the rest of dataset by cell type
dict_of_cell_type = {k: v for k, v in temp_dataset.astype(str).groupby('cell type')}
for i in dict_of_cell_type.keys():
    # discard if sample number < 10
    sample_size = dict_of_cell_type[i].shape[0]
    if i == 'nan' or sample_size < 10:
        continue
    sample_list[i] = [sample_size, 2]

'''
print out sample list
'''
# convert sample_list to a sorted list
sample_list = sorted(sample_list.items(), key=lambda kv: kv[1])
total = 0
with open('rnaseq_sample_stat.txt', 'w') as f:
    for item in sample_list:
         _ = f.write(item[0] + ', ' + str(item[1][0]) + '\n')
         total += item[1][0]
    _ = f.write('TOTAL: ' + str(total) + '\n')
    f.close

'''
split expression data into different tissue types
'''
# function for generating expression datasets for GENIE3
def generate_files_for_GENIE3(metadata, key, ind, dir_name, exp, prior):
    os.system('mkdir -p ' + dir_name)
    if ind == 0:
        keys = key.split('/')
        temp_list = metadata.loc[(metadata['tissue'] == keys[0]) & (metadata['developmental stage'] == keys[1]), 
                'sample_name'].tolist()
    elif ind == 1:
        temp_list = metadata.loc[metadata['tissue'] == key, 'sample_name'].tolist()
    else:
        temp_list = metadata.loc[metadata['cell type'] == key, 'sample_name'].tolist()    
    sub_exp = exp[temp_list]
    pri = prior
    # intersect for expression data and remove all zero counts
    index = list(pri.index)
    sub_exp = sub_exp.loc[sub_exp.index.isin(index), :]
    sub_exp = sub_exp[(sub_exp.T != 0).any()]
    # get final intersected genes for prior
    index2 = list(sub_exp.index)
    pri = pri.loc[index2, pri.columns.isin(index2)]
    # remove TFs having no targets
    pri = pri.loc[:, (pri != 0).any()]
    # write output
    sub_exp.to_csv(dir_name + '/' + key.replace('/', '.') + '_exp.tsv', header=False, index=True, sep='\t')
    pri.to_csv(dir_name + '/' + key.replace('/', '.') + '_prior.tsv', header=True, index=True, sep='\t')
    # return log
    return key.replace('/', '.') + ', ' + str(sub_exp.shape[1])

# select groups of samples
np.random.seed(123)
groups = [sample_list[i] for i in np.sort(np.random.randint(0, len(sample_list) - 1, 10))]
# import expression and prior
exp = pd.read_csv('rnaseq_data_tpm_from_metadata.tsv', sep='\t', header=0, index_col=0)
prior = pd.read_csv('flynet_supervised_0.6.txt', sep='\t', header=None)
# making prior file
prior = prior.pivot(columns=0, index=1, values=2).rename_axis(None).fillna(0)

with open('input_for_NetREX_log.txt', 'w') as f:
    for i in range(0, len(groups)):
        log = generate_files_for_NetREX(metadata, groups[i][0], groups[i][1][1], str(i + 1), exp, prior)
        _ = f.write(log + '\n')
    f.close

