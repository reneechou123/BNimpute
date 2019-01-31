## Modify the metadata

import pandas as pd

metadata = pd.read_table('rnaseq_metadata.tsv', header=0)

# Modify tissue and developmetal stage according to cell lines
cell_line_samples = metadata[metadata['cell type'].notna()]

def modify_metadata(cell_types, tissue_name, dev_stage_name, metadata):
    sample_list = cell_line_samples[cell_line_samples['cell type'].str.contains('|'.join(cell_types))]['sample_name'].tolist()
    metadata.loc[metadata['sample_name'].isin(sample_list), 'tissue'] = tissue_name
    metadata.loc[metadata['sample_name'].isin(sample_list), 'developmental stage'] = dev_stage_name

# BG series
BG = ['BG']
modify_metadata(BG, 'cell line-CNS', 'third instar larval stage', metadata)

# D10 and D20
D10_D20 = ['D10', 'D20']
modify_metadata(D10_D20, 'cell line-antennal disc', 'third instar larval stage', metadata)

# Kc series
Kc = ['Kc']
modify_metadata(Kc, 'cell line-embryo', 'dorsal closure stage', metadata)

# embryo
em = ['Jupiter', 'G2', 'G1', 'Ras\[V12\]', 'Sg4', 'E-OR', 'E-CS', 'Pten-X', 'CCa', 'DX', 'D1', 'GM2', 'GM3', 'PR8']
modify_metadata(em, 'cell line-embryo', 'embryonic stage', metadata)

# late embryonic
ltem = ['S2', 's2', 'S3', 'S1']
modify_metadata(ltem, 'cell line-embryo', 'late embryonic stage', metadata)

# ovary
ov = ['OSS', 'OSC']
modify_metadata(ov, 'cell line-ovary', 'adult stage', metadata)

# wing disc
wd = ['D16', 'D1', 'Cl.8', 'D21', 'D23', 'D32', 'D8', 'D9', 'MCW12']
modify_metadata(wd, 'cell line-wing disc', 'third instar larval stage', metadata)

# cell immortalization
ci = ['R1', 'R4', 'R5']
modify_metadata(ci, 'cell line-embryo', 'cell immortalization', metadata)

# These may not be cell lines
# neuron
neuron = ['neuron']
modify_metadata(neuron, 'neuron', None, metadata)

# embryonic
embryonic = ['embryonic']
modify_metadata(embryonic, 'embryo', None, metadata)

# output new metadata
metadata.to_csv('rnaseq_metadata_copy.tsv', sep='\t', header=True, index=False)

