import rawk as rk
import pathlib
import pandas as pd
import numpy as np
import time



input_dir = 'tutorial_data/'


out_dir = 'tutorial_output/example_human_data_analysis/'
pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

# Analyze human L1000 assay dataset --------------------------------------------

# Read and transform gene expression zscore table
h_gene_prop_df = (
    pd.read_csv(input_dir + 'withaferin_a_HA1E_level5_zscore_all_genes.csv')
    .set_index('gene')
    .apply(lambda x: x + 10, axis=0)
    .copy()
)

# Read metabolic reaction network edge table
h_rxn_edge_df = pd.read_table(input_dir + 'recon3d_edges.tsv.gz')

# Read metabolic gene reaction association table
h_rxn_gene_df = pd.read_table(input_dir + 'recon3d_rxn_gene.tsv')


# Pre-process input data
node_df, edge_df = rk.get_met_net_dfs(
    h_rxn_gene_df, h_rxn_edge_df, h_gene_prop_df,
    205, prop_na_handling='remove')



# Save pre-rocesses node and edge tables
node_df.to_csv(out_dir + 'h_rxn_props.csv', index=False)
edge_df.to_csv(out_dir + 'h_rxn_edges.csv', index=False)


##### Run Rawk on one sample #####

# Select a sample to run Rawk
selected_sample = 'LJP005_HA1E_24H_G13___HA1E___withaferin-a___10_0_um___24_h'
ss_node_df = (
    node_df
    .loc[:, ['rxn', 'pathway', 'gene', selected_sample]]
    .copy()
    .rename(columns={selected_sample: 'property'})
)


h_rawk = rk.Rawk(
    rk.RawkSample(
        ss_node_df, edge_df,
        name='LJP005_HA1E_24H_G13___HA1E___withaferin-a___10_0_um___24_h',
        workers=1),
    rk.RawkSample(
        ss_node_df, edge_df, name='uniform_background',
        uniform_property=True, workers=1)
)


# Plot metabolic network nodes and edges
h_rawk.sample.plot_graph('{}/h_nodes.png'.format(out_dir), draw_edges=False)

# Plot 'Fatty acid synthesis' pathway neighborhood
for i in ['property', 'rw_ps', 'log10p1_rw_ps', 'rw_ps_diff', 'rw_ps_log2fc']:
    if i in ['rw_ps_log2fc', 'rw_ps_diff']:
        j = 0
    else:
        j = None
    h_rawk.plot_pw_neighborhood(
        'Fatty acid synthesis',
        '{}/h_{}_neighborhood.png'.format(out_dir, i),
        i, node_color_vmin=None, node_color_vmax=None,
        node_color_center=j, n_cutoff=10, node_alpha=0.5)


# Plot node2vec random walk related matrices
h_rawk.sample.plot_n2v_mat_dict(out_dir + 'h_sample')

# Run enrichment tests
h_rawk_test_res = h_rawk.test_res_dict_to_df(
    h_rawk.test_num_steps(cmp_norm_fac=10000), '_rawk')

# Save enrichment test results
h_rawk_test_res[0].to_csv('{}/h_rawk_ss_padj_df.csv'.format(out_dir))
h_rawk_test_res[1].to_csv('{}/h_rawk_ss_es_df.csv'.format(out_dir))



##### Run Rawk on mutliple samples #####
ms_node_df = (
    node_df
    .loc[:, ['rxn', 'pathway', 'gene',
             'LJP005_HA1E_24H_G13___HA1E___withaferin-a___10_0_um___24_h',
             'LJP005_HA1E_24H_G14___HA1E___withaferin-a___3_33_um___24_h',
             'LJP005_HA1E_24H_G15___HA1E___withaferin-a___1_11_um___24_h',
             'LJP005_HA1E_24H_G16___HA1E___withaferin-a___0_37_um___24_h',
             'LJP005_HA1E_24H_G17___HA1E___withaferin-a___0_12_um___24_h',
             'LJP005_HA1E_24H_G18___HA1E___withaferin-a___0_04_um___24_h',]]
    .copy()
)

h_ms_rawk = rk.MultiSampleRawk(ms_node_df, edge_df, workers=1)

# Run enrichment tests
h_ms_rawk_ns_test_res = h_ms_rawk.test_num_steps(cmp_norm_fac=10000)

# Save enrichment test results
h_ms_rawk_ns_test_res[0].to_csv('{}/h_rawk_ms_padj_df.csv'.format(out_dir))
h_ms_rawk_ns_test_res[1].to_csv('{}/h_rawk_ms_es_df.csv'.format(out_dir))
