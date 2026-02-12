import rawk as rk
import pathlib
import pandas as pd
import numpy as np
import time



input_dir = 'tutorial_data/'


out_dir = 'tutorial_output/example_mouse_data_analysis/'
pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)



# Analyze mouse immunotherapy scRNA-seq dataset --------------------------------

# Read cluster average normalized expression table
m_gene_prop_df = pd.read_table(
    input_dir + 'immune_cell_cluster_average_normalized_expression.tsv')

# Read metabolic reaction network edge table
m_rxn_edge_df = pd.read_table(input_dir + 'imm1415_edges.tsv')

# Read metabolic gene reaction association table
m_rxn_gene_df = pd.read_table(input_dir + 'imm1415_rxn_gene.tsv')

# Pre-process input data
node_df, edge_df = rk.get_met_net_dfs(
    m_rxn_gene_df, m_rxn_edge_df, m_gene_prop_df,
    15, prop_na_handling='replace_with_0',
    rxn_gene_prop_agg_func=lambda x: np.log1p(np.mean(np.expm1(x.values))))


# Save pre-rocesses node and edge tables
node_df.to_csv(out_dir + 'm_rxn_props.csv', index=False)
edge_df.to_csv(out_dir + 'm_rxn_edges.csv', index=False)


##### Run Rawk on one sample #####

# Select a sample to run Rawk
selected_sample = 'ActivatedMyeloidCells___aPD1aCTLA4'
ss_node_df = (
    node_df
    .loc[:, ['rxn', 'pathway', 'gene', selected_sample]]
    .copy()
    .rename(columns={selected_sample: 'property'})
)


m_rawk = rk.Rawk(
    rk.RawkSample(
        ss_node_df, edge_df, name='ActivatedMyeloidCells___aPD1aCTLA4',
        workers=1),
    rk.RawkSample(
        ss_node_df, edge_df, name='uniform_background',
        uniform_property=True,
        workers=1)
)

# Plot metabolic network nodes and edges
m_rawk.sample.plot_graph('{}/m_nodes.png'.format(out_dir), draw_edges=False)


# Plot 'Glycolysis/Gluconeogenesis' pathway neighborhood
for i in ['property', 'rw_ps', 'log10p1_rw_ps', 'rw_ps_diff', 'rw_ps_log2fc']:
    if i in ['rw_ps_log2fc', 'rw_ps_diff']:
        j = 0
    else:
        j = None
    m_rawk.plot_pw_neighborhood(
        'Glycolysis/Gluconeogenesis',
        '{}/m_{}_neighborhood.png'.format(out_dir, i),
        i, node_color_vmin=None, node_color_vmax=None,
        node_color_center=j, n_cutoff=10, node_alpha=0.5)


# Plot node2vec random walk related matrices
m_rawk.sample.plot_n2v_mat_dict(out_dir + 'm_sample')


# Run enrichment tests
m_rawk_test_res = m_rawk.test_res_dict_to_df(
    m_rawk.test_num_steps(cmp_norm_fac=100), '_rawk')

# Save enrichment test results
m_rawk_test_res[0].to_csv('{}/m_rawk_ss_padj_df.csv'.format(out_dir))
m_rawk_test_res[1].to_csv('{}/m_rawk_ss_es_df.csv'.format(out_dir))



##### Run Rawk on mutliple samples #####

ms_node_df = (
    node_df
    .loc[:, ['rxn', 'pathway', 'gene',
             'ActivatedMyeloidCells___aPD1aCTLA4',
             'HighlyProliferativeTCells___aPD1aCTLA4',
             'OtherMyeloidCells___aPD1aCTLA4',]]
    .copy()
)

m_ms_rawk = rk.MultiSampleRawk(ms_node_df, edge_df, workers=1)

# Run enrichment tests
m_ms_rawk_ns_test_res = m_ms_rawk.test_num_steps(cmp_norm_fac=100)

# Save enrichment test results
m_ms_rawk_ns_test_res[0].to_csv('{}/m_rawk_ms_padj_df.csv'.format(out_dir))
m_ms_rawk_ns_test_res[1].to_csv('{}/m_rawk_ms_es_df.csv'.format(out_dir))
