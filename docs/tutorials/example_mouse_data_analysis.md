# Run Rawk on an example mouse dataset

This tutorial shows how to use Rawk to analyze the pseudobulk average
normalized UMI counts of the mouse immune checkpoint cancer therapy
scRNA-seq dataset.

## Import python dependencies

```python
import rawk as rk
import pathlib
import os
import pandas as pd
import numpy as np
```



## Set up input and output directories

```python
input_dir = "tutorial_data"

out_dir = os.path.join(
    "tutorial_output",
    "example_mouse_data_analysis")
pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
```

## Read input data

```python
original_data_fn = (
    "mouse_ict_scrna_seq"
    "_pseudobulk_average_normalized_uc.tsv"
)

original_data_df = pd.read_table(
    os.path.join(input_dir, original_data_fn))

mrn_node_df = pd.read_table(
    os.path.join(
        input_dir,
        "imm1415_unfiltered_nodes.tsv"))

mrn_edge_df = pd.read_table(
    os.path.join(
        input_dir,
        "imm1415_wg0_edges.tsv.gz"),
    low_memory=False)
```

## Set edge pruning weight threshold

```python
mrn_edge_weights = mrn_edge_df["mn_weight"].values

assert np.all(mrn_edge_weights > 0)

mrn_edge_weight_threshold = np.round(
    np.percentile(
        mrn_edge_weights,
        20),
    3)
```

## Preprocess input data

```python
# Select the first gene column and one or more other
# columns for running rawk
#
# Here, select columns gene and
# ActivatedMyeloidCells___aPD1aCTLA4
columns_for_rawk = [
    "gene",
    "ActivatedMyeloidCells___aPD1aCTLA4",
]

df_for_rawk_analysis = original_data_df.loc[:, columns_for_rawk]

def preprocess_1d_arr(x):
    p = np.e ** rk.qn_transform(
        x, log1p=True, collapse_0s=True,
        center=False)
    return p


pp_node_df, pp_edge_df = rk.get_met_net_dfs(
    mrn_node_df, mrn_edge_df,
    df_for_rawk_analysis,
    mrn_edge_weight_threshold,
    fill_missing_gene_prop=0,
    transform_gene_prop_func=preprocess_1d_arr)
```

## Run Rawk

```python
msrk = rk.MultiSampleRawk(
    pp_node_df.drop(columns=["rxn_name", "equation"]),
    pp_edge_df,
    n_jobs=1,
    n2v_walk_length=20,
    n2v_num_walks=8000,
    seed=42,
)

mrn_n_genes_by_pathway = (
    mrn_node_df
    .groupby("pathway")["gene"]
    .apply(lambda x: len(set(x)))
    .to_dict()
)

pathways_for_testing = [
    p for p, n in mrn_n_genes_by_pathway.items()
    if n >= 5 and p != "Miscellaneous"]

rk_fdr_df, rk_es_df = msrk.test_num_steps(
    pw_subset=pathways_for_testing)
```

## Plot Rawk results

```python
plot_pathway = "Fatty acid oxidation"
plot_rawk = msrk.rawk_list[0]
# each Rawk instance has a sample and a background sample
assert plot_rawk.sample.name == df_for_rawk_analysis.columns[1]
assert plot_rawk.bg_sample.name == "uniform_background"
```

```python
rk.plot_pw_neighborhood(
    plot_rawk,
    plot_pathway,
    os.path.join(
        out_dir,
        "rawk_pathway_nbr_rwv_cpm_p1_log2fc.png"),
    "rwv_cpm_p1_log2fc",
    node_color_center=0,
    n_cutoff=6, node_alpha=0.7,
    node_color_title=(
        "Sample vs background\n"
        "RWVCPM+1 log2 fold change"
    ),
    title=f"Source pathway:\n{plot_pathway}")
```

![](tutorial_output/example_mouse_data_analysis/rawk_pathway_nbr_rwv_cpm_p1_log2fc.png)

```python
rk.plot_graph(
    plot_rawk.sample.rxn_graph,
    plot_rawk.sample.rxn_pos,
    os.path.join(
        out_dir,
        "rawk_mrn_reaction_tsne.png"),
    draw_edges=True,
    non_s_pw_node_size=4,
    title="Rawk metabolic reactions",
    dim_name="t-SNE")
```

![](tutorial_output/example_mouse_data_analysis/rawk_mrn_reaction_tsne.png)

```python
pathway_pos = {
    n: d["pw_pos"]
    for n, d in plot_rawk.sample.pw_graph.nodes.items()
}
rk.plot_graph(
    plot_rawk.sample.pw_graph,
    pathway_pos,
    os.path.join(
        out_dir,
        "rawk_mrn_pathway_tsne.png"),
    node_color_attr="pw_property",
    non_s_pw_node_size=30,
    draw_edges=True,
    title="Rawk metabolic pathways",
    dim_name="t-SNE")
```

![](tutorial_output/example_mouse_data_analysis/rawk_mrn_pathway_tsne.png)

```python
rk.plot_rawk_sample_mtx(
    plot_rawk.sample,
    os.path.join(
        out_dir, "rawk_test_sample"))
```

![](tutorial_output/example_mouse_data_analysis/rawk_test_sample_pp_rwvcpm_mtx.png)

## Save Rawk results

```python
rk_fdr_df.to_csv(
    os.path.join(out_dir, "rawk_fdr.tsv"),
    sep="\t", index=False)

rk_es_df.to_csv(
    os.path.join(out_dir, "rawk_es.tsv"),
    sep="\t", index=False)
```
