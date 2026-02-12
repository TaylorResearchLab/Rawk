import pandas as pd
import numpy as np



def get_met_net_dfs(rxn_gene_df, rxn_edge_df, gene_prop_df, mn_weight_cutoff,
                    prop_na_handling="replace_with_0",
                    rxn_gene_prop_agg_func=lambda x: x.mean()):
    """
    Prepare input node and edge dataframes for Rawk

    Parameters
    ----------
    rxn_gene_df : dataframe
        A dataframe of reactions and their associated genes. Following are the
        required columns: 'rxn' (reaction ID), 'rxn_name' (reaction name),
        'equation' (reaction equation), 'pathway' (reaction pathway), and
        'gene' (gene symbol). If a reaction is associated with multiple genes,
        one row lists one associated gene.
    rxn_edge_df : dataframe
        A dataframe of reactions edges. Following are the required columns:
        'src' (source node reaction ID), 'dest' (destination node reaction ID),
        'mn_weight' (metabolic network edge weight). The edges are undirected,
        with the src <= dest in alphabetical order.
    gene_prop_df : dataframe
        A dataframe of gene properties, such as log fold changes, normalized
        read counts, and z-scores. The index of the dataframe is a list of gene
        symbols. Each column of the dataframe is a sample.
    mn_weight_cutoff : float
        The cutoff of metabolic network weights. Keep only edges with weights >
        mn_weight_cutoff.
    prop_na_handling : string
        Supported NA handling methods: 1) 'remove': remove missing
        (NA) property values; 2) 'replace_with_0': replace missing property
        values with 0.
    rxn_gene_prop_agg_func : function
        The function used to aggregate the property values of multiple genes
        that are associated with each reaction.

    Returns
    -------
    (dataframe, dataframe)
        A 2-tuple of a reaction node property dataframe and a reaction edge
        dataframe
    """
    if gene_prop_df.isnull().values.sum() != 0:
        raise ValueError("gene_prop_df contains one or more NA/NaN... values")
    if rxn_edge_df.isnull().values.sum() != 0:
        raise ValueError("rxn_edge_df contains one or more NA/NaN... values")
    if rxn_gene_df.isnull().values.sum() != 0:
        raise ValueError("rxn_gene_df contains one or more NA/NaN... values")
    #
    if (rxn_edge_df[["src", "dest"]].drop_duplicates().shape[0] !=
            rxn_edge_df.shape[0]):
        raise ValueError("rxn_edge_df contains duplicated (src, dest) pairs")
    #
    if not all(rxn_edge_df.src <= rxn_edge_df.dest):
        raise ValueError(
            "rxn_edge_df requres src <= dest in alphabetical order")
    #
    if (rxn_gene_df[["rxn", "gene"]].drop_duplicates().shape[0] !=
            rxn_gene_df.shape[0]):
        raise ValueError("rxn_gene_df contains duplicated (rxn, gene) pairs")
    #
    if len(gene_prop_df.index.unique()) != gene_prop_df.shape[0]:
        raise ValueError("gene_prop_df contains duplicated rownames.")
    #
    rxn_edge_df = rxn_edge_df.loc[
        rxn_edge_df["mn_weight"] > mn_weight_cutoff, :].copy()
    #
    edge_rxn_set = set(
        rxn_edge_df["src"].tolist() + rxn_edge_df["dest"].tolist())
    rxn_gene_prop_df = (
        rxn_gene_df
        .loc[rxn_gene_df["rxn"].isin(edge_rxn_set), :]
        .merge(gene_prop_df, how="left", on="gene")
    )
    #
    if prop_na_handling == "remove":
        rxn_gene_prop_df = rxn_gene_prop_df.dropna()
    elif prop_na_handling == "replace_with_0":
        rxn_gene_prop_df = rxn_gene_prop_df.fillna(0)
    else:
        raise ValueError("Unknown rxn_gene_prop_df {}".format(rxn_gene_prop_df))
    #
    rxn_edge_df = rxn_edge_df.loc[
        np.logical_and(
            rxn_edge_df["src"].isin(set(rxn_gene_prop_df["rxn"].tolist())),
            rxn_edge_df["dest"].isin(set(rxn_gene_prop_df["rxn"].tolist()))), :]
    #
    edge_rxn_list = list(
        set(rxn_edge_df["src"].tolist() + rxn_edge_df["dest"].tolist()))
    self_loop_only_nodes = []
    for i in edge_rxn_list:
        i_df = rxn_edge_df.loc[
            np.logical_or(rxn_edge_df["src"] == i, rxn_edge_df["dest"] == i), :]
        if i_df.shape[0] == 1:
            i_src = rxn_edge_df["src"].tolist()[0]
            i_dest = rxn_edge_df["dest"].tolist()[0]
            if i_src == i_dest:
                self_loop_only_nodes.append(i)
    #
    self_loop_only_nodes = set(self_loop_only_nodes)
    #
    rxn_edge_df = rxn_edge_df.loc[
        np.logical_not(np.logical_or(
            rxn_edge_df["src"].isin(self_loop_only_nodes),
            rxn_edge_df["dest"].isin(self_loop_only_nodes))), :].copy()
    #
    def unique_one(x):
        x_set = set(x.tolist())
        assert len(x_set) == 1, str(x_set)
        return list(x_set)[0]
    #
    agg_func_dict = {
        "rxn_name": unique_one, "equation": unique_one, "pathway": unique_one,
        "gene": lambda x: ";;;".join(x.tolist())
    }
    for i in rxn_gene_prop_df.columns:
        if i != "rxn" and i not in agg_func_dict:
            agg_func_dict[i] = rxn_gene_prop_agg_func
    #
    rxn_prop_df = (
        rxn_gene_prop_df
        .groupby("rxn")
        .agg(agg_func_dict)
        .reset_index(names=["rxn"])
    )
    rxn_prop_df.shape
    assert rxn_prop_df.isnull().values.sum() == 0
    #
    common_rxns = (
        set(rxn_prop_df["rxn"].tolist())
        .intersection(set(rxn_edge_df["src"].tolist()))
        .intersection(set(rxn_edge_df["dest"].tolist()))
    )
    #
    rxn_edge_df = rxn_edge_df.loc[
        np.logical_and(
            rxn_edge_df["src"].isin(common_rxns),
            rxn_edge_df["dest"].isin(common_rxns)), :].copy()
    #
    rxn_prop_df = rxn_prop_df.loc[
        rxn_prop_df["rxn"].isin(common_rxns), :].copy()
    return rxn_prop_df, rxn_edge_df
