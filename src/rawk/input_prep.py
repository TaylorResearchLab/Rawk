import pandas as pd
import numpy as np
import networkx as nx
from sklearn.preprocessing import QuantileTransformer



def transform_gene_prop(gene_prop_df, transform):
    """
    Apply transform function to each property column

    Parameters
    ----------
    gene_prop_df : dataframe
        A dataframe of gene properties, such as log fold changes, normalized
        read counts, and z-scores. The 'gene' column of the dataframe is a
        list of gene symbols. Other columns of the dataframe are the gene
        properties of samples.
    transform : function
        A function to transform each property column.

    Returns
    -------
    dataframe after transformation
    """
    assert not gene_prop_df['gene'].duplicated().any()

    if transform is not None:
        tdf = (
            gene_prop_df
            .copy()
            .set_index('gene')
            .apply(transform, axis=0)
            .reset_index(names='gene')
        )

    return tdf



def qn_transform(s, sigma=0.367879, log1p=False, collapse_0s=False,
                 center=False, seed=42):
    """
    Quantile transform gene properties to normal distribution

    Parameters
    ----------
    s : series
        A pandas.Series of gene properties, such as log fold changes,
        normalized read counts, and z-scores.
    sigma : float
        Output normal distribution sigma. Default to
        np.round(1 / np.e, 6), so np.e ** (sigma * (3 - (-3))) ~= 9.
    log1p : bool
        Apply log1p transform on the properties or not, before
        quantile normalization.
    collapse_0s : bool
        If True, all 0s will be collapse into one 0, and the 0's
        quantile normalized value will be assigned to all 0s.
    center : bool
        If True, the tranformed values will be centered at the input 0s.
    seed : int
        Random state.

    Returns
    -------
    Series after transformation
    """
    x = s.values.copy()
    if log1p:
        x = np.log1p(x)

    assert len(x.shape) == 1
    n = len(x)

    qt = QuantileTransformer(
        output_distribution='normal',
        random_state=seed)

    if collapse_0s:
        # collapse all 0s when transform
        non0_idc = x != 0
        non0_x = x[non0_idc].copy()

        c0_x = np.concatenate((np.array([0]), non0_x))

        c0_x = c0_x.reshape(-1, 1)
        assert c0_x[0] == 0
        assert c0_x[1] != 0
        assert c0_x.shape == (sum(non0_idc) + 1, 1)

        c0_x = qt.fit_transform(c0_x)
        assert c0_x.shape == (sum(non0_idc) + 1, 1)

        c0_x = c0_x.flatten()
        assert c0_x.shape == (sum(non0_idc) + 1,)

        x[np.logical_not(non0_idc)] = c0_x[0]

        x[non0_idc] = c0_x[1:]

    else:
        x = x.reshape(-1, 1)
        assert x.shape == (n, 1)

        x = qt.fit_transform(x)
        assert x.shape == (n, 1)

        x = x.flatten()

    assert len(x.shape) == 1

    if center:
        t0 = qt.transform([[0.0]])[0][0]
        x = x - t0

    x = x * sigma

    return pd.Series(x, index=s.index.copy())



def get_mrn_gp_df(gene_prop_df, rxn_gene_df,
                  fill_missing_gene_prop=None,
                  transform_gene_prop_func=None):
    """
    Prepare metabolic reaction network gene property dataframe

    Parameters
    ----------
    gene_prop_df : dataframe
        A dataframe of gene properties, such as log fold changes, normalized
        read counts, and z-scores. The 'gene' column of the dataframe is a
        list of gene symbols. Other columns of the dataframe are the gene
        properties of samples.
    rxn_gene_df : dataframe
        A dataframe of reactions and their associated genes. Following are the
        required columns: 'rxn' (reaction ID), 'rxn_name' (reaction name),
        'equation' (reaction equation), 'pathway' (reaction pathway), and
        'gene' (gene symbol). If a reaction is associated with multiple genes,
        one row lists one associated gene.
    fill_missing_gene_prop : int, or float, or None
        If int or float, replace missing gene properties with this value.
        If None, drop genes with missing properties.
    transform_gene_prop_func : function
        A function to transform each property column. If None, no
        transformation will be applied.

    Returns
    -------
    dataframe
    """
    if fill_missing_gene_prop is not None:
        all_genes = pd.concat(
            [gene_prop_df["gene"], rxn_gene_df["gene"]]
        ).drop_duplicates().tolist()

        gene_prop_df = (
            gene_prop_df
            .copy()
            .set_index("gene")
            .reindex(
                all_genes,
                fill_value=fill_missing_gene_prop)
            .reset_index()
        )

        assert not np.any(gene_prop_df.isnull().values)

    if transform_gene_prop_func is not None:
        gene_prop_df = transform_gene_prop(
            gene_prop_df, transform_gene_prop_func)

        assert not np.any(gene_prop_df.isnull().values)

    return gene_prop_df



def get_met_net_dfs(rxn_gene_df, rxn_edge_df, gene_prop_df, mn_weight_cutoff,
                    fill_missing_gene_prop=0,
                    transform_gene_prop_func=None,
                    rxn_gene_prop_agg_func=None):
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
        read counts, and z-scores. The 'gene' column of the dataframe is a
        list of gene symbols. Other columns of the dataframe are the gene
        properties of samples.
    mn_weight_cutoff : float
        The cutoff of metabolic network weights. Keep only edges with weights >
        mn_weight_cutoff.
    fill_missing_gene_prop : int, or float, or None
        If int or float, replace missing gene properties with this value.
        If None, drop genes with missing properties.
    transform_gene_prop_func : function
        A function to transform each property column.
    rxn_gene_prop_agg_func : function or None
        The function used to aggregate the property values of multiple genes
        that are associated with each reaction. If None, use
        lambda x: x.mean() to aggregate.

    Returns
    -------
    (dataframe, dataframe)
        A 2-tuple of a reaction node property dataframe and a reaction edge
        dataframe
    """
    if np.any(gene_prop_df.isnull().values):
        raise ValueError("gene_prop_df contains one or more NA/NaN... values")
    if np.any(rxn_edge_df.isnull().values):
        raise ValueError("rxn_edge_df contains one or more NA/NaN... values")
    if np.any(rxn_gene_df.isnull().values):
        raise ValueError("rxn_gene_df contains one or more NA/NaN... values")

    if rxn_edge_df[["src", "dest"]].duplicated().any():
        raise ValueError("rxn_edge_df contains duplicated (src, dest) pairs")

    if not all(rxn_edge_df.src <= rxn_edge_df.dest):
        raise ValueError(
            "rxn_edge_df requires src <= dest in alphabetical order")

    if rxn_gene_df[["rxn", "gene"]].duplicated().any():
        raise ValueError("rxn_gene_df contains duplicated (rxn, gene) pairs")

    if gene_prop_df["gene"].duplicated().any():
        raise ValueError("gene_prop_df contains duplicated genes.")

    gene_prop_df = get_mrn_gp_df(
        gene_prop_df, rxn_gene_df,
        fill_missing_gene_prop=fill_missing_gene_prop,
        transform_gene_prop_func=transform_gene_prop_func)

    rxn_pw_set_dict = (
        rxn_gene_df
        .groupby("rxn")["pathway"]
        .apply(set)
        .to_dict()
    )
    rxn_pw_dict = {}
    for k, v in rxn_pw_set_dict.items():
        if len(v) != 1:
            raise ValueError(
                "rxn_gene_df contains one reaction to "
                "multiple pathways mappings.")
        rxn_pw_dict[k] = list(v)[0]

    pw_rxn_set_dict = (
        rxn_gene_df
        .groupby("pathway")["rxn"]
        .apply(set)
        .to_dict()
    )

    rxn_edge_df = rxn_edge_df.loc[
        rxn_edge_df["mn_weight"] > mn_weight_cutoff, :].copy()

    edge_rxn_set = set(
        rxn_edge_df["src"].tolist() + rxn_edge_df["dest"].tolist())

    rxn_gene_prop_df = (
        rxn_gene_df
        .loc[rxn_gene_df["rxn"].isin(edge_rxn_set), :]
        .merge(
            gene_prop_df, how="left", on="gene",
            validate="many_to_one")
    )

    if fill_missing_gene_prop is None:
        rxn_gene_prop_df = rxn_gene_prop_df.dropna()
    else:
        assert not np.any(rxn_gene_prop_df.isnull().values)

    # weight filtered; property exists
    wf_pe_rxn_set = set(rxn_gene_prop_df["rxn"].tolist())
    rxn_edge_df = rxn_edge_df.loc[
        np.logical_and(
            rxn_edge_df["src"].isin(wf_pe_rxn_set),
            rxn_edge_df["dest"].isin(wf_pe_rxn_set)),
        :].copy()

    # If a rxn in a pathway cannot reach other pathways, remove the rxn.
    # Such rxns will have all random walk steps within their own pathways
    # regardless of the property values.
    f_graph = nx.from_pandas_edgelist(
        rxn_edge_df,
        source="src",
        target="dest")
    assert not f_graph.is_directed()
    # rxn to reachable rxn
    r_rrs_dict = {
        k: nx.node_connected_component(f_graph, k)
        for k in list(f_graph.nodes())
    }

    rm_rxns = set([
        k for k, v in r_rrs_dict.items()
        if len(v - pw_rxn_set_dict[rxn_pw_dict[k]]) <= 0
    ])

    # If a reaction in a pathway can reach other pathways,
    # the reaction cannot reach any reaction that cannot.
    rxn_edge_df = rxn_edge_df.loc[
        np.logical_not(np.logical_or(
            rxn_edge_df["src"].isin(rm_rxns),
            rxn_edge_df["dest"].isin(rm_rxns))),
        :].copy()

    def unique_one(x):
        x_set = set(x.tolist())
        assert len(x_set) == 1, str(x_set)
        return list(x_set)[0]

    agg_func_dict = {
        "rxn_name": unique_one,
        "equation": unique_one,
        "pathway": unique_one,
        "gene": (
            lambda x: ";;;".join(sorted(x.tolist()))
        )
    }

    if rxn_gene_prop_agg_func is None:
        rxn_gene_prop_agg_func = lambda x: x.mean()

    for i in rxn_gene_prop_df.columns:
        if i != "rxn" and i not in agg_func_dict:
            agg_func_dict[i] = rxn_gene_prop_agg_func

    rxn_prop_df = (
        rxn_gene_prop_df
        .groupby("rxn")
        .agg(agg_func_dict)
        .reset_index(names=["rxn"])
    )
    assert rxn_prop_df.isnull().values.sum() == 0

    common_rxns = (
        set(rxn_prop_df["rxn"].tolist())
        .intersection(set(rxn_edge_df["src"].tolist()))
        .intersection(set(rxn_edge_df["dest"].tolist()))
    )

    rxn_edge_df = rxn_edge_df.loc[
        np.logical_and(
            rxn_edge_df["src"].isin(common_rxns),
            rxn_edge_df["dest"].isin(common_rxns)),
        :].reset_index(drop=True).copy()

    rxn_prop_df = rxn_prop_df.loc[
        rxn_prop_df["rxn"].isin(common_rxns),
        :].reset_index(drop=True).copy()

    return rxn_prop_df, rxn_edge_df
