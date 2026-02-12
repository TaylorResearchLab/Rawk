from .rawk_sample import RawkSample
from .rawk import Rawk

from joblib import Parallel, delayed



class MultiSampleRawk:
    """
    Run Rawk on multiple samples

    Parameters
    ----------
    node_df : dataframe
        Required columns: 'rxn', 'pathway', 'gene', 'sample1', 'sample2',
        'sample3', .... A 'gene' column value is a ;;; joined list of gene
        symbols.
    edge_df : dataframe
        Required columns: 'src', 'dest', 'mn_weight'.
    bg_sample_col : str or None
        The column in node_df used as background. If bg_sample_col is None, a
        uniform background sample will be created.
    workers : int
        The number of workers for running each sample parallely. Note that
        each sample is analyzed using one worker.
    seed : int
        Random number seed.
    learn_emb : str or dict
        Learn node embeddings via the following options:
        - 'fastrp': Learn TSNE embeddings from fastrp node embeddings.
        - 'node2vec': Learn TSNE embeddings from node2vec node representation.
          This procedure may take several hours.
        - 'spring_layout': Learn node spring layout using
          networkx.spring_layout.
        - dict: A dictionary that maps node IDs to pre-computed 2D embedding
          coordinates, i.e., {node1: (dim1, dim2), node2: (dim1, dim2), ...}.
    n2v_walk_length : int
        node2vec walk length parameter.
    n2v_num_walks : int
        node2vec number of walks parameter.
    n2v_p : float
        node2vec p parameter.
    n2v_q : float
        node2vec q parameter.
    n2v_dimensions : int
        node2vec number of representation dimensions.
    n2v_window : int
        node2vec window parameter.
    n2v_min_count : int
        node2vec min count parameter.
    n2v_epochs : int
        node2vec number of epochs parameter.
    random_walk_n_steps_method : str
        Generate random walk number of steps with one of the following methods:
        - 'simulation': Simulate random walks. Aggregate number of steps.
        - 'calculation': Calculate number of steps using transition probabity
          matrix. NOTE: This method is only implemented for n2v_p=1 and n2v_q=1.
          If learn_emb='node2vec', the random walks still need to be simulated.
        - None: When n2v_p=1 and n2v_q=1, use the calculation method. Otherwise,
          use  the simulation method.
    """

    def __init__(self, node_df, edge_df, bg_sample_col=None,
                 workers=1, seed=17, learn_emb="fastrp", n2v_walk_length=20,
                 n2v_num_walks=5000, n2v_p=1, n2v_q=1, n2v_dimensions=50,
                 n2v_window=10, n2v_min_count=1, n2v_epochs=10,
                 random_walk_n_steps_method=None):

        node_cols = ['rxn', 'pathway', 'gene']
        all_sample_cols = [x for x in node_df.columns.tolist()
                           if x not in node_cols]

        if bg_sample_col is None:
            bg_node_df = node_df.loc[:, node_cols + all_sample_cols[:1]].copy()
            bg_node_df.columns = node_cols + ["property"]
            bg_node_df["property"] = 1e6 / bg_node_df.shape[0]
            sample_cols = all_sample_cols
            bg_sample_col = "uniform_background"
        else:
            bg_node_df = node_df.loc[:, node_cols + [bg_sample_col]].copy()
            bg_node_df.columns = node_cols + ["property"]
            sample_cols = [x for x in all_sample_cols if x != bg_sample_col]

        if len(sample_cols) == 0:
            raise ValueError("sample_columns should be non-empty")

        def create_rawk_sample(node_df, name, sample_type):
            mnea_sample = RawkSample(
                node_df, edge_df.copy(), name, workers=1,
                seed=seed, learn_emb=learn_emb,
                n2v_walk_length=n2v_walk_length, n2v_num_walks=n2v_num_walks,
                n2v_p=n2v_p, n2v_q=n2v_q, n2v_dimensions=n2v_dimensions,
                n2v_window=n2v_window, n2v_min_count=n2v_min_count,
                n2v_epochs=n2v_epochs,
                random_walk_n_steps_method=random_walk_n_steps_method)
            # Empty n2v_walks to save memory space
            mnea_sample.n2v_walks = []
            return (mnea_sample, sample_type)

        mnea_sample_arg_list = [
            {
                "node_df": bg_node_df,
                "name": bg_sample_col,
                "sample_type": "background"
            }
        ]
        for i in sample_cols:
            i_df = node_df.loc[:, node_cols + [i]].copy()
            i_df.columns = node_cols + ["property"]
            mnea_sample_arg_list.append(
                {
                    "node_df": i_df,
                    "name": i,
                    "sample_type": "sample"
                }
            )

        if workers == 1:
            mnea_sample_list = [
                create_rawk_sample(**x)
                for x in mnea_sample_arg_list
            ]
        else:
            mnea_sample_list = Parallel(n_jobs=workers, verbose=10)(
                delayed(create_rawk_sample)(**x)
                for x in mnea_sample_arg_list
            )

        mnea_list = []
        bg_sample = [x for x in mnea_sample_list if x[1] == "background"]
        assert len(bg_sample) == 1
        bg_sample = bg_sample[0][0]
        for i in mnea_sample_list:
            if i[1] == "sample":
                mnea_list.append(Rawk(i[0], bg_sample))

        self.mnea_sample_list = mnea_sample_list
        self.mnea_list = mnea_list
        return

    def test_num_steps(self, h_exponent=1, cmp_norm_fac=300):
        """
        Enrichment test using the number of random walk steps

        Parameters
        ----------
        h_exponent : float
            The heuristic factor for exponentiating sample vs background
            fold changes of the number of random walk steps. This heuristic
            factor intends to distinguish pathways with relatively small and
            large effect sizes.

        cmp_norm_fac : float
            Normaliziation factor for the number of random walk steps.

        Returns
        -------
        A 2-tuple of (fdr_dataframe, enrichment_score_dataframe).
        """
        test_res_list = []
        for i in self.mnea_list:
            test_res_list.append(
                i.test_res_dict_to_df(
                    i.test_num_steps(
                        h_exponent=h_exponent, cmp_norm_fac=cmp_norm_fac),
                    "_rawk"))

        adj_pval_df = test_res_list[0][0]
        es_df = test_res_list[0][1]
        for i in test_res_list[1:]:
            adj_pval_df = adj_pval_df.merge(i[0], how="inner", on="pathway")
            es_df = es_df.merge(i[1], how="inner", on="pathway")

        return adj_pval_df, es_df

    def test_property_values(self, rw_s_prop_cutoff=0.1):
        """
        Enrichment test using the node property values

        Parameters
        ----------
        rw_s_prop_cutoff : float
            The cutoff proportion of random walk steps used for defininig local
            pathways.

        Returns
        -------
        A 2-tuple of (fdr_dataframe, enrichment_score_dataframe).
        """
        test_res_list = []
        for i in self.mnea_list:
            test_res_list.append(
                i.test_res_dict_to_df(
                    i.test_property_values(rw_s_prop_cutoff=rw_s_prop_cutoff),
                    "_m2_tlw"))

        adj_pval_df = test_res_list[0][0]
        es_df = test_res_list[0][1]
        for i in test_res_list[1:]:
            adj_pval_df = adj_pval_df.merge(i[0], how="inner", on="pathway")
            es_df = es_df.merge(i[1], how="inner", on="pathway")

        return adj_pval_df, es_df
