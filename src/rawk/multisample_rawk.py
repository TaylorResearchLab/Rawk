from .rawk_sample import RawkSample
from .rawk import Rawk, RawkTest

from joblib import Parallel, delayed
import json



class MultiSampleRawk:
    """
    Run Rawk on multiple samples

    Parameters
    ----------
    node_df : dataframe
        The reaction node dataframe with the following columns:
        'rxn', 'pathway', 'gene', 'sample1', 'sample2', 'sample3', ....
        Each value of the 'gene' column is a ;;; joined list of gene
        symbols.
    edge_df : dataframe
        The reaction edge dataframe with the following columns:
        'src', 'dest', 'mn_weight'.
    bg_sample_col : str or None
        The column in node_df used as background. If bg_sample_col is None, a
        uniform background sample will be created.
    n_jobs : int
        The number of jobs for running each sample parallely. Note that
        the number of CPUs for running each sample is set by
        n_workers_per_sample.
    seed : int
        Random number seed.
    node_emb_method : str or dict
        The method to set node embeddings, which can be one of the following
        options:
        - 'fastrp': Compute 2D t-SNE from fastrp node embeddings.
        - 'node2vec': Compute 2D t-SNE from node2vec node embeddings.
          This procedure may take several hours.
        - dict: Reuse pre-computed 2D embeddings of nodes stored in a dict
          with the following format:
          {node_name1: (dim1_coord, dim2_coord), ...}.
    n2v_walk_length : int
        node2vec walk length parameter.
    n2v_num_walks : int
        node2vec number of walks parameter.
    n2v_chunk_num_walks : int
        node2vec number of walks per chunk. Run node2vec in chunks to reduce
        mem usage.
    n2v_p : float
        node2vec p parameter.
    n2v_q : float
        node2vec q parameter.
    n2v_dimensions : int
        node2vec number of representation dimensions.
    n2v_window : int
        node2vec window parameter.
    n2v_epochs : int
        node2vec number of epochs parameter.
    random_walk_n_steps_method : str
        Generate random walk number of steps with one of the following methods:
        - 'simulation': Simulate random walks. Aggregate number of steps.
        - 'calculation': Calculate number of steps using transition probabity
          matrix. NOTE: This method is only implemented for n2v_p=1 and n2v_q=1.
          If node_emb_method='node2vec', the random walks still need to be
          simulated.
        - None: When n2v_p=1 and n2v_q=1, use the calculation method. Otherwise,
          use  the simulation method.
    n_workers_per_sample : int
        The number of workers to run each sample.
    verbose : bool
        Print runtime messages or not.
    """

    def __init__(self, node_df, edge_df, bg_sample_col=None,
                 n_jobs=1, seed=17, node_emb_method="fastrp",
                 n2v_walk_length=20,
                 n2v_num_walks=8000, n2v_chunk_num_walks=1000,
                 n2v_p=1.0, n2v_q=1.0, n2v_dimensions=50,
                 n2v_window=10, n2v_epochs=1,
                 random_walk_n_steps_method=None, n_workers_per_sample=1,
                 verbose=False):

        node_df = node_df.copy()
        edge_df = edge_df.copy()

        node_cols = ['rxn', 'pathway', 'gene']
        all_sample_cols = [x for x in node_df.columns.tolist()
                           if x not in node_cols]
        if len(all_sample_cols) == 0:
            raise ValueError("node_df should have > 0 sample column")

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
            rawk_sample = RawkSample(
                node_df, edge_df.copy(), name,
                workers=n_workers_per_sample,
                seed=seed, node_emb_method=node_emb_method,
                n2v_walk_length=n2v_walk_length,
                n2v_num_walks=n2v_num_walks,
                n2v_chunk_num_walks=n2v_chunk_num_walks,
                n2v_p=n2v_p, n2v_q=n2v_q, n2v_dimensions=n2v_dimensions,
                n2v_window=n2v_window,
                n2v_epochs=n2v_epochs,
                random_walk_n_steps_method=random_walk_n_steps_method,
                verbose=verbose)
            return (rawk_sample, sample_type)

        rawk_sample_arg_list = [
            {
                "node_df": bg_node_df,
                "name": bg_sample_col,
                "sample_type": "background"
            }
        ]
        for i in sample_cols:
            i_df = node_df.loc[:, node_cols + [i]].copy()
            i_df.columns = node_cols + ["property"]
            rawk_sample_arg_list.append(
                {
                    "node_df": i_df,
                    "name": i,
                    "sample_type": "sample"
                }
            )

        if n_jobs == 1:
            rawk_sample_list = [
                create_rawk_sample(**x)
                for x in rawk_sample_arg_list
            ]
        else:
            rawk_sample_list = Parallel(n_jobs=n_jobs, verbose=10)(
                delayed(create_rawk_sample)(**x)
                for x in rawk_sample_arg_list
            )

        rawk_list = []
        bg_sample = [x for x in rawk_sample_list if x[1] == "background"]
        assert len(bg_sample) == 1
        bg_sample = bg_sample[0][0]
        for i in rawk_sample_list:
            if i[1] == "sample":
                rawk_list.append(Rawk(i[0], bg_sample))

        self.rawk_sample_list = rawk_sample_list
        self.test = MultiSampleRawkTest(
            [x.rawk_test for x in rawk_list])
        self.rawk_list = rawk_list
        self.verbose = verbose

        self.node_df = node_df
        self.edge_df = edge_df
        self.node_emb_method = node_emb_method
        self.n2v_walk_length = n2v_walk_length
        self.n2v_num_walks = n2v_num_walks
        self.n2v_p = n2v_p
        self.n2v_q = n2v_q
        return


    def set_node_emb(self, node_emb_method,
                     pca_ndim=None, kmeans_n_clust=None):
        """
        Set node embeddings

        Parameters
        ----------
        node_emb_method : str or dict
            The method to set node embeddings, which can be one of the
            following options:
            - 'fastrp': Compute 2D t-SNE from fastrp node embeddings.
            - 'node2vec': Compute 2D t-SNE from node2vec node embeddings.
              This procedure may take several hours.
            - dict: Reuse pre-computed 2D embeddings of nodes stored in a dict
              with the following format:
              {node_name1: (dim1_coord, dim2_coord), ...}. In this option,
              pca_ndim and kmeans_n_clust should be None.
        pca_ndim : int or None
            When node_emb_method in ['fastrp', 'node2vec'], the number of PCA
            dimensions of the fastrp or node2vec embeddings to compute 2D t-SNE
            and kmeans clusters. If pca_ndim is None, default pca_ndim=6.
        kmeans_n_clust : int or None
            When node_emb_method in ['fastrp', 'node2vec'], the number of
            clusters parameter to run kmeans clustering. If None, default
            pca_ndim=10.

        Returns
        -------
        None
        """
        for r in self.rawk_list:
            r.set_node_emb(
                node_emb_method=node_emb_method,
                pca_ndim=pca_ndim, kmeans_n_clust=kmeans_n_clust)


    def test_num_steps(self, h_exponent=1, cmp_norm_fac=5000, pw_subset=None):
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
        pw_subset : set, or list, or None
            The subset of pathways to keep in the result dataframes.
            If is None, keep all pathways with not-NA results.

        Returns
        -------
        A 2-tuple of (fdr_dataframe, enrichment_score_dataframe).
        """
        return self.test.test_num_steps(
            h_exponent=h_exponent, cmp_norm_fac=cmp_norm_fac,
            pw_subset=pw_subset)


    def test_property_values(self, rw_s_prop_cutoff=0.005, s_ppd_dict=None,
                             pw_subset=None):
        """
        Enrichment test using the node property values

        Parameters
        ----------
        rw_s_prop_cutoff : float
            The cutoff proportion of random walk steps used for defininig local
            pathways.
        s_ppd_dict : dict or None
            {sample1: {pathway1: [property_value1, ...], ...}, ...} dict used
            for extracting pathway property values for testing. If None, use
            the foreground sample node properties.
        pw_subset : set, or list, or None
            The subset of pathways to keep in the result dataframes.
            If is None, keep all pathways with not-NA results.

        Returns
        -------
        A 2-tuple of (fdr_dataframe, enrichment_score_dataframe).
        """
        return self.test.test_property_values(
            rw_s_prop_cutoff=rw_s_prop_cutoff,
            s_ppd_dict=s_ppd_dict,
            pw_subset=pw_subset)



class MultiSampleRawkTest:
    """
    Run RawkTest on multiple samples

    Parameters
    ----------
    rawk_test_list : list of RawkTest
    """

    def __init__(self, rawk_test_list):
        self.rt_list = rawk_test_list

    @staticmethod
    def merge_test_res(test_res_list):
        """
        Merge a list of test results generated by Rawk.test_res_dict_to_df
        """
        adj_pval_df = test_res_list[0][0]
        es_df = test_res_list[0][1]
        for i in test_res_list[1:]:
            adj_pval_df = adj_pval_df.merge(i[0], how="inner", on="pathway")
            es_df = es_df.merge(i[1], how="inner", on="pathway")

        return adj_pval_df, es_df

    def test_num_steps(self, h_exponent=1, cmp_norm_fac=5000, pw_subset=None):
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

        pw_subset : set, or list, or None
            The subset of pathways to keep in the result dataframes.
            If is None, keep all pathways with not-NA results.

        Returns
        -------
        A 2-tuple of (fdr_dataframe, enrichment_score_dataframe).
        """
        test_res_list = []
        for i in self.rt_list:
            test_res_list.append(
                i.test_res_dict_to_df(
                    i.test_num_steps(
                        h_exponent=h_exponent,
                        cmp_norm_fac=cmp_norm_fac),
                    "_rawk", pw_subset))

        return self.merge_test_res(test_res_list)

    def test_property_values(self, rw_s_prop_cutoff=0.005, s_ppd_dict=None,
                             pw_subset=None):
        """
        Enrichment test using the node property values

        Parameters
        ----------
        rw_s_prop_cutoff : float
            The cutoff proportion of random walk steps used for defininig local
            pathways.
        s_ppd_dict : dict or None
            {sample1: {pathway1: [property_value1, ...], ...}, ...} dict used
            for extracting pathway property values for testing. If None, use
            the foreground sample node properties.

        pw_subset : set, or list, or None
            The subset of pathways to keep in the result dataframes.
            If is None, keep all pathways with not-NA results.

        Returns
        -------
        A 2-tuple of (fdr_dataframe, enrichment_score_dataframe).
        """
        test_res_list = []
        for i in self.rt_list:
            if s_ppd_dict is None:
                pw_pvs_dict = None
            else:
                pw_pvs_dict = s_ppd_dict[i.sample.name]

            test_res_list.append(
                i.test_res_dict_to_df(
                    i.test_property_values(
                        rw_s_prop_cutoff=rw_s_prop_cutoff,
                        pw_pvs_dict=pw_pvs_dict),
                    "_tlw", pw_subset))

        return self.merge_test_res(test_res_list)

    def to_list(self):
        """
        Convert data to a list
        """
        return [x.to_dict() for x in self.rt_list]

    @classmethod
    def from_list(cls, dl):
        """
        Construct a MultiSampleRawkTest from a list
        """
        return cls([RawkTest.from_dict(x) for x in dl])

    def save_to_json(self, file_path):
        """
        Save data to a json file

        Parameters
        ----------
        file_path : str
            Save json file path.

        Returns
        -------
        None
        """
        dl = self.to_list()

        with open(file_path, 'w') as f:
            json.dump(dl, f, indent=4)

        return

    @classmethod
    def load_json(cls, file_path):
        """
        Read json files

        Parameters
        ----------
        file_path : str
            Read file prefix.

        Returns
        -------
        A MultiSampleRawkTest constructed from the read data.
        """
        with open(file_path, 'r') as f:
            jl = json.load(f)

        return cls.from_list(jl)
