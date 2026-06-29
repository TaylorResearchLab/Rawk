from .rawk_sample import RawkSample

import pandas as pd
import numpy as np
import copy
import scipy.stats as sps
import json



class Rawk:
    """
    Run Rawk on one sample vs background

    Parameters
    ----------
    sample : RawkSample
        An RawkSample instance to test for enrichment.
    bg_sample : RawkSample
        An RawkSample instance used as background.
    """

    def __init__(self, sample, bg_sample):
        self.assert_rs_compatible(sample, bg_sample)

        self.sample = sample
        self.bg_sample = bg_sample

        self.n2v_walk_length = sample.n2v_walk_length
        self.n2v_num_walks = sample.n2v_num_walks
        self.n2v_spr = self.n2v_num_walks * self.n2v_walk_length

        self._contrast_pw_graph()

        self.rawk_test = RawkTest(
            self.sample.node_df,
            self.sample.rw_cspr_rtp_sc,
            self.bg_sample.rw_cspr_rtp_sc,
            self.n2v_num_walks,
            self.sample.name)

        assert self.n2v_spr == self.rawk_test.spr


    @classmethod
    def assert_graph_compatible(cls, left, right):
        """
        Assert that two networkx graphs are compatible with Rawk
        """
        left_node_set = set(left.nodes)
        left_edge_set = set(left.edges)

        assert len(left_node_set) == len(left.nodes)
        assert len(left_edge_set) == len(left.edges)

        right_node_set = set(right.nodes)
        right_edge_set = set(right.edges)

        assert len(right_node_set) == len(right.nodes)
        assert len(right_edge_set) == len(right.edges)

        assert left_node_set == right_node_set
        assert left_edge_set == right_edge_set


    @classmethod
    def assert_rs_compatible(cls, left, right):
        """
        Assert that two RawkSample instances are compatible with Rawk
        """
        assert left.n2v_num_walks == right.n2v_num_walks
        assert left.n2v_walk_length == right.n2v_walk_length
        cls.assert_graph_compatible(left.rxn_graph, right.rxn_graph)
        cls.assert_graph_compatible(left.pw_graph, right.pw_graph)


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
            and kmeans clusters. If None, default pca_ndim=6.
        kmeans_n_clust : int or None
            When node_emb_method in ['fastrp', 'node2vec'], the number of
            clusters parameter to run kmeans clustering. If None, default
            pca_ndim=10.

        Returns
        -------
        None
        """
        for rs in [self.sample, self.bg_sample]:
            rs.set_node_emb(
                node_emb_method=node_emb_method,
                pca_ndim=pca_ndim, kmeans_n_clust=kmeans_n_clust)

        self._contrast_pw_graph()


    @staticmethod
    def _diff_x_y(x, xlab, y, ylab):
        """
        Computed labeled x values - labeled y values
        """
        assert len(xlab) == len(set(xlab))
        assert all([i in ylab for i in xlab]), str(xlab) + str(ylab)
        x_dict = dict(zip(xlab, x))
        y_dict = dict(zip(ylab, y))
        r = [x_dict[i] - y_dict[i] for i in xlab]
        return r


    def _contrast_pw_graph(self):
        """
        Create sample vs background pathway graph
        """
        contrast_pw_graph = copy.deepcopy(self.sample.pw_graph)
        for k, v in contrast_pw_graph.edges.items():
            assert (
                self.bg_sample.pw_graph.edges[k]['pw_n2v_transition_prob'] >= 0
            )
            x_ev = self.sample.pw_graph.edges[k]
            y_ev = self.bg_sample.pw_graph.edges[k]
            v['pw_n2v_transition_prob_diff'] = (
                x_ev['pw_n2v_transition_prob'] - y_ev['pw_n2v_transition_prob']
            )
            v['rw_pw_spct_diff'] = x_ev['rw_pw_spct'] - y_ev['rw_pw_spct']

            v['n2v_transition_prob_diff'] = self._diff_x_y(
                x_ev['n2v_transition_prob'], x_ev['n2v_transition_src'],
                y_ev['n2v_transition_prob'], y_ev['n2v_transition_src']
            )
            v['rw_spct_diff'] = self._diff_x_y(
                x_ev['rw_spct'], x_ev['n2v_transition_src'],
                y_ev['rw_spct'], y_ev['n2v_transition_src']
            )

        self.contrast_pw_graph = contrast_pw_graph


    def test_num_steps(self, h_exponent=1, cmp_norm_fac=5000):
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
        A dictionary of testing results, with the following structure:
        {
            sample_id: {
                pathway1: {'test_nes': test_nes, 'test_pval': test_pval},
                ...
            }
        }
        """
        mnea_res_dict = self.rawk_test.test_num_steps(
            h_exponent=h_exponent,
            cmp_norm_fac=cmp_norm_fac)

        return mnea_res_dict

    def test_property_values(self, rw_s_prop_cutoff=0.005, pw_pvs_dict=None):
        """
        Enrichment test using the node property values

        Parameters
        ----------
        rw_s_prop_cutoff : float
            The cutoff proportion of random walk steps used for defininig local
            pathways.
        pw_pvs_dict : dict or None
            {pathway1: [property_value1, ...], ...} dict used for extracting
            pathway property values for testing. If None, use the foreground
            sample node properties.

        Returns
        -------
        A dictionary of testing results, with the following structure:
        {
            sample_id: {
                pathway1: {'test_nes': test_nes, 'test_pval': test_pval},
                ...
            }
        }
        """
        res = self.rawk_test.test_property_values(
            rw_s_prop_cutoff=rw_s_prop_cutoff,
            pw_pvs_dict=pw_pvs_dict)

        return res



class RawkTest:
    """
    Rawk test number of random walk steps

    Parameters
    ----------
    fg_node_df : pd.DataFrame
        A foreground dataframe of reactions. Following are the required columns:
        'rxn' (reaction ID), 'pathway' (reaction pathway),
        'property' (reaction property value),
        'gene' (genes associated with the reaction). The 'gene' column lists
        the genes associated with each reaction in a ;;; joined list of gene
        symbols.
    fg_df : pd.DataFrame
        A foreground (pathways, pathway_____reaction)
        dataframe. (i, j) is the number of steps of all
        random walks starting from j in pathway i.
    bg_df : pd.DataFrame
        A background (pathways, pathway_____reaction)
        dataframe. (i, j) is the number of steps of all
        random walks starting from j in pathway i.
    num_walks : int
        Number of random walks starting from each reaction.
    test_id : str
        ID for this test.
    """

    def __init__(self, fg_node_df, fg_df, bg_df, num_walks,
                 test_id):
        self.fg_node_df = fg_node_df
        self.spr = self.get_spr(fg_df, bg_df)
        self.fg_df = fg_df
        self.bg_df = bg_df
        self.num_walks = num_walks
        self.test_id = test_id

    @staticmethod
    def get_spr(fg_df, bg_df):
        assert fg_df.columns.tolist() == bg_df.columns.tolist()
        assert fg_df.index.tolist() == bg_df.index.tolist()

        fg_uniq_colsum = fg_df.sum(axis=0).unique()
        assert len(fg_uniq_colsum) == 1

        bg_uniq_colsum = bg_df.sum(axis=0).unique()
        assert len(bg_uniq_colsum) == 1

        assert fg_uniq_colsum[0] == bg_uniq_colsum[0]

        return fg_uniq_colsum[0]

    @staticmethod
    def get_pw_ns(xdf, ydf, pw):
        pw_rxns = [c for c in xdf.columns.tolist()
                   if c.split("_____")[0] == pw]
        x = xdf.loc[pw, pw_rxns].values
        y = ydf.loc[pw, pw_rxns].values
        return x, y

    @staticmethod
    def get_a_test_ns_tbls(x, y, norm_fac, total_ns, h_exponent):
        assert len(x.shape) == 1
        assert x.shape == y.shape
        ns_tbl = np.array([x.sum(), y.sum()])
        ns_tbl = np.array([ns_tbl, total_ns * len(x) - ns_tbl])

        n_ns_tbl = np.array([
            ns_tbl[0, 0] / ns_tbl[:, 0].sum() * norm_fac,
            ns_tbl[0, 1] / ns_tbl[:, 1].sum() * norm_fac
        ])
        n_ns_tbl = np.array([n_ns_tbl, norm_fac - n_ns_tbl])

        h_fac = (n_ns_tbl[0, 0] / n_ns_tbl[0, 1]) ** h_exponent
        n_ns_tbl[0, 0] = np.clip(n_ns_tbl[0, 1] * h_fac, 0, norm_fac)
        n_ns_tbl[1, 0] = norm_fac - n_ns_tbl[0, 0]
        return ns_tbl, n_ns_tbl

    def test_num_steps(self, h_exponent=1, cmp_norm_fac=5000):
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
        A dictionary of testing results, with the following structure:
        {
            sample_id: {
                pathway1: {'test_nes': test_nes, 'test_pval': test_pval},
                ...
            }
        }
        """
        n2v_spr = self.spr
        assert n2v_spr == self.get_spr(self.fg_df, self.bg_df)

        mnea_res_dict = {}
        sid = self.test_id
        mnea_res_dict[sid] = {}

        # (pw, pw_rxn)
        # col sums to n2v_spr
        x_df = self.fg_df
        y_df = self.bg_df
        for i in x_df.index.tolist():
            mnea_res_dict[sid][i] = {}

            tx, ty = self.get_pw_ns(x_df, y_df, i)

            assert tx.shape == ty.shape

            if len(tx) == 0:
                mnea_res_dict[sid][i]["test_pval"] = np.nan
                mnea_res_dict[sid][i]["test_nes"] = np.nan
            else:
                tx = tx - self.num_walks
                ty = ty - self.num_walks
                mnea_res_dict[sid][i]["tx"] = tx
                mnea_res_dict[sid][i]["ty"] = ty

                cmp_ns_tbls = self.get_a_test_ns_tbls(
                    tx, ty, cmp_norm_fac, n2v_spr,
                    h_exponent)

                test_res = sps.fisher_exact(
                    np.round(cmp_ns_tbls[1]),
                    alternative="two-sided")

                ap_nns_tbl = np.round(cmp_ns_tbls[1]) + 0.5

                mnea_res_dict[sid][i]["test_nes"] = np.log(
                    (ap_nns_tbl[0, 0] / ap_nns_tbl[1, 0])
                    / (ap_nns_tbl[0, 1] / ap_nns_tbl[1, 1])
                )

                mnea_res_dict[sid][i]["test_pval"] = test_res.pvalue

        return mnea_res_dict

    def test_property_values(self, rw_s_prop_cutoff=0.005,
                             pw_pvs_dict=None):
        """
        Enrichment test using the node property values

        Parameters
        ----------
        rw_s_prop_cutoff : float
            The cutoff proportion of random walk steps used for defininig local
            pathways.
        pw_pvs_dict : dict or None
            {pathway1: [property_value1, ...], ...} dict used for extracting
            pathway property values for testing. If None, use the foreground
            sample node properties.

        Returns
        -------
        A dictionary of testing results, with the following structure:
        {
            sample_id: {
                pathway1: {'test_nes': test_nes, 'test_pval': test_pval},
                ...
            }
        }
        """
        n2v_spr = self.spr
        assert n2v_spr == self.get_spr(self.fg_df, self.bg_df)

        mnea_res_dict = {}
        sid = self.test_id
        mnea_res_dict[sid] = {}

        if pw_pvs_dict is None:
            pw_pvs_dict = (
                self.fg_node_df
                .groupby("pathway")["property"]
                .apply(lambda x: x.values)
                .to_dict()
            )

        # Use background to define local neighborhood.
        ln_df = self.bg_df
        for i_pw in ln_df.index.tolist():
            mnea_res_dict[sid][i_pw] = {}
            tp_rxns = [
                xx for xx in ln_df.columns.tolist()
                if xx.split("_____")[0] == i_pw]

            # pw fraction series
            cmp_rsf = (
                ln_df.loc[:, tp_rxns].sum(axis=1)
                / (len(tp_rxns) * n2v_spr)
            )

            # greater than cutoff
            gc_pws = cmp_rsf[
                cmp_rsf > rw_s_prop_cutoff].index.tolist()

            l_pws = [cr for cr in gc_pws if cr != i_pw]

            mnea_res_dict[sid][i_pw]["rows"] = l_pws

            if len(l_pws) == 0:
                # no local neighbor pw to compare to.
                mnea_res_dict[sid][i_pw]["test_pval"] = np.nan
                mnea_res_dict[sid][i_pw]["test_nes"] = np.nan
            else:
                etx = pw_pvs_dict[i_pw]
                elx = np.concatenate(
                    [pw_pvs_dict[x] for x in l_pws])

                test_res = sps.ranksums(etx, elx)

                mnea_res_dict[sid][i_pw]["test_nes"] = np.log(
                    (np.mean(etx) + 1) / (np.mean(elx) + 1)
                )

                mnea_res_dict[sid][i_pw]["test_pval"] = test_res.pvalue

        return mnea_res_dict

    @staticmethod
    def test_res_dict_to_df(mnea_res_dict, m_sfx, pw_subset=None):
        """
        Convert a test result dictionary to a 2-tuple of dataframes

        Parameters
        ----------
        mnea_res_dict : dict
            Rawk test result dictionary.

        m_sfx : str
            Method suffix for dataframe columns.

        pw_subset : set, or list, or None
            The subset of pathways to keep in the result dataframes.
            If is None, keep all pathways with not-NA results.

        Returns
        -------
        A 2-tuple of (fdr_dataframe, enrichment_score_dataframe).
        """
        if len(mnea_res_dict) != 1:
            raise ValueError(
                "Input mnea_res_dict should"
                "only have one sample")

        # k is sample_id
        # v is {pw1: {nes, pval}, ...}
        k = list(mnea_res_dict.keys())[0]
        v = mnea_res_dict[k]
        if pw_subset is None:
            pw_subset = set(v.keys())

        kcol = k + m_sfx
        mnea_padj_df = pd.DataFrame(
            [(pw, r["test_pval"])
             for pw, r in v.items()
             if pw in pw_subset],
            columns=["pathway", kcol]).dropna()

        mnea_padj_df.loc[:, kcol] = sps.false_discovery_control(
            mnea_padj_df.loc[:, kcol].values)

        # nes df
        mnea_nes_df = pd.DataFrame(
            [(pw, r["test_nes"])
             for pw, r in v.items()
             if pw in pw_subset],
            columns=["pathway", k + m_sfx]).dropna()

        return mnea_padj_df, mnea_nes_df

    def to_dict(self):
        """
        Convert data to a dict
        """
        assert self.spr == self.get_spr(
            self.fg_df, self.bg_df)

        dd = {
            'num_walks': self.num_walks,
            'test_id': self.test_id,
            'fg_node_df': self.fg_node_df.to_dict(
                orient='split'),
            'fg_df': self.fg_df.to_dict(orient='split'),
            'bg_df': self.bg_df.to_dict(orient='split'),
        }

        return dd

    @staticmethod
    def df_dict_to_df(jd):
        """
        Construct a dataframe from a dict
        """
        df = pd.DataFrame(
            data=jd['data'],
            index=jd['index'],
            columns=jd['columns'])
        return df

    @classmethod
    def from_dict(cls, dd):
        """
        Construct a RawkTest from a dict
        """
        return cls(
            cls.df_dict_to_df(dd['fg_node_df']),
            cls.df_dict_to_df(dd['fg_df']),
            cls.df_dict_to_df(dd['bg_df']),
            dd['num_walks'],
            dd['test_id'])

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
        dd = self.to_dict()

        with open(file_path, 'w') as f:
            json.dump(dd, f, indent=4)

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
        A RawkTest constructed from the read data.
        """
        with open(file_path, 'r') as f:
            jl = json.load(f)

        return cls.from_dict(jl)
