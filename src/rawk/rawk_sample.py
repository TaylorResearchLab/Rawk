import networkx as nx
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import numpy as np
from collections import Counter
import time
from gensim.models import Word2Vec


from .fastrp import fastrp_wrapper
from .n2v import RawkNode2Vec



class RawkSample:
    """
    Create a Rawk sample

    Parameters
    ----------
    node_df : dataframe
        A dataframe of reactions. Following are the required columns:
        'rxn' (reaction ID), 'pathway' (reaction pathway),
        'property' (reaction property value),
        'gene' (genes associated with the reaction). The 'gene' column lists
        the genes associated with each reaction in a ;;; joined list of gene
        symbols.
    edge_df : dataframe
        A dataframe of reactions edges. Following are the required columns:
        'src' (source node reaction ID), 'dest' (destination node reaction ID),
        'mn_weight' (metabolic network edge weight). The edges are undirected,
        with the src <= dest in alphabetical order.
    name : str
        The name of the sample.
    uniform_property : bool
        Whether the property values should be converted to the same value.
    workers : int
        The number of workers for running node2vec.
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
          If node_emb_method='node2vec', the random walks still need to be simulated.
        - None: When n2v_p=1 and n2v_q=1, use the calculation method. Otherwise,
          use  the simulation method.
    verbose : bool
        Print runtime messages or not.
    """

    def __init__(self, node_df, edge_df, name="sample", uniform_property=False,
                 workers=1, seed=17, node_emb_method="fastrp",
                 n2v_walk_length=20, n2v_num_walks=8000,
                 n2v_chunk_num_walks=1000,
                 n2v_p=1.0, n2v_q=1.0, n2v_dimensions=50,
                 n2v_window=10, n2v_epochs=1,
                 random_walk_n_steps_method=None,
                 keep_walks=False, verbose=False):

        if verbose:
            print("Running RawkSample {}...".format(name))

        step_start_time = time.perf_counter()

        if uniform_property:
            node_df = node_df.copy()
            node_df["property"] = 1e6 / node_df.shape[0]

        if not isinstance(name, str):
            raise ValueError("name parameter needs to be a string.")

        self.name = name

        # to be compatible with json
        # np.int64 is not compatible
        n2v_walk_length = int(n2v_walk_length)
        n2v_num_walks = int(n2v_num_walks)
        n2v_p = float(n2v_p)
        n2v_q = float(n2v_q)

        self.node_df = node_df
        self.edge_df = edge_df
        self.n2v_walk_length = n2v_walk_length
        self.n2v_num_walks = n2v_num_walks
        self.n2v_p = n2v_p
        self.n2v_q = n2v_q
        self.keep_walks = keep_walks
        self.verbose = verbose

        if random_walk_n_steps_method is None:
            if (n2v_p == 1 and n2v_q == 1) and (node_emb_method != "node2vec"):
                random_walk_n_steps_method = "calculation"
            else:
                random_walk_n_steps_method = "simulation"

        if random_walk_n_steps_method == "calculation":
            if not (n2v_p == 1 and n2v_q == 1):
                raise ValueError(
                    "random_walk_n_steps_method='calculation' is only"
                    "implemented for n2v_p == 1 and n2v_q == 1.")

            if node_emb_method == "node2vec":
                raise ValueError(
                    "node_emb_method='node2vec' requires "
                    "random_walk_n_steps_method='simulation'.")
        elif random_walk_n_steps_method != "simulation":
            raise ValueError(
                "Unknown random_walk_n_steps_method {}.".format(
                    random_walk_n_steps_method))

        self.random_walk_n_steps_method = random_walk_n_steps_method

        self.seed = seed

        self._init_n2v_edge_weight(node_df, edge_df)
        self._init_rxn_graph()
        self._init_pw_rxn_mappings()

        step_time = time.perf_counter() - step_start_time
        if verbose:
            print("Data init. took {:.2f} seconds.".format(step_time))

        step_start_time = time.perf_counter()

        self._run_n2v(
            workers=workers, seed=seed, node_emb_method=node_emb_method,
            n2v_walk_length=n2v_walk_length, n2v_num_walks=n2v_num_walks,
            n2v_chunk_num_walks=n2v_chunk_num_walks,
            n2v_p=n2v_p, n2v_q=n2v_q, n2v_dimensions=n2v_dimensions,
            n2v_window=n2v_window,
            n2v_epochs=n2v_epochs,
            random_walk_n_steps_method=random_walk_n_steps_method,
            keep_walks=keep_walks)

        self._init_rw_sc_mtx()

        self._init_pw_graph()

        step_time = time.perf_counter() - step_start_time
        if verbose:
            print("Random walks took {:.2f} seconds.".format(step_time))

        step_start_time = time.perf_counter()

        self.set_node_emb(node_emb_method)

        step_time = time.perf_counter() - step_start_time
        if verbose:
            print(
                "Node embedding took {:.2f} seconds.".format(
                    step_time))
            print("RawkSample {} completed\n------".format(name))

    def _init_n2v_edge_weight(self, node_df, edge_df):
        assert node_df.notna().all().all()
        assert edge_df.notna().all().all()

        undirected_pairs = pd.DataFrame({
            "low": edge_df[["src", "dest"]].min(axis=1),
            "high": edge_df[["src", "dest"]].max(axis=1)
        })
        assert not undirected_pairs.duplicated().any()

        reverse_edge_df = (
            edge_df
            .loc[edge_df["src"] != edge_df["dest"], :]
            .copy()
            .rename(columns={"src": "dest", "dest": "src"})
        )

        directed_df = pd.concat(
            [edge_df.copy(), reverse_edge_df],
            ignore_index=True)

        assert not directed_df.duplicated(subset=["src", "dest"]).any()
        assert directed_df.notna().all().all()

        assert not node_df.duplicated(subset="rxn").any()
        property_by_node = dict(zip(node_df["rxn"], node_df["property"]))

        directed_df["dest_property"] = directed_df["dest"].map(
            property_by_node, na_action="ignore")

        assert directed_df.notna().all().all()

        dp = directed_df["dest_property"]

        assert dp.min() >= 0

        if dp.min() == 0:
            directed_df["n2v_weight"] = (
                directed_df["dest_property"]
                + dp[dp > 0].min() / np.e
            )
        else:
            directed_df["n2v_weight"] = directed_df["dest_property"]

        self.ea_df = directed_df

    def _init_rxn_graph(self):
        node_list = []
        for i, i_row in self.node_df.iterrows():
            i_row = i_row.to_dict()

            i_node = (
                i_row["rxn"],
                {
                    "pathway": i_row["pathway"],
                    "gene": tuple(
                        sorted(i_row["gene"].split(";;;"))),
                    "property": i_row["property"],
                    "pos": np.zeros(2)
                },
            )

            node_list.append(i_node)

        edge_list = []
        for i, i_row in self.ea_df.iterrows():
            i_row = i_row.to_dict()

            i_edge = (
                i_row["src"],
                i_row["dest"],
                {
                    "mn_weight": i_row["mn_weight"],
                    "n2v_weight": i_row["n2v_weight"]
                },
            )
            edge_list.append(i_edge)

        rxn_graph = nx.DiGraph()
        rxn_graph.add_nodes_from(node_list)
        rxn_graph.add_edges_from(edge_list)

        self.rxn_graph = rxn_graph
        self._add_n2v_prob_to_rxn_g()


    def _init_pw_rxn_mappings(self):
        pw_by_rxn = {}
        for n, nd in self.rxn_graph.nodes.items():
            assert n not in pw_by_rxn
            pw_by_rxn[n] = nd["pathway"]

        pr_by_r = {}
        r_by_pr = {}
        p_by_pr = {}
        for rxn, pw in pw_by_rxn.items():
            pr = pw + "_____" + rxn
            assert rxn not in pr_by_r
            pr_by_r[rxn] = pr

            assert pr not in r_by_pr
            r_by_pr[pr] = rxn

            assert pr not in p_by_pr
            p_by_pr[pr] = pw

        self.pw_by_rxn = pw_by_rxn
        self.pr_by_r = pr_by_r
        self.r_by_pr = r_by_pr
        self.p_by_pr = p_by_pr


    def _run_n2v(self, workers=1, seed=17, node_emb_method="node2vec",
                 n2v_walk_length=20, n2v_num_walks=5000,
                 n2v_chunk_num_walks=500,
                 n2v_p=1, n2v_q=1, n2v_dimensions=50,
                 n2v_window=10, n2v_epochs=10,
                 random_walk_n_steps_method=None,
                 keep_walks=False):

        if random_walk_n_steps_method == "simulation":
            assert (
                sorted(self.node_df["rxn"].tolist())
                == sorted(self.rxn_graph.nodes())
            )

            n2v = RawkNode2Vec(
                self.node_df, self.ea_df,
                p=n2v_p, q=n2v_q,
                num_walks=n2v_num_walks,
                chunk_num_walks=n2v_chunk_num_walks,
                walk_length=n2v_walk_length,
                workers=workers, seed=seed,
                keep_walks=keep_walks)

            if node_emb_method == "node2vec":
                model = Word2Vec(
                    sentences=n2v,
                    vector_size=n2v_dimensions,
                    window=n2v_window,
                    min_count=0,
                    sg=1,
                    workers=workers,
                    epochs=n2v_epochs,
                    seed=seed
                )

                n2v_emb_dict = {
                    n2v.idx_to_node[int(x)]: model.wv.get_vector(x)
                    for x in model.wv.key_to_index.keys()
                }

                n2v_emb_df = pd.DataFrame(
                    [[k] + v.tolist() for k, v in n2v_emb_dict.items()])
            else:
                for _ in n2v:
                    pass
                model = None
                n2v_emb_dict = None
                n2v_emb_df = None

            n2v_walks = n2v.walks

        else:
            assert node_emb_method != "node2vec"
            n2v = None
            n2v_walks = None
            model = None
            n2v_emb_dict = None
            n2v_emb_df = None

        self.n2v = n2v
        self.n2v_walks = n2v_walks
        self.model = model
        self.n2v_emb_dict = n2v_emb_dict
        self.n2v_emb_df = n2v_emb_df

    @staticmethod
    def run_kmeans(X, k, seed):
        km = KMeans(
            n_clusters=k, random_state=seed,
            n_init="auto")

        clust = km.fit_predict(X)

        return km, clust


    @staticmethod
    def most_frequent(x):
        xc = sorted(
            list(Counter(x).items()),
            key=lambda t: (-t[1], t[0]))
        return xc[0][0]


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
        if isinstance(node_emb_method, dict):
            assert pca_ndim is None
            assert kmeans_n_clust is None

            pos = {
                k: np.array(v[:2])
                for k, v in node_emb_method.items()
            }

            for i in self.rxn_graph.nodes:
                if i not in node_emb_method:
                    err_msg = (
                        "node_emb_method dict should contain "
                        "all node ids. {} is not in the "
                        "node_emb_method dict."
                    ).format(i)
                    raise ValueError(err_msg)

            emb_arr = None
            emb_keys = None
            pca = None
            pca_df = None
            tsne = None
            km = None
            rxn_clusters = None

        else:
            if node_emb_method == "node2vec":
                assert self.n2v_emb_df is not None
                emb_df = self.n2v_emb_df.copy()

                emb_keys = emb_df.iloc[:, 0].tolist()
                emb_arr = emb_df.iloc[:, 1:].values

            elif node_emb_method == "fastrp":
                fastrp_conf = {
                    'projection_method': 'sparse',
                    'input_matrix': 'adj',
                    'weights': [0.0, 1.0, 0.5, 0.25],
                    'normalization': True,
                    'dim': 32,
                    'alpha': -0.5
                }

                emb_keys = list(self.rxn_graph.nodes)

                adj_mat = nx.to_scipy_sparse_array(
                    self.rxn_graph, nodelist=emb_keys, weight="n2v_weight",
                    format="csc")

                emb_arr = fastrp_wrapper(adj_mat, fastrp_conf)

            else:
                raise ValueError(
                    "Unknown emb {}".format(self.emb))

            assert len(emb_keys) == len(set(emb_keys))

            pca = PCA(
                n_components=min(emb_arr.shape[1], 30),
                svd_solver="full", random_state=self.seed)

            pca_arr = pca.fit_transform(
                StandardScaler().fit_transform(emb_arr))

            pca_df = pd.DataFrame(pca_arr, index=emb_keys)

            tsne = TSNE(
                n_components=2, random_state=self.seed,
                perplexity=min(30, pca_arr.shape[0] - 1))

            if pca_ndim is None:
                pca_ndim = 6

            tsne_arr = tsne.fit_transform(
                pca_arr[:, :pca_ndim])

            pos = {
                emb_keys[i]: tsne_arr[i, :]
                for i in range(len(emb_keys))
            }

            if kmeans_n_clust is None:
                kmeans_n_clust = 10

            km, rxn_clusters = self.run_kmeans(
                pca_arr[:, :pca_ndim],
                kmeans_n_clust, self.seed)


        if km is None:
            assert rxn_clusters is None
            rxn_pw_clust_df = None
            pw_clust_df = None

        else:
            rxn_pw_clust_df = []
            for i in range(len(emb_keys)):
                rxn_pw_clust_df.append(
                    [emb_keys[i],
                     self.pw_by_rxn[emb_keys[i]],
                     str(rxn_clusters[i])])

            rxn_pw_clust_df = pd.DataFrame(
                sorted(rxn_pw_clust_df, key=lambda r: r[0]),
                columns=["rxn", "pw", "clust"])

            clust_by_rxn = dict(
                zip(rxn_pw_clust_df["rxn"], rxn_pw_clust_df["clust"]))

            pw_clust_df = rxn_pw_clust_df.groupby("pw")["clust"].agg(
                self.most_frequent).reset_index()

            assert pw_clust_df["pw"].is_unique

            clust_by_pw = dict(
                zip(pw_clust_df["pw"], pw_clust_df["clust"]))


        for rxn, rxn_d in self.rxn_graph.nodes.items():
            rxn_d["pos"] = pos[rxn]

            if km is None:
                rxn_d["cluster"] = None
            else:
                rxn_d["cluster"] = clust_by_rxn[rxn]

        for pw, pw_d in self.pw_graph.nodes.items():
            pw_d["rxn_pos"] = np.array(
                [pos[x] for x in pw_d["rxn"]])

            pw_d["pw_pos"] = np.mean(pw_d["rxn_pos"], axis=0)
            if km is None:
                pw_d["cluster"] = None
            else:
                pw_d["cluster"] = clust_by_pw[pw]

        self.emb_arr = emb_arr
        self.emb_keys = emb_keys

        self.pca = pca
        self.pca_df = pca_df
        self.pca_ndim = pca_ndim

        self.tsne = tsne
        self.rxn_pos = pos
        self.rxn_pos_df = pd.DataFrame(
            [(k, v[0], v[1]) for k, v in sorted(pos.items())],
            columns=["rxn", "dim1", "dim2"])

        self.km = km
        self.kmeans_n_clust = kmeans_n_clust
        self.rxn_pw_clust_df = rxn_pw_clust_df
        self.pw_clust_df = pw_clust_df

        self.node_emb_method = node_emb_method


    def _add_n2v_prob_to_rxn_g(self):
        rxn_graph = self.rxn_graph
        for node in rxn_graph.nodes():
            # Get all outgoing edges and their weights
            # If no weight, None will be returned, and the sum below
            # will be stopped.
            out_edges = rxn_graph.out_edges(node, data="n2v_weight")
            # All nodes should have possible next steps, so all random walks
            # will have the same length as specified, which is required
            # in the counting of visit steps.
            assert len(out_edges) > 0
            # Calculate total sum of outgoing weights
            total_weight = sum(weight for _, _, weight in out_edges)
            assert total_weight > 0
            for _, nbr, weight in out_edges:
                assert weight > 0
                # Normalize weights to create a probability distribution
                prob = weight / total_weight
                rxn_graph[node][nbr]["n2v_transition_prob"] = prob

    def save_n2v_res(self, pfx):
        if self.model is not None:
            self.model.save(pfx + ".model")
            self.model.wv.save_word2vec_format(pfx + ".emb")
            self.n2v_emb_df.to_csv(pfx + "_emb.csv")

        if self.rxn_pos_df is not None:
            self.rxn_pos_df.to_csv(pfx + "_rxn_pos.csv")

        if self.n2v_walks is not None:
            pd.DataFrame(self.n2v_walks, copy=False).to_csv(pfx + "_walks.csv")

        self.edge_attr_df.to_csv(pfx + "_edge_attr.csv")
        nx.write_gml(self.rxn_graph, pfx + "_graph.gml")

    def read_n2v_res(self, pfx):
        self.rxn_graph = nx.read_gml(pfx + "_graph.gml")

        try:
            self.n2v_walk_df = pd.read_csv(
                pfx + "_walks.csv", low_memory=False, index_col=0)
            self.n2v_walks = self.n2v_walk_df.values.tolist()
        except FileNotFoundError:
            self.n2v_walks = None

        try:
            self.rxn_pos_df = pd.read_csv(pfx + "_rxn_pos.csv", index_col=0)
            self.rxn_pos = dict(zip(
                self.rxn_pos_df.iloc[:, 0].values.tolist(),
                self.rxn_pos_df.iloc[:, 1:3].values.tolist()
            ))
        except FileNotFoundError:
            self.rxn_pos_df = None
            self.rxn_pos = None
        return rd

    def rxn_to_pw_rxn(self, df, index=True, columns=True):
        df = df.copy()

        if index:
            df = (
                df
                .rename(index=self.pr_by_r)
                .sort_index(axis=0)
            )

        if columns:
            df = (
                df
                .rename(columns=self.pr_by_r)
                .sort_index(axis=1)
            )

        return df.copy()


    @staticmethod
    def _calculate_rw_ns_mat(tpm, walk_length=20, num_walks=5000):
        p = np.diag([1] * tpm.shape[0])
        res = p * num_walks
        # tpm is (src, tgt) prob. tpm[:, 0] is prob to first node
        # p[src, tgt] is the count of walks starting from src
        # at tgt, so p[src, :] * tpm[:, src] is the count of walks
        # at the next step.
        for i in range(walk_length - 1):
            p = p @ tpm
            res = res + (p * num_walks)

        res = res.T

        # normalize (tgt, src) number of steps matrix
        # - col sums to num_walks * walk_length
        # - start src rxn needs to have >= num_walks
        #   steps to be compatible with simulation

        ns_per_rxn = num_walks * walk_length
        min_start_ns = num_walks

        n = res.shape[0]
        assert len(res.shape) == 2
        assert res.shape[0] == res.shape[1]
        assert not np.isnan(res).any()
        assert (res >= 0).all()

        for j in range(n):
            col_sum = np.sum(res[:, j])

            if col_sum == 0:
                res[:, j] = 0
                res[j, j] = ns_per_rxn
                continue

            scaled = res[:, j] * (ns_per_rxn / col_sum)

            # make res[j, j] >= min_start_ns
            if scaled[j] < min_start_ns:
                # scaled sums to ns_per_rxn
                # res sums to col_sum
                remaining_sum = col_sum - res[j, j]

                if remaining_sum > 0:
                    not_start_mask = np.ones(n, dtype=bool)
                    not_start_mask[j] = False

                    scaled[not_start_mask] = (
                        res[not_start_mask, j]
                        * ((ns_per_rxn - min_start_ns) / remaining_sum)
                    )
                else:
                    scaled[:] = 0

                scaled[j] = min_start_ns

            # make res sum to ns_per_rxn
            base = np.floor(scaled).astype(int)

            diff = ns_per_rxn - base.sum()
            if diff > 0:
                remainders = np.round(scaled - base, 12)

                if scaled[j] == min_start_ns:
                    # remainders all >= 0
                    remainders[j] = -1.0

                largest_remainder_idx = np.argsort(
                    remainders, kind="stable")[::-1]

                assert diff <= len(largest_remainder_idx)
                base[largest_remainder_idx[:diff]] += 1

            elif diff < 0:
                # add back at the end of this branch
                base[j] = base[j] - min_start_ns

                largest_base_idx = np.argsort(
                    base, kind="stable")[::-1]

                largest_g0_base_idx = largest_base_idx[
                    base[largest_base_idx] > 0]

                n_g0_base = len(largest_g0_base_idx)

                assert base[largest_g0_base_idx].sum() >= abs(diff)

                n_steps = abs(diff)

                for step in range(n_steps):
                    d_idx = largest_g0_base_idx[(step % n_g0_base)]
                    base[d_idx] -= 1
                    diff += 1
                    if (step + 1) % n_g0_base == 0:
                        # update after iter. over all g0 base idx
                        if diff == 0:
                            assert step == n_steps - 1
                            break
                        largest_g0_base_idx = largest_base_idx[
                            base[largest_base_idx] > 0]
                        n_g0_base = len(largest_g0_base_idx)
                        assert n_g0_base > 0

                base[j] = base[j] + min_start_ns

            res[:, j] = base

        assert not np.isnan(res).any()
        assert (res >= 0).all()
        assert (res.sum(axis=0) == ns_per_rxn).all()
        return res

    def _init_rw_sc_mtx(self):
        rxn_graph = self.rxn_graph

        rxns = sorted(rxn_graph.nodes)
        assert len(rxns) == len(set(rxns))

        if self.random_walk_n_steps_method == "simulation":
            # random walk
            # (row target rxn, col source rxn)
            # step count
            rw_csr_rtr_sc_arr = self.n2v.get_visit_arr()

            n2v_nodes = [
                self.n2v.idx_to_node[x]
                for x in range(self.n2v.num_nodes)]

            assert n2v_nodes == rxns

        elif self.random_walk_n_steps_method == "calculation":
            rw_csr_rtr_sc_arr = self._calculate_rw_ns_mat(
                # (src, tgt) probs
                # returns (tgt, src) step counts
                nx.to_scipy_sparse_array(
                    rxn_graph, nodelist=rxns,
                    weight="n2v_transition_prob",
                    format="csc").toarray(),
                self.n2v_walk_length, self.n2v_num_walks)

        else:
            raise ValueError(
                "Unknown random_walk_n_steps_method {}.".format(
                    self.random_walk_n_steps_method))

        rw_csr_rtr_sc = pd.DataFrame(
            rw_csr_rtr_sc_arr, index=rxns, columns=rxns)

        # pr is pathway_____reaction
        rw_cspr_rtpr_sc = self.rxn_to_pw_rxn(
            rw_csr_rtr_sc,
            index=True, columns=True)

        assert rw_cspr_rtpr_sc.index.is_unique
        assert rw_cspr_rtpr_sc.columns.is_unique
        assert (
            rw_cspr_rtpr_sc.index.tolist()
            == rw_cspr_rtpr_sc.columns.tolist()
        )

        # (tgt_pw sum, src_rxn) n steps
        rw_cspr_rtp_sc = (
            rw_cspr_rtpr_sc
            .copy()
            .groupby(self.p_by_pr)
            .sum()
        )

        rw_csr_rtp_sc = (
            rw_cspr_rtp_sc
            .copy()
            .rename(columns=self.r_by_pr)
            .sort_index(axis=1)
        )

        # random walk steps per rxn
        spr = self.n2v_num_walks * self.n2v_walk_length

        np.testing.assert_allclose(
            rw_cspr_rtp_sc.sum().values,
            spr, rtol=0, atol=1e-6)
        np.testing.assert_allclose(
            rw_csr_rtp_sc.sum().values,
            spr, rtol=0, atol=1e-6)

        # (tgt pathway sum, src pathway mean) steps %
        rw_csp_rtp_spct = (
            rw_cspr_rtp_sc
            .copy()
            .T
            .groupby(self.p_by_pr)
            .mean()
            .T
        )

        np.testing.assert_allclose(
            rw_csp_rtp_spct.sum().values,
            spr, rtol=0, atol=1e-6)

        rw_csp_rtp_spct = rw_csp_rtp_spct * 100 / spr

        self.rw_csr_rtr_sc = rw_csr_rtr_sc
        self.rw_cspr_rtpr_sc = rw_cspr_rtpr_sc

        self.rw_cspr_rtp_sc = rw_cspr_rtp_sc
        self.rw_csr_rtp_sc = rw_csr_rtp_sc

        self.rw_csp_rtp_spct = rw_csp_rtp_spct


    def _init_pw_graph(self):
        rxn_g = self.rxn_graph
        rw_csp_rtp_spct = self.rw_csp_rtp_spct
        rw_csr_rtp_sc = self.rw_csr_rtp_sc

        pw_node_dict = {}
        for k, v in rxn_g.nodes.items():
            pw = v["pathway"]
            if pw not in pw_node_dict:
                pw_node_dict[pw] = {
                    "rxn": [],
                    "property": [],
                    "rxn_pos": [],
                    "pathway": pw,
                    "gene": set()
                }

            pw_node_dict[pw]["rxn"].append(k)
            pw_node_dict[pw]["property"].append(v["property"])
            pw_node_dict[pw]["rxn_pos"].append(v["pos"])
            pw_node_dict[pw]["gene"].update(v["gene"])

        for v in pw_node_dict.values():
            v["rxn_pos"] = np.array(v["rxn_pos"])
            v["pw_pos"] = np.mean(v["rxn_pos"], axis=0)
            v["pw_property"] = np.mean(v["property"])
            v["gene"] = sorted(v["gene"])

        pw_edge_dict = {}
        for src, src_v in rxn_g.nodes.items():
            # accumulated prob from the src rxn to dest pw
            pw_p_dict = {}

            for dest in nx.neighbors(rxn_g, src):
                dest_v = rxn_g.nodes[dest]

                pw_e = (src_v["pathway"], dest_v["pathway"])
                if pw_e not in pw_p_dict:
                    pw_p_dict[pw_e] = 0

                pw_p_dict[pw_e] += rxn_g.edges[(src, dest)][
                    "n2v_transition_prob"]

                if pw_e not in pw_edge_dict:
                    pw_edge_dict[pw_e] = {
                        "rxn": [],
                        "n2v_weight": [],
                        "n2v_transition_src": [],
                        "n2v_transition_prob": [],
                    }

                pw_edge_dict[pw_e]["rxn"].append((src, dest))
                pw_edge_dict[pw_e]["n2v_weight"].append(
                    rxn_g.edges[(src, dest)]["n2v_weight"])

            # pw_p_dict = {
            #    (src_pw, tgt_pw): sum(prob_from_src_rxn_to_tgt_pw),
            #    ...}
            for k, v in pw_p_dict.items():
                pw_edge_dict[k]["n2v_transition_src"].append(src)
                pw_edge_dict[k]["n2v_transition_prob"].append(v)

        for k, v in pw_edge_dict.items():
            assert len(v["n2v_transition_prob"]) != 0
            v["pw_n2v_weight"] = np.mean(v["n2v_weight"])
            v["pw_n2v_transition_prob"] = np.mean(v["n2v_transition_prob"])
            v["dest_pw_property"] = pw_node_dict[k[1]]["pw_property"]
            v["dest_property"] = pw_node_dict[k[1]]["property"]
            v["rw_pw_spct"] = rw_csp_rtp_spct.loc[k[1], k[0]]
            v["rw_spct"] = [
                (
                    rw_csr_rtp_sc.loc[k[1], x]
                    / (self.n2v_num_walks * self.n2v_walk_length) * 100
                )
                for x in v["n2v_transition_src"]
            ]

        pw_graph = nx.DiGraph()
        pw_graph.add_nodes_from(
            [(k, v) for k, v in pw_node_dict.items()])
        pw_graph.add_edges_from(
            [(k[0], k[1], v) for k, v in pw_edge_dict.items()])

        self.pw_graph = pw_graph
