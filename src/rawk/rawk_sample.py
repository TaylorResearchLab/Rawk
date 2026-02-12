import networkx as nx
import pandas as pd
from sklearn.manifold import TSNE
from node2vec import Node2Vec
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import time
import warnings
import textwrap


from .fastrp import fastrp_wrapper
from .n2v import NoWalkNode2Vec



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
    learn_emb : str or dict
        Learn node embeddings with one of the following options:
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

    def __init__(self, node_df, edge_df, name="sample", uniform_property=False,
                 workers=1, seed=17, learn_emb="fastrp", n2v_walk_length=20,
                 n2v_num_walks=5000, n2v_p=1, n2v_q=1, n2v_dimensions=50,
                 n2v_window=10, n2v_min_count=1, n2v_epochs=10,
                 random_walk_n_steps_method=None):

        if uniform_property:
            node_df = node_df.copy()
            node_df["property"] = 1e6 / node_df.shape[0]

        if not isinstance(name, str):
            raise ValueError("name parameter needs to be a string.")

        self.name = name

        self.node_df = node_df
        self.edge_df = edge_df
        self.learn_emb = learn_emb
        self.n2v_walk_length = n2v_walk_length
        self.n2v_num_walks = n2v_num_walks

        if random_walk_n_steps_method is None:
            if n2v_p == 1 and n2v_q == 1:
                random_walk_n_steps_method = "calculation"
            else:
                random_walk_n_steps_method = "simulation"

        if random_walk_n_steps_method == "calculation":
            if not (n2v_p == 1 and n2v_q == 1):
                raise ValueError(
                    "random_walk_n_steps_method='calculation' is only"
                    "implemented for n2v_p == 1 and n2v_q == 1.")

        if random_walk_n_steps_method not in ["calculation", "simulation"]:
            raise ValueError(
                "Unknown random_walk_n_steps_method {}.".format(
                    random_walk_n_steps_method))

        self.random_walk_n_steps_method = random_walk_n_steps_method

        self.seed = seed

        self.ea_df = self._add_edge_df_attr(node_df, edge_df)

        self._run_n2v(
            workers=workers, seed=seed, learn_emb=learn_emb,
            n2v_walk_length=n2v_walk_length, n2v_num_walks=n2v_num_walks,
            n2v_p=n2v_p, n2v_q=n2v_q, n2v_dimensions=n2v_dimensions,
            n2v_window=n2v_window, n2v_min_count=n2v_min_count,
            n2v_epochs=n2v_epochs,
            random_walk_n_steps_method=random_walk_n_steps_method)

        self._add_n2v_emb_dim_red()
        self._add_n2v_prob_to_graph()
        self._get_n2v_mat_dict()
        self._rxn_to_pw_graph()

    def _add_edge_df_attr(self, node_df, edge_df):
        adf = node_df.copy()
        # Convert undirected edges to directed.
        rdf = edge_df.copy()
        rdf2 = edge_df.copy()
        rdf2[["src", "dest"]] = rdf2[["dest", "src"]]
        rdf2 = rdf2.loc[(rdf2["src"] != rdf2["dest"]).values, :]
        assert rdf.columns.tolist() == rdf2.columns.tolist()
        rdf = pd.concat([rdf, rdf2], ignore_index=True)
        assert rdf[["src", "dest"]].drop_duplicates().shape[0] == rdf.shape[0]

        rdf = rdf.join(
            adf.set_index("rxn").rename(columns=lambda x: "src_" + str(x)),
            on = "src")
        rdf = rdf.join(
            adf.set_index("rxn").rename(columns=lambda x: "dest_" + str(x)),
            on = "dest")

        rdf["norm_mn_weight"] = rdf["mn_weight"] * 1e6 / rdf["mn_weight"].sum()

        assert min(rdf["dest_property"]) >= 0
        if min(rdf["dest_property"]) == 0:
            rdf["norm_weight"] = (
                rdf["dest_property"]
                + min(rdf["dest_property"][rdf["dest_property"] > 0]) / 1e6
            )
        else:
            rdf["norm_weight"] = rdf["dest_property"]
        return rdf

    def _run_n2v(self, workers=1, seed=17, learn_emb="node2vec",
                 n2v_walk_length=20, n2v_num_walks=5000,
                 n2v_p=1, n2v_q=1, n2v_dimensions=50,
                 n2v_window=10, n2v_min_count=1, n2v_epochs=10,
                 random_walk_n_steps_method=None):
        node_list = []
        node_set = set()
        for i in range(self.ea_df.shape[0]):
            i_row = self.ea_df.iloc[i, :].to_dict()
            for j in ["src", "dest"]:
                j_key = i_row[j]
                if j_key not in node_set:
                    node_set.add(j_key)
                    k = (
                        j_key,
                        {
                            "pathway": i_row[j + "_pathway"],
                            "gene": tuple(
                                sorted(i_row[j + "_gene"].split(";;;"))),
                            "property": i_row[j + "_property"],
                        }
                    )
                    node_list.append(k)

        edge_list = []
        edge_set = set()
        for i in range(self.ea_df.shape[0]):
            i_row = self.ea_df.iloc[i, :].to_dict()
            i_key = (i_row["src"], i_row["dest"])
            if i_key not in edge_set:
                edge_set.add(i_key)
                # ({},) is a tuple
                edge = i_key + (
                    {
                        "mn_weight": i_row["mn_weight"],
                        "norm_weight": i_row["norm_weight"]
                    },
                )
                edge_list.append(edge)
            else:
                raise ValueError(str(edge))

        G = nx.DiGraph()
        G.add_nodes_from(node_list)
        G.add_edges_from(edge_list)

        if (learn_emb == "node2vec" or
                random_walk_n_steps_method == "simulation"):
            n2v = Node2Vec(
                G, dimensions=n2v_dimensions,
                walk_length=n2v_walk_length, num_walks=n2v_num_walks,
                weight_key="norm_weight", workers=workers, seed=seed,
                p=n2v_p, q=n2v_q, quiet=True)
        else:
            n2v = NoWalkNode2Vec(
                G, dimensions=n2v_dimensions,
                walk_length=n2v_walk_length, num_walks=n2v_num_walks,
                weight_key="norm_weight", workers=workers, seed=seed,
                p=n2v_p, q=n2v_q, quiet=True)

        n2v_walks = n2v.walks

        if learn_emb == "node2vec":
            model = n2v.fit(
                window=n2v_window, min_count=n2v_min_count, epochs=n2v_epochs,
                workers=workers)

            n2v_emb_dict = dict([
                (x, model.wv.get_vector(x))
                for x in model.wv.key_to_index.keys()
            ])

            n2v_emb_df = pd.DataFrame(
                [[k] + v.tolist() for k, v in n2v_emb_dict.items()])
        else:
            model = None
            n2v_emb_dict = None
            n2v_emb_df = None

        self.G = G
        self.n2v = n2v
        self.model = model
        self.n2v_walks = n2v_walks
        self.n2v_emb_dict = n2v_emb_dict
        self.n2v_emb_df = n2v_emb_df
        return

    def _add_n2v_emb_dim_red(self):
        emb_df = self.n2v_emb_df
        if emb_df is None:
            emb_model = None
            if self.learn_emb == "spring_layout":
                pos = nx.spring_layout(
                    self.G, weight="norm_weight", seed=self.seed)
                dim_red_df = [(k, v[0], v[1]) for k, v in sorted(pos.items())]

            elif self.learn_emb == "fastrp":
                fastrp_conf = {
                    'projection_method': 'sparse',
                    'input_matrix': 'adj',
                    'weights': [0.0, 1.0, 1.0, 1.0],
                    'normalization': True,
                    'dim': 32,
                    'alpha': 0
                }
                node_order = list(self.G.nodes())
                with warnings.catch_warnings():
                    adj_mat = nx.to_scipy_sparse_array(
                        self.G, nodelist=node_order, weight="norm_weight",
                        format="csc")

                emb_arr = fastrp_wrapper(adj_mat, fastrp_conf)

                emb_model = TSNE(n_components=2, random_state=7)
                emb_model_arr = emb_model.fit_transform(emb_arr)

                dim_red_df = []
                pos = {}
                for i in range(len(node_order)):
                    pos[node_order[i]] = emb_model_arr[i, :2]
                    dim_red_df.append(
                        [node_order[i]] + emb_model_arr[i].tolist())

            elif isinstance(self.learn_emb, dict):
                pos = {
                    k: np.array(v[:2]) for k, v in self.learn_emb.items()
                }
                for i in self.G.nodes:
                    if i not in self.learn_emb:
                        err_msg = (
                            "learn_emb dict should contain all node ids. "
                            "{} is not in the learn_emb dict."
                        ).format(i)
                        raise ValueError("{} not in learn_emb dict.".format(i))
                dim_red_df = [
                    (k, v[0], v[1]) for k, v in sorted(pos.items())]
            else:
                raise ValueError("Unknown learn_emb {}".format(self.learn_emb))
        else:
            keys = emb_df.iloc[:, 0].tolist()
            emb_arr = emb_df.iloc[:, 1:].values
            emb_model = TSNE(
                n_components=2, random_state=7,
                perplexity=min(30, emb_arr.shape[0] - 1))
            emb_model_arr = emb_model.fit_transform(emb_arr)

            dim_red_df = []
            pos = {}
            for i in range(len(keys)):
                pos[keys[i]] = emb_model_arr[i, :2]
                dim_red_df.append([keys[i]] + emb_model_arr[i].tolist())

        self.dim_red_model = emb_model
        self.dim_red_plot_pos = pos
        self.dim_red_df = pd.DataFrame(dim_red_df)
        return

    def _add_n2v_prob_to_graph(self):
        G = self.G
        n2v_d_graph = self.n2v.d_graph
        for k, v in G.nodes.items():
            neighbor_pvec = n2v_d_graph[k]["first_travel_key"]
            neighbors = n2v_d_graph[k]["neighbors"]
            assert neighbor_pvec.shape == (len(neighbors),)
            if k in neighbors:
                prob = neighbor_pvec[neighbors.index(k)]
            else:
                prob = 0

            v["n2v_prob"] = prob

        for k, v in G.edges.items():
            neighbor_pvec = n2v_d_graph[k[0]]["first_travel_key"]
            neighbors = n2v_d_graph[k[0]]["neighbors"]
            assert neighbor_pvec.shape == (len(neighbors),)
            assert k[1] in neighbors
            prob = neighbor_pvec[neighbors.index(k[1])]
            v["n2v_transition_prob"] = prob
        return

    def save_n2v_res(self, pfx):
        if self.model is not None:
            self.model.save(pfx + ".model")
            self.model.wv.save_word2vec_format(pfx + ".emb")
            self.n2v_emb_df.to_csv(pfx + "_emb.csv")
            self.dim_red_df.to_csv(pfx + "_dim_red.csv")

        if self.n2v_walks is not None:
            pd.DataFrame(self.n2v_walks, copy=False).to_csv(pfx + "_walks.csv")

        self.edge_attr_df.to_csv(pfx + "_edge_attr.csv")
        nx.write_gml(self.G, pfx + "_graph.gml")
        return

    def read_n2v_res(self, pfx):
        self.G = nx.read_gml(pfx + "_graph.gml")

        try:
            self.n2v_walk_df = pd.read_csv(
                pfx + "_walks.csv", low_memory=False, index_col=0)
            self.n2v_walks = self.n2v_walk_df.values.tolist()
        except FileNotFoundError:
            self.n2v_walks = None

        try:
            self.dim_red_df = pd.read_csv(pfx + "_dim_red.csv", index_col=0)
            self.dim_red_plot_pos = dict(zip(
                self.dim_red_df.iloc[:, 0].values.tolist(),
                self.dim_red_df.iloc[:, 1:3].values.tolist()
            ))
        except FileNotFoundError:
            self.dim_red_df = None
            self.dim_red_plot_pos = None
        return rd

    def _rxn_df_to_pw_rxn_df(self, df, rxn_pw_dict, change=0):
        rxns = []
        pw_rxns = []
        for k, v in rxn_pw_dict.items():
            rxns.append(k)
            pw_rxns.append(v + "_____" + k)

        pw_rxn_df = pd.DataFrame(
            {
                "rxn": rxns,
                "pathway_rxn": pw_rxns,
            }
        ).sort_values("pathway_rxn")

        if change == 0:
            pg_df = df.loc[pw_rxn_df["rxn"].tolist(), pw_rxn_df["rxn"].tolist()]
            pg_df.index = pw_rxn_df["pathway_rxn"].tolist()
            pg_df.columns = pw_rxn_df["pathway_rxn"].tolist()
        elif change == 1:
            pg_df = df.loc[:, pw_rxn_df["rxn"].tolist()]
            pg_df.columns = pw_rxn_df["pathway_rxn"].tolist()
        else:
            raise ValueError(str(change))

        return pg_df


    def _calculate_rw_ns_mat(self, tpm, walk_length=20, num_walks=5000):
        p = np.diag([1] * tpm.shape[0])
        res = p * num_walks
        for i in range(walk_length - 1):
            p = p @ tpm
            res = res + (p * num_walks)

        res = np.round(res.T)

        diff = num_walks * walk_length - res.sum(axis=0)
        max_ind = np.argmax(res, axis=0)
        for i in range(len(diff)):
            res[max_ind[i], i] += diff[i]

        return res

    def _get_n2v_mat_dict(self):
        G = self.G

        rxns = list(G.nodes)
        assert len(rxns) == len(set(rxns))

        edge_weight_df = pd.DataFrame(np.float64(0), index=rxns, columns=rxns)
        for k, v in G.edges.items():
            edge_weight_df.loc[k[0], k[1]] = v["mn_weight"]

        rxn_pw_dict = dict([(k, v["pathway"]) for k, v in G.nodes.items()])

        pw_rxn_dict = {}
        for k, v in G.nodes.items():
            vp = v["pathway"]
            if vp not in pw_rxn_dict:
                pw_rxn_dict[vp] = []
            pw_rxn_dict[vp].append(k)

        uniq_pathways = sorted(set(list(rxn_pw_dict.values())))

        if self.random_walk_n_steps_method == "simulation":
            random_walks = sorted(self.n2v_walks, key=lambda x: x[0])

            pathway_index_dict = dict(
                zip(uniq_pathways, range(len(uniq_pathways))))
            rxn_index_dict = dict(zip(rxns, range(len(rxns))))

            assert (len(random_walks) / self.n2v_num_walks) == len(rxns)
            rw_rxn_rxn_ns_df = np.zeros((len(rxns), len(rxns)))
            for i in range(int(len(random_walks) / self.n2v_num_walks)):
                i_start = i * self.n2v_num_walks
                i_end = (i+1) * self.n2v_num_walks

                i_rw = random_walks[i_start:i_end]
                i_c = Counter()
                i_rw_start_rxn = i_rw[0][0]
                i_rw_start_rxn_index = rxn_index_dict[i_rw_start_rxn]

                for j in i_rw:
                    assert j[0] == i_rw_start_rxn
                    i_c.update(j)

                for k, v in i_c.items():
                    k_rxn_index = rxn_index_dict[k]
                    rw_rxn_rxn_ns_df[k_rxn_index, i_rw_start_rxn_index] += v

            assert (
                tuple([pathway_index_dict[x] for x in uniq_pathways]) ==
                    tuple(range(len(pathway_index_dict)))
            )
            assert (
                tuple([rxn_index_dict[x] for x in rxns]) ==
                    tuple(range(len(rxn_index_dict)))
            )

            rw_rxn_rxn_ns_df = pd.DataFrame(
                rw_rxn_rxn_ns_df, index=rxns, columns=rxns)

        elif self.random_walk_n_steps_method == "calculation":
            rw_ns_mat = self._calculate_rw_ns_mat(
                nx.to_scipy_sparse_array(
                    G, nodelist=rxns, weight="n2v_transition_prob",
                    format="csc").toarray(),
                self.n2v_walk_length, self.n2v_num_walks)

            rw_rxn_rxn_ns_df = pd.DataFrame(rw_ns_mat, index=rxns, columns=rxns)

        else:
            raise ValueError(
                "Unknown random_walk_n_steps_method {}.".format(
                    self.random_walk_n_steps_method))

        rw_rxn_rxn_ps_df = (
            rw_rxn_rxn_ns_df / (self.n2v_num_walks * self.n2v_walk_length) * 100
        )

        rw_pw_n_df = pd.DataFrame(
            np.float64(0), index=uniq_pathways, columns=rxns)
        for i in uniq_pathways:
            i_sum = rw_rxn_rxn_ns_df.loc[pw_rxn_dict[i], :].sum(axis=0)
            rw_pw_n_df.loc[i, :] = i_sum

        rw_pw_n_pr_df = self._rxn_df_to_pw_rxn_df(
            rw_pw_n_df, rxn_pw_dict, change=1)

        rw_pw_pw_ns_df = pd.DataFrame(
            np.float64(0), index=uniq_pathways, columns=uniq_pathways)
        for i in uniq_pathways:
            i_pr = [x for x in rw_pw_n_pr_df.columns
                    if x.split("_____")[0] == i]
            i_pw_pr_sum = (
                rw_pw_n_pr_df.loc[:, i_pr].values /
                    (self.n2v_num_walks * self.n2v_walk_length) * 100
            ).mean(axis=1)
            rw_pw_pw_ns_df[i] = i_pw_pr_sum

        self.edge_weight_df = edge_weight_df
        self.pr_edge_weight_df = self._rxn_df_to_pw_rxn_df(
            edge_weight_df, rxn_pw_dict)
        self.rw_pw_n_df = rw_pw_n_df
        self.rw_pw_n_pr_df = rw_pw_n_pr_df
        self.rw_pw_pw_ns_df = rw_pw_pw_ns_df
        self.rw_rxn_rxn_ns_df = rw_rxn_rxn_ns_df
        self.rw_rxn_rxn_ps_df = rw_rxn_rxn_ps_df
        return

    def plot_n2v_mat_dict(self, pfx):
        """
        Plot matrices
        """
        for k_str, k in [("pr_edge_weight_df", self.pr_edge_weight_df)]:
            uniq_pws = []
            for i in k.index.tolist():
                i_pw = i.split("_____")[0]
                if i_pw not in uniq_pws:
                    uniq_pws.append(i_pw)

            cols = sns.mpl_palette("Set1", 9)
            pw_col_dict = {}
            for i in range(len(uniq_pws)):
                pw_col_dict[uniq_pws[i]] = cols[i % 9]

            pw_cols = [
                pw_col_dict[x.split("_____")[0]] for x in k.index.tolist()]

            chm = sns.clustermap(
                np.log10(k + 1),
                row_cluster=False, col_cluster=False,
                row_colors=pw_cols, col_colors=pw_cols,
                cmap="viridis")
            chm.figure.set_size_inches(15, 15)
            chm.savefig(pfx + "_" + k_str + "_heatmap.png", dpi=300)
            plt.close()

        for k_str, k in [("rw_pw_pw_ns_df", self.rw_pw_pw_ns_df)]:
            chm = sns.clustermap(
                np.log10(k + 1),
                row_cluster=False, col_cluster=False,
                cmap="viridis")
            chm.figure.set_size_inches(15, 15)
            chm.savefig(pfx + "_" + k_str + "_heatmap.png", dpi=300)
            plt.close()

        return

    def _get_pw_graph(self, rxn_g, rxn_pos_dict, pw_pw_ns_df, pw_rxn_ns_df):
        pw_node_dict = {}
        for k, v in rxn_g.nodes.items():
            pw = v["pathway"]
            if pw not in pw_node_dict:
                pw_node_dict[pw] = {
                    "rxn": [],
                    "property": [],
                    "n2v_prob": [],
                    "rxn_pos": [],
                    "pathway": pw
                }

            pw_node_dict[pw]["rxn"].append(k)
            pw_node_dict[pw]["property"].append(v["property"])
            pw_node_dict[pw]["n2v_prob"].append(v["n2v_prob"])
            pw_node_dict[pw]["rxn_pos"].append(rxn_pos_dict[k])

        for v in pw_node_dict.values():
            v["rxn_pos"] = np.array(v["rxn_pos"])
            v["pw_pos"] = np.mean(v["rxn_pos"], axis=0)
            v["pw_property"] = np.mean(v["property"])
            v["pw_n2v_prob"] = np.mean(v["n2v_prob"])

        pw_edge_dict = {}
        for src, src_v in rxn_g.nodes.items():
            pw_p_dict = {}
            for dest in nx.neighbors(rxn_g, src):
                dest_v = rxn_g.nodes[dest]
                pw = (src_v["pathway"], dest_v["pathway"])
                if pw not in pw_p_dict:
                    pw_p_dict[pw] = 0

                pw_p_dict[pw] += rxn_g.edges[(src, dest)]["n2v_transition_prob"]

                if pw not in pw_edge_dict:
                    pw_edge_dict[pw] = {
                        "rxn": [],
                        "norm_weight": [],
                        "n2v_transition_src": [],
                        "n2v_transition_prob": [],
                    }

                pw_edge_dict[pw]["rxn"].append((src, dest))
                pw_edge_dict[pw]["norm_weight"].append(
                    rxn_g.edges[(src, dest)]["norm_weight"])

            for k, v in pw_p_dict.items():
                pw_edge_dict[k]["n2v_transition_src"].append(src)
                pw_edge_dict[k]["n2v_transition_prob"].append(v)

        for k, v in pw_edge_dict.items():
            assert len(v["n2v_transition_prob"]) != 0
            v["pw_norm_weight"] = np.mean(v["norm_weight"])
            v["pw_n2v_transition_prob"] = np.mean(v["n2v_transition_prob"])
            v["pw_property"] = pw_node_dict[k[1]]["pw_property"]
            v["property"] = pw_node_dict[k[1]]["property"]
            v["pw_rw_n"] = pw_pw_ns_df.loc[k[1], k[0]]
            v["rw_n"] = [
                (
                    pw_rxn_ns_df.loc[k[1], x]
                    / (self.n2v_num_walks * self.n2v_walk_length) * 100
                )
                for x in v["n2v_transition_src"]
            ]

        # Completes self-loops.
        for k, v in pw_node_dict.items():
            if (k, k) not in pw_edge_dict:
                pw_edge_dict[(k, k)] = {
                    'rxn': [],
                    'norm_weight': [],
                    'n2v_transition_src': [],
                    'n2v_transition_prob': [],
                    'pw_norm_weight': 0,
                    'pw_n2v_transition_prob': 0,
                    'pw_property': np.mean(v['property']),
                    "property": v['property'],
                    'pw_rw_n': pw_pw_ns_df.loc[k, k],
                    'rw_n': [
                        (
                            pw_rxn_ns_df.loc[k, x]
                            / (self.n2v_num_walks * self.n2v_walk_length) * 100
                        )
                        for x in v['rxn']
                    ]
                }
        G = nx.DiGraph()
        G.add_nodes_from([(k, v) for k, v in pw_node_dict.items()])
        G.add_edges_from([(k[0], k[1], v) for k, v in pw_edge_dict.items()])
        return G

    def _rxn_to_pw_graph(self):
        self.pw_graph = self._get_pw_graph(
            self.G, self.dim_red_plot_pos, self.rw_pw_pw_ns_df, self.rw_pw_n_df)

    def plot_graph(self, filename, s_pw=None, node_color_attr="property",
                   draw_edges=False, draw_labels=False, non_s_pw_node_size=10,
                   s_pw_node_size=150):
        G = self.G
        pos = self.dim_red_plot_pos
        if s_pw is None:
            s_pw = []

        node_list = [(k, v) for k, v in G.nodes.items()]
        s_pw_nodes = [k for k, v in G.nodes.items() if v["pathway"] in s_pw]
        other_nodes = [k for k, v in G.nodes.items() if k not in s_pw_nodes]

        nodes = other_nodes + s_pw_nodes
        node_colors = [G.nodes[x][node_color_attr] for x in nodes]
        node_sizes =  (
            [non_s_pw_node_size] * len(other_nodes) +
                [s_pw_node_size] * len(s_pw_nodes)
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            npc = nx.draw_networkx_nodes(
                G, pos, nodelist=nodes, node_color=node_colors,
                node_size=node_sizes, alpha=0.5, cmap="viridis")

        if draw_edges:
            edge_list = [x for x in list(G.edges) if x[0] != x[1]]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                nx.draw_networkx_edges(
                    G,
                    pos,
                    edgelist=edge_list,
                    alpha=0.01, arrowstyle="-",
                    nodelist=nodes, node_size=node_sizes
                )

        if draw_labels:
            label_dict = {}
            for n in G:
                if n in s_pw:
                    label_dict[n] = n

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                nx.draw_networkx_labels(
                    G, pos, label_dict, font_size=5, font_color="black")

        plt.title(
            "\n".join(
                textwrap.wrap(
                    filename[(filename.rfind("/") + 1):].replace("_", " "),
                    80, break_long_words=False, break_on_hyphens=False)
            ),
            fontdict={"fontsize": 6}
        )

        plt.tick_params(
            left=True, bottom=True, labelleft=True, labelbottom=True,
            labelsize=5)
        plt.colorbar(npc)
        plt.tight_layout()
        plt.savefig(filename, dpi=600, format="png")
        plt.close()
        return
