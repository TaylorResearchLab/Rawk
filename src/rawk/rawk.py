from .rawk_sample import RawkSample

import networkx as nx
import pandas as pd
import numpy as np
import warnings
from adjustText import adjust_text
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import textwrap
import seaborn as sns
import scipy.stats as sps



class Rawk:
    """
    Run Rawk on one sample

    Parameters
    ----------
    sample : RawkSample
        An RawkSample instance to test for enrichment.
    bg_sample: RawkSample
        An RawkSample instance used as background.
    """

    def __init__(self, sample, bg_sample):
        if sample.n2v_num_walks != bg_sample.n2v_num_walks:
            raise ValueError("Not implemented.")

        if sample.n2v_walk_length != bg_sample.n2v_walk_length:
            raise ValueError("Not implemented.")

        self.sample = sample
        self.bg_sample = bg_sample

        self.n2v_walk_length = sample.n2v_walk_length
        self.n2v_num_walks = sample.n2v_num_walks
        self.n2v_spr = self.n2v_num_walks * self.n2v_walk_length

        self._contrast_pw_graph()

    def _diff_x_y(self, x, xlab, y, ylab):
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
        Create sample - background pathway graph
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
            v['pw_rw_n_diff'] = x_ev['pw_rw_n'] - y_ev['pw_rw_n']

            v['n2v_transition_prob_diff'] = self._diff_x_y(
                x_ev['n2v_transition_prob'], x_ev['n2v_transition_src'],
                y_ev['n2v_transition_prob'], y_ev['n2v_transition_src']
            )
            v['rw_n_diff'] = self._diff_x_y(
                x_ev['rw_n'], x_ev['n2v_transition_src'],
                y_ev['rw_n'], y_ev['n2v_transition_src']
            )

        self.contrast_pw_graph = contrast_pw_graph

    def _expand_v0range(self, vmin, vmax):
        """
        If vmin == vmax, make them different
        """
        if vmin == vmax:
            if vmin == 0:
                vmin = -0.1
                vmax = 0.1
            else:
                vmin = vmin - 0.1 * abs(vmin)
                vmax = vmax + 0.1 * abs(vmax)
        return vmin, vmax

    def _center_cmap(self, cmap, vmin, vmax, center, return_norm=False):
        """
        Center color map on a value
        """
        # Adapted from seaborn._HeatMapper._determine_cmap_params
        bad = cmap(np.ma.masked_invalid([np.nan]))[0]

        under = cmap(-np.inf)
        over = cmap(np.inf)
        under_set = under != cmap(0)
        over_set = over != cmap(cmap.N - 1)

        vrange = max(vmax - center, center - vmin)
        normlize = mpl.colors.Normalize(center - vrange, center + vrange)
        cmin, cmax = normlize([vmin, vmax])
        cc = np.linspace(cmin, cmax, 256)
        cmap = mpl.colors.ListedColormap(cmap(cc))
        cmap.set_bad(bad)
        if under_set:
            cmap.set_under(under)
        if over_set:
            cmap.set_over(over)
        if return_norm:
            return cmap, normlize
        else:
            return cmap

    def plot_pw_neighborhood(self, s_pw, filename, node_color_attr,
                             node_color_vmin=None, node_color_vmax=None,
                             node_color_center=None,
                             n_cutoff=10,
                             node_alpha=0.5):
        """
        Plot pathway reachable local neighborhood

        Parameters
        ----------
        s_pw : str
            A pathway to plot.
        filename : str
            Saved plot filename
        node_color_attr : str
            Node addtribute used for coloring:
            - 'property': node property
            - 'rw_ps': random walk percent steps
            - 'log10p1_rw_ps': log10(rw_ps + 1)
            - 'rw_ps_diff': sample - background rw_ps
            - 'rw_ps_log2fc': log2((sample_rw_ps + 1) / (background_rw_ps + 1))
        node_color_vmin : float, optional
            Minimum value of the node color map.
        node_color_vmax : float, optional
            Maximum value of the node color map.
        node_color_center : float, optional
            Center value of the node color map.
        n_cutoff : int
            The number of top reached pathways to plot.
        node_alpha : float
            The alpha value of plotted nodes


        Returns
        -------
        None
        """
        G = self.contrast_pw_graph
        pos = dict([(k, v["pw_pos"]) for k, v in G.nodes.items()])

        # {pw: [rxn1, rxn2, ...], ...}
        pw_rxns_dict = {}
        for i in range(self.sample.node_df.shape[0]):
            i_rxn = self.sample.node_df.iloc[i, :].loc["rxn"]
            i_pw = self.sample.node_df.iloc[i, :].loc["pathway"]

            if i_pw not in pw_rxns_dict:
                pw_rxns_dict[i_pw] = []

            pw_rxns_dict[i_pw].append(i_rxn)

        rw_pw_pw_df = self.sample.rw_pw_pw_ns_df
        rw_rxn_rxn_df = self.sample.rw_rxn_rxn_ps_df
        bg_rw_rxn_rxn_df = self.bg_sample.rw_rxn_rxn_ps_df
        rxn_g = self.sample.G

        # Only implemented one source pw plot
        assert isinstance(s_pw, str)
        s_pw = [s_pw]
        s_pw_rxns_dict = {}
        for k, v in pw_rxns_dict.items():
            s_pw_rxns_dict[k] = [x for x in v if x in rw_rxn_rxn_df.columns]

        s_pw_rxns = [s_pw_rxns_dict[x] for x in s_pw]

        def get_topn_a0_inds(s, n):
            r = [k for k, v in s.sort_values(ascending=False)[:n].items()
                if v > 0]
            return r

        rw_reached_nodes = get_topn_a0_inds(rw_pw_pw_df[s_pw[0]], n_cutoff)

        # The first step should always be in the pw.
        assert s_pw[0] in rw_reached_nodes
        edges = []
        for i in rw_reached_nodes:
            for j in rw_reached_nodes:
                if (i, j) in G.edges:
                    if j in get_topn_a0_inds(rw_pw_pw_df.loc[:, i], n_cutoff):
                        edges.append((i, j))

        # Keep only nodes with one or more edges that are not self-loops.
        nodes = set()
        for i in edges:
            if i[0] != i[1]:
                nodes.add(i[0])
                nodes.add(i[1])

        nodes = sorted(nodes)
        edges = sorted([x for x in edges
                        if (x[0] in nodes) and (x[1] in nodes) and (x[0] != x[1])])
        if len(nodes) == 0:
            nodes = [s_pw[0]]
        assert len(nodes) == len(set(nodes))
        assert len(edges) == len(set(edges))

        # Prepare colors
        node_colors = []
        node_color_dict = {}
        rw_pw_rxn_s = rw_rxn_rxn_df.loc[:, s_pw_rxns[0]].mean(axis=1)
        bg_rw_pw_rxn_s = bg_rw_rxn_rxn_df.loc[
            rw_rxn_rxn_df.index.tolist(), s_pw_rxns[0]].mean(axis=1)
        for i in nodes:
            if node_color_attr == "property":
                i_col = G.nodes[i][node_color_attr]
            else:
                if node_color_attr == "rw_ps":
                    # rw_pw_rxn_s is pw->rxn percent steps
                    nc_s = rw_pw_rxn_s
                elif node_color_attr == "log10p1_rw_ps":
                    # log10 steps per million + 1
                    nc_s = np.log10(rw_pw_rxn_s + 1)
                elif node_color_attr == "rw_ps_diff":
                    nc_s = rw_pw_rxn_s - bg_rw_pw_rxn_s
                elif node_color_attr == "rw_ps_log2fc":
                    nc_s = np.log2((rw_pw_rxn_s + 1) / (bg_rw_pw_rxn_s + 1))
                else:
                    raise ValueError(
                        "Unknown node_color_attr: {}".format(node_color_attr))
                i_col = []
                for j in s_pw_rxns_dict[i]:
                    if rw_pw_rxn_s[j] > 0:
                        i_col.append(nc_s[j])
            i_col = sorted(i_col)
            if len(i_col) == 0:
                i_col = [0]
            node_colors.append(i_col)
            node_color_dict[i] = i_col

        if node_color_vmin is None:
            node_color_vmin = np.min([np.min(x) for x in node_colors])
        if node_color_vmax is None:
            node_color_vmax = np.max([np.max(x) for x in node_colors])

        node_color_vmin, node_color_vmax = self._expand_v0range(
            node_color_vmin, node_color_vmax)

        if node_color_center is None:
            node_color_cmap = self._center_cmap(
                sns.color_palette('viridis', as_cmap=True),
                node_color_vmin, node_color_vmax,
                (node_color_vmin + node_color_vmax) / 2)
        else:
            node_color_cmap = self._center_cmap(
                sns.color_palette("blend:blue,white,red", as_cmap=True),
                node_color_vmin, node_color_vmax, node_color_center)

        node_color_cmap_alpha = node_color_cmap(np.arange(node_color_cmap.N))
        node_color_cmap_alpha[:, -1] = node_alpha
        node_color_cmap_alpha = mpl.colors.ListedColormap(node_color_cmap_alpha)

        node_sizes = [800 if x in s_pw else 380 for x in nodes]

        # Adjust node positions to avoid overlapping
        adj_pos = {x: tuple(pos[x]) for x in nodes}
        pos_x_range = (
            max([x[0] for x in adj_pos.values()]) -
                min([x[0] for x in adj_pos.values()])
        )
        pos_y_range = (
            max([x[1] for x in adj_pos.values()]) -
                min([x[1] for x in  adj_pos.values()])
        )
        overlap_xd = pos_x_range / 7
        overlap_yd = pos_y_range / 7
        def is_overlap(p1, p2, xd, yd):
            if abs(p1[0] - p2[0]) < overlap_xd and abs(p1[1] - p2[1]) < overlap_yd:
                return True
            else:
                return False

        for k in range(3):
            for i in nodes:
                for j in nodes:
                    if i == j:
                        continue
                    if is_overlap(adj_pos[i], adj_pos[j], overlap_xd, overlap_yd):
                        # Move j.
                        if adj_pos[j][0] > adj_pos[i][0]:
                            # j on the right of i
                            ajx = adj_pos[i][0] + overlap_xd
                        else:
                            ajx = adj_pos[i][0] - overlap_xd
                        if adj_pos[j][1] > adj_pos[i][1]:
                            # j on the up of i
                            ajy = adj_pos[i][1] + overlap_yd
                        else:
                            ajy = adj_pos[i][1] - overlap_yd
                        adj_pos[j] = (ajx, ajy)

        pos = adj_pos
        edge_widths = []
        for i in edges:
            i_rank = rw_pw_pw_df.loc[:, i[0]].sort_values(
                ascending=False).index.get_indexer([i[1]])[0]
            assert i_rank >= 0
            edge_widths.append(0.1 + 0.05 * max(0, 11 - i_rank))

        label_dict = {}
        for n in G:
            if n in nodes:
                n_pw_rxns = s_pw_rxns_dict[n]
                label_dict[n] = "{}".format(n)

        node_color_norm = mpl.colors.Normalize(node_color_vmin, node_color_vmax)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            npc = nx.draw_networkx_nodes(
                G, pos, nodelist=nodes,
                node_color=[np.max(x) for x in node_colors],
                node_size=0, alpha=node_alpha, cmap=node_color_cmap,
                vmin=node_color_vmin, vmax=node_color_vmax, edgecolors="k")
            npc = nx.draw_networkx_nodes(
                G, pos, nodelist=nodes,
                node_color=[np.min(x) for x in node_colors],
                node_size=0, alpha=node_alpha, cmap=node_color_cmap,
                vmin=node_color_vmin, vmax=node_color_vmax, edgecolors="k")

            ax = plt.gca()
            xlim = plt.xlim()
            xlim = (xlim[0] - 0.08 * (xlim[1] - xlim[0]),
                    xlim[1] + 0.08 * (xlim[1] - xlim[0]))
            xw = xlim[1] - xlim[0]
            ylim = plt.ylim()
            ylim = (ylim[0] - 0.08 * (ylim[1] - ylim[0]),
                    ylim[1] + 0.08 * (ylim[1] - ylim[0]))
            yw = ylim[1] - ylim[0]
            ax.set_aspect("equal")
            assert xw != 0
            assert yw != 0
            pie_pos = {}
            if xw > yw:
                o_ylim = ylim
                ylim = (np.mean(ylim) - xw / 2, np.mean(ylim) + xw / 2)
                for k in nodes:
                    v = pos[k]
                    pie_pos[k] = (v[0], np.interp(v[1], o_ylim, ylim))
            else:
                o_xlim = xlim
                xlim = (np.mean(xlim) - yw / 2, np.mean(xlim) + yw / 2)
                for k in nodes:
                    v = pos[k]
                    pie_pos[k] = (np.interp(v[0], o_xlim, xlim), v[1])

            ax.set(xlim=xlim, ylim=ylim)
            nx.draw_networkx_edges(
                G,
                pie_pos,
                edgelist=edges,
                width=edge_widths,
                alpha=0.7,
                arrowsize=5,
                nodelist=nodes, node_size=node_sizes,
                edge_color="gray"
            )
            nx_labels = nx.draw_networkx_labels(
                G, pie_pos, label_dict, font_size=5, font_color="black",
                horizontalalignment="left")
            unit_pie_size = (
                (2 * (ax.transData.transform([1, 0.95]) -
                            ax.transData.transform([0, 0])) *
                                (72 / ax.figure.dpi))[1] ** 2
            )
            pie_sizes = [(x / unit_pie_size) ** 0.5 for x in node_sizes]
            for i in range(len(nodes)):
                pie_r = ax.pie(
                    [1] * len(node_colors[i]),
                    colors=[node_color_cmap_alpha(node_color_norm(x))
                            for x in node_colors[i]],
                    center=pie_pos[nodes[i]],
                    radius=pie_sizes[i],
                    startangle=90, frame=True)
                ax.add_patch(
                    mpl.patches.Circle(
                        pie_r[0][0].center, pie_r[0][0].r,
                        fill=False, edgecolor="gray", linewidth=1)
                )

        adjust_text(list(nx_labels.values()))

        plt.tick_params(
            left=True, bottom=True, labelleft=True, labelbottom=True,
            labelsize=5)
        plt.title(
            "\n".join(
                textwrap.wrap(
                    filename[(filename.rfind("/") + 1):].replace("_", " "),
                    60, break_long_words=False, break_on_hyphens=False)
            ),
            fontdict={"fontsize": 6}
        )
        plt.colorbar(npc)
        fig = plt.gcf()
        fig.set_size_inches(5, 3.9)
        plt.tight_layout()
        plt.savefig(filename, dpi=300, format="png")
        plt.close()
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
        A dictionary of testing results.
        """
        n2v_spr = self.n2v_spr
        assert np.all(self.sample.rw_pw_n_pr_df.sum(axis=0) == n2v_spr)
        assert (
            self.sample.rw_pw_n_pr_df.index.tolist() ==
                self.bg_sample.rw_pw_n_pr_df.index.tolist()
        )
        assert (
            self.sample.rw_pw_n_pr_df.columns.tolist() ==
                self.bg_sample.rw_pw_n_pr_df.columns.tolist()
        )

        def get_pw_ns(xdf, ydf, pw, total_ns):
            pw_rxns = [c for c in xdf.columns.tolist()
                       if c.split("_____")[0] == pw]
            x = xdf.loc[pw, pw_rxns].values
            y = ydf.loc[pw, pw_rxns].values
            # Exclude rxns that do not go outside of the pathway
            f = np.logical_not(np.logical_and(x == y, x == total_ns))
            return x[f], y[f]

        def get_a_test_ns_tbls(x, y, norm_fac, total_ns):
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

        mnea_res_dict = {}
        sid = self.sample.name
        mnea_res_dict[sid] = {}
        x_df = self.sample.rw_pw_n_pr_df
        y_df = self.bg_sample.rw_pw_n_pr_df
        for i in x_df.index.tolist():
            mnea_res_dict[sid][i] = {}
            tx, ty = get_pw_ns(x_df, y_df, i, n2v_spr)
            assert tx.shape == ty.shape
            if len(tx) == 0:
                mnea_res_dict[sid][i]["test_pval"] = np.nan
                mnea_res_dict[sid][i]["test_nes"] = np.nan
            else:
                tx = tx - self.n2v_num_walks
                ty = ty - self.n2v_num_walks
                mnea_res_dict[sid][i]["tx"] = tx
                mnea_res_dict[sid][i]["ty"] = ty

                cmp_ns_tbls = get_a_test_ns_tbls(tx, ty, cmp_norm_fac, n2v_spr)
                test_res = sps.fisher_exact(
                    np.round(cmp_ns_tbls[1]) + 1, alternative="two-sided")

                mnea_res_dict[sid][i]["test_nes"] = np.log2(test_res.statistic)
                mnea_res_dict[sid][i]["test_pval"] = test_res.pvalue
        return mnea_res_dict

    def test_res_dict_to_df(self, mnea_res_dict, m_sfx):
        """
        Convert a test result dictionary to a dataframe

        Parameters
        ----------
        mnea_res_dict : dict
            Rawk test result dictionary.

        m_sfx : str
            Method suffix for dataframe columns.

        Returns
        -------
        A 2-tuple of (fdr_dataframe, enrichment_score_dataframe).
        """
        # padj df
        mnea_padj_df_list = []
        for k, v in mnea_res_dict.items():
            kcol = k + m_sfx
            cmp_df = pd.DataFrame(
                [(pw, r["test_pval"]) for pw, r in v.items()],
                columns=["pathway", kcol]).dropna()
            cmp_df[kcol] = sps.false_discovery_control(cmp_df[kcol].values)
            mnea_padj_df_list.append(cmp_df)

        mnea_padj_df = mnea_padj_df_list[0]
        for i in mnea_padj_df_list[1:]:
            mnea_padj_df = mnea_padj_df.join(
                i.set_index("pathway"), on = "pathway", how="outer")

        # nes df
        mnea_nes_df_list = []
        for k, v in mnea_res_dict.items():
            cmp_df = pd.DataFrame(
                [(pw, r["test_nes"]) for pw, r in v.items()],
                columns=["pathway", k + m_sfx]).dropna()
            mnea_nes_df_list.append(cmp_df)

        mnea_nes_df = mnea_nes_df_list[0]
        for i in mnea_nes_df_list[1:]:
            mnea_nes_df = mnea_nes_df.join(
                i.set_index("pathway"), on="pathway", how="outer")

        return mnea_padj_df, mnea_nes_df

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
        A dictionary of testing results.
        """
        n2v_spr = self.n2v_spr
        mnea_res_dict = {}
        sid = self.sample.name
        mnea_res_dict[sid] = {}

        # Used to define local neighborhood.
        ln_df = self.bg_sample.rw_pw_n_pr_df
        x_df = self.bg_sample.rw_pw_n_pr_df
        for i in x_df.index.tolist():
            mnea_res_dict[sid][i] = {}
            tp_rxns = [
                xx for xx in x_df.columns.tolist()
                if xx.split("_____")[0] == i]
            cmp_rsf = (
                ln_df.loc[:, tp_rxns].sum(axis=1) / (len(tp_rxns) * n2v_spr)
            )
            l_pws = [
                cr for cr in cmp_rsf[cmp_rsf > rw_s_prop_cutoff].index.tolist()
                if cr != i
            ]

            mnea_res_dict[sid][i]["rows"] = l_pws
            if len(l_pws) == 0:
                # no local neighbor pw to compare to.
                mnea_res_dict[sid][i]["test_pval"] = np.nan
                mnea_res_dict[sid][i]["test_nes"] = np.nan
            else:
                etx = self.sample.node_df.loc[
                    self.sample.node_df["pathway"] == i, "property"].values
                elx = self.sample.node_df.loc[
                    self.sample.node_df["pathway"].isin(l_pws),
                    "property"].values

                test_res = sps.ranksums(etx, elx)

                mnea_res_dict[sid][i]["test_nes"] = np.log(
                    (etx.mean() + 1) / (elx.mean() + 1))
                mnea_res_dict[sid][i]["test_pval"] = test_res.pvalue

        return mnea_res_dict
