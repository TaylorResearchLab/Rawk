import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
import textwrap
from adjustText import adjust_text
import random



def assert_frame_equal(x, y):
    pd.testing.assert_frame_equal(
        x, y, rtol=0, atol=1e-8)



def plot_nw_stats(uf_ntbl, uf_etbl,
                  mn_weight_cutoff,
                  mn_weight_cutoff_label,
                  title_pfx, fn_pfx):

    uf_etbl_wg0 = uf_etbl.loc[
        uf_etbl['mn_weight'] > 0, :].copy()
    uf_g0_w = uf_etbl_wg0['mn_weight']
    uf_g0_w_s = uf_g0_w.sort_values(ascending=False)

    layout = '''
        A
        A
        A
        B
        C
        D
        E
    '''

    fig, ax_dict = plt.subplot_mosaic(
        layout, figsize=(6, 8.5),
        constrained_layout=True
    )

    for label, ax in ax_dict.items():
        if label in []:
            lpx = -0.05
        else:
            lpx = -0.1

        ax.text(
            lpx, 1.00, label,
            transform=ax.transAxes,
            fontsize=16, fontweight='bold',
            va='bottom', ha='right')

    # A
    ax_dict['A'].plot(
        np.arange(len(uf_g0_w_s)), uf_g0_w_s, 'b',
        label='edges with weight > 0')
    ax_dict['A'].hlines(
        mn_weight_cutoff,
        ax_dict['A'].get_xlim()[0],
        ax_dict['A'].get_xlim()[1],
        linestyles='--',
        label=mn_weight_cutoff_label)

    ax_dict['A'].set_xlabel('Edge weight ranking')
    ax_dict['A'].set_ylabel('Edge weight')
    ax_dict['A'].set_yscale('log', base=10)
    ax_dict['A'].legend(loc='best')

    # B
    n_hist_bins = 100
    all_rxns = set(
        uf_etbl['src'].tolist()
        + uf_etbl['dest'].tolist())
    n_rxns = len(all_rxns)
    mg = nx.from_pandas_edgelist(
        uf_etbl_wg0,
        source='src',
        target='dest',
        edge_attr='mn_weight')

    e_wg0_rxns = list(mg.nodes)
    assert len(set(e_wg0_rxns)) == len(e_wg0_rxns)
    e_wg0_rxns = set(e_wg0_rxns)

    no_e_wg0_rxns = all_rxns - e_wg0_rxns
    mg.add_nodes_from(list(no_e_wg0_rxns))

    uf_wg0_g_dd = dict(mg.degree())

    curr_ax = ax_dict['B']
    curr_ax.hist(
        list(uf_wg0_g_dd.values()),
        bins=n_hist_bins)
    txp = 0.25
    curr_ax.text(
        txp, 0.95,
        '{:,d} edges with weight > 0'.format(
            mg.number_of_edges()),
        transform=curr_ax.transAxes,
        va='top', ha='left')
    curr_ax.set_xlabel('Node degree')
    curr_ax.set_ylabel('N nodes')
    assert len(uf_wg0_g_dd) == mg.number_of_nodes()
    assert uf_etbl_wg0.shape[0] == mg.number_of_edges()

    # C
    gc_el = [
        (u, v) for u, v, d in mg.edges(data=True)
        if d['mn_weight'] <= mn_weight_cutoff
    ]
    mg.remove_edges_from(gc_el)
    cg_dd = dict(mg.degree())

    curr_ax = ax_dict['C']
    curr_ax.hist(
        list(cg_dd.values()), bins=n_hist_bins)
    curr_ax.text(
        txp, 0.95,
        '{:,d} edges with weight > {:.3f}'.format(
            mg.number_of_edges(), mn_weight_cutoff),
        transform=curr_ax.transAxes,
        va='top', ha='left')
    curr_ax.set_xlabel('Node degree')
    curr_ax.set_ylabel('N nodes')
    assert len(cg_dd) == mg.number_of_nodes() == n_rxns

    # D
    uf_pw_rc = (
        uf_ntbl
        .groupby('pathway')['rxn']
        .apply(lambda x: len(set(x)))
        .to_dict()
    )
    uf_pw_gc = (
        uf_ntbl
        .groupby('pathway')['gene']
        .apply(lambda x: len(set(x)))
        .to_dict()
    )
    curr_ax = ax_dict['D']
    curr_ax.hist(
        list(uf_pw_rc.values()),
        bins=n_hist_bins)
    curr_ax.text(
        txp, 0.95,
        '{:,d} pathways\n{:,d} reactions'.format(
            len(uf_pw_rc),
            uf_ntbl['rxn'].nunique()),
        transform=curr_ax.transAxes,
        va='top', ha='left')
    curr_ax.set_xlabel('N reactions per pathway')
    curr_ax.set_ylabel('N pathways')

    # E
    curr_ax = ax_dict['E']
    curr_ax.hist(
        list(uf_pw_gc.values()),
        bins=n_hist_bins)
    curr_ax.text(
        txp, 0.95,
        '{:,d} pathways\n{:,d} genes'.format(
            len(uf_pw_rc),
            uf_ntbl['gene'].nunique()),
        transform=curr_ax.transAxes,
        va='top', ha='left')
    curr_ax.set_xlabel('N genes per pathway')
    curr_ax.set_ylabel('N pathways')

    #
    fig.suptitle(title_pfx)
    plt.savefig(
        fn_pfx + 'mrn_stats.png',
        dpi=300, bbox_inches='tight',
        pad_inches=0.1)
    plt.close()



def plot_elbow(metric_by_n, selected_n, save_pfx, xlab, ylab,
               annot):
    plt.figure(figsize=(5, 2))
    plt.plot(
        sorted(metric_by_n.keys()),
        [metric_by_n[x] for x in sorted(metric_by_n.keys())],
        marker='o', markersize=5,
        linestyle='--', color='black')
    plt.axvline(x=selected_n, color='blue', linestyle=':')
    plt.text(
        x=selected_n + 0.3,
        y=metric_by_n[sorted(metric_by_n.keys())[2]],
        s=annot,
        ha='left', va='center')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.savefig(
        f'{save_pfx}_elbow.png',
        dpi=300, bbox_inches='tight',
        pad_inches=0.1)
    plt.close()



# From R package ggsci::pal_jco
ggsci_pal_jco = [
    "#0073C2", "#EFC000", "#868686",
    "#CD534C", "#7AA6DC", "#003C67",
    "#8F7700", "#3B3B3B", "#A73030",
    "#4A6990"
]



def plot_graph(G, pos, filename, s_pw=None, node_color_attr="property",
               draw_edges=False, draw_labels=False, non_s_pw_node_size=10,
               s_pw_node_size=150, title=None,
               figsize=(3, 3), seed=42, dim_name="Dim",
               discrete_color=False, cmap=None, node_alpha=1,
               edge_width=None, dpi=300, discrete_col_legend_ncol=3,
               edge_alpha=None):
    """
    Plot networkx graph

    Parameters
    ----------
    cmap : str, dict, optional
        If discrete_color, cmap is a dictionary to map
        node_color_attr values to colors. If not discrete_color,
        cmap is a string. If None, discrete_color cmap is from
        R package ggsci::pal_jco, and not discrete_color cmap is
        'viridis'.

    Returns
    -------
    None
    """
    if s_pw is None:
        s_pw = []

    node_list = [
        (k, v) for k, v in G.nodes.items()]
    s_pw_nodes = [
        k for k, v in G.nodes.items() if v["pathway"] in s_pw]
    other_nodes = [
        k for k, v in G.nodes.items() if k not in s_pw_nodes]

    nodes = other_nodes + s_pw_nodes
    node_colors = [G.nodes[x][node_color_attr] for x in nodes]
    node_sizes =  (
        [non_s_pw_node_size] * len(other_nodes) +
            [s_pw_node_size] * len(s_pw_nodes)
    )

    if cmap is None:
        if discrete_color:
            uniq_col_vals = sorted(set(node_colors))
            cmap = {
                x: ggsci_pal_jco[i % len(ggsci_pal_jco)]
                for i, x in enumerate(uniq_col_vals)
            }
        else:
            cmap = "viridis"

    fig, ax = plt.subplots(figsize=figsize)

    if draw_edges:
        edge_list = [x for x in list(G.edges) if x[0] != x[1]]
        if edge_alpha is None:
            if len(edge_list) < 5000:
                edge_alpha = 0.05
            elif len(edge_list) < 100_000:
                edge_alpha = 0.01
            else:
                edge_alpha = 0.002
        else:
            edge_alpha = edge_alpha

        if edge_width is None:
            edge_width = 1

        nx.draw_networkx_edges(
            G,
            pos, ax=ax,
            edgelist=edge_list,
            alpha=edge_alpha, arrows=False,
            nodelist=nodes, node_size=node_sizes,
            width=edge_width
        )

    if draw_labels:
        label_dict = {}
        for n in G:
            if n in s_pw:
                label_dict[n] = n

        nx.draw_networkx_labels(
            G, pos, label_dict, font_size=5, font_color="black")

    if discrete_color:
        node_colors = [cmap[x] for x in node_colors]
        nc = nx.draw_networkx_nodes(
            G, pos, ax=ax, nodelist=nodes,
            node_color=node_colors,
            node_size=node_sizes, alpha=node_alpha,
            edgecolors="none")
        legend_handles = [
            mpl.patches.Patch(color=cmap[x], label=x)
            for x in sorted(cmap.keys())
        ]
        plt.legend(
            handles=legend_handles,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.02),
            ncol=discrete_col_legend_ncol,
            title=title, title_fontsize=12,
            alignment="left")
    else:
        nc = nx.draw_networkx_nodes(
            G, pos, ax=ax, nodelist=nodes,
            node_color=node_colors,
            node_size=node_sizes,
            alpha=node_alpha, cmap=cmap,
            edgecolors="none")
        cbar = fig.colorbar(
            nc, ax=ax, orientation="horizontal",
            location="top", pad=0.1)
        cbar.ax.xaxis.set_ticks_position("bottom")
        cbar.ax.text(
            x=-0.05,
            y=1.22,
            s=title,
            transform=ax.transAxes,
            va="bottom",
            ha="left",
        )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)

    ax.set_axisbelow(True)

    ax.set_xticks([])
    ax.set_yticks([])

    ax.set_xlabel(f'{dim_name} 1', labelpad=1)
    ax.set_ylabel(f'{dim_name} 2', labelpad=0)
    ax.margins(0.005)

    plt.savefig(
        filename,
        dpi=dpi, bbox_inches='tight',
        pad_inches=0.1)
    plt.close()



def hist(data, n_bins, title, xlab, ylab,
         filename, figsize=(1.8, 1.5)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(data, bins=n_bins)
    plt.title(title)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel(
        xlab, ha='left',
        multialignment='left')
    ax.set_ylabel(ylab)
    ax.xaxis.set_label_coords(
        x=-0.3, y=-0.2, transform=ax.transAxes)
    plt.savefig(
        filename,
        dpi=300, bbox_inches='tight',
        pad_inches=0.05)
    plt.close()



def expand_v0range(vmin, vmax):
    """
    If vmin == vmax, make them different
    """
    if vmin == vmax:
        if vmin == 0:
            vmin = -0.1
            vmax = 0.1
        else:
            vmin = vmin - 0.01 * abs(vmin)
            vmax = vmax + 0.01 * abs(vmax)
    return vmin, vmax



def get_color_mapper(vmin, vmax, vcenter=None, cmap_name=None):
    """
    Creates a color mapper function.

    Parameters:
    -----------
    vmin, vmax : float
        The minimum and maximum limits of the data.
    vcenter : float, optional
        The center point (e.g., 0). If None, standard linear
        mapping is used.
    cmap_name : str, optional
        The name of the Matplotlib colormap. If None, centered
        use seismic, otherwise viridis.

    Returns:
    --------
    mapper : matplotlib.cm.ScalarMappable
        An object with a .to_rgba(value) method to fetch colors.
    """
    if vcenter is not None:
        hr = max(abs(vmax - vcenter), abs(vcenter - vmin))
        norm = mcolors.CenteredNorm(
            vcenter=vcenter, halfrange=hr)
    else:
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    if vcenter is None:
        cmap_name = "viridis"
    else:
        cmap_name = "bwr"

    cmap = mpl.colormaps.get_cmap(cmap_name)
    return cm.ScalarMappable(norm=norm, cmap=cmap)



def adjust_plot_pos(pos, node_sizes_in, fig, ax,
                    padding_in=0.1, iterations=200, force_step_in=0.05,
                    seed=42):

    rng = np.random.default_rng(seed)
    adjusted_pos = {
        node: np.array(coord, dtype=float)
        for node, coord in pos.items()
    }
    nodes = list(adjusted_pos.keys())

    fig.canvas.draw()

    def data_to_inches(coord_pair):
        pixels = ax.transData.transform(coord_pair)
        return fig.dpi_scale_trans.inverted().transform(pixels)

    def inches_to_data(inch_pair):
        pixels = fig.dpi_scale_trans.transform(inch_pair)
        return ax.transData.inverted().transform(pixels)

    for _ in range(iterations):
        overlap_found = False
        rng.shuffle(nodes)

        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                n1, n2 = nodes[i], nodes[j]

                p1_in = data_to_inches(adjusted_pos[n1])
                p2_in = data_to_inches(adjusted_pos[n2])

                min_safe_dist_in = (
                    (node_sizes_in[n1] / 2.0)
                    + (node_sizes_in[n2] / 2.0)
                    + padding_in
                )

                direction_in = p2_in - p1_in
                distance_in = np.linalg.norm(direction_in)

                if distance_in == 0:
                    direction_in = rng.uniform(-0.5, 0.5, size=2)
                    distance_in = np.linalg.norm(direction_in)
                    if distance_in == 0:
                        direction_in = np.array([0.01, 0.00])
                        distance_in = 0.01

                if distance_in < min_safe_dist_in:
                    overlap_found = True
                    push_vector_in = (
                        direction_in
                        / distance_in
                        * force_step_in
                    )

                    p2_in += push_vector_in
                    p1_in -= push_vector_in

                    adjusted_pos[n2] = inches_to_data(p2_in)
                    adjusted_pos[n1] = inches_to_data(p1_in)

        if not overlap_found:
            break

    return adjusted_pos


def reset_xy_lim(ax, pos, ratio):
    all_x = [p[0] for p in pos.values()]
    all_y = [p[1] for p in pos.values()]
    x_expand = (max(all_x) - min(all_x)) * ratio
    y_expand = (max(all_y) - min(all_y)) * ratio
    ax.set_xlim(min(all_x) - x_expand, max(all_x) + x_expand)
    ax.set_ylim(min(all_y) - y_expand, max(all_y) + y_expand)



def plot_pw_subgraph(G, pos,
                     nodes, node_sizes,
                     node_data,
                     node_cmap, node_alpha,
                     node_cbar_title,
                     node_labels,
                     edges, edge_widths,
                     title,
                     figsize=None, seed=42):
    # node_sizes is {n: s}
    pos = {n: pos[n] for n in nodes}

    if figsize is None:
        figsize = (3, 4)

    fig, (ax_cbar, ax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=figsize,
        gridspec_kw={
            'height_ratios': (1, 12.5)
        }
    )

    reset_xy_lim(ax, pos, 0.5)

    adjusted_pos = adjust_plot_pos(
        pos, node_sizes, fig, ax,
        padding_in=0.8, seed=42
    )
    reset_xy_lim(ax, adjusted_pos, 0.1)

    for node in nodes:
        node_data_coord = adjusted_pos[node]
        node_pixel_coord = ax.transData.transform(node_data_coord)
        node_inch_coord = fig.dpi_scale_trans.inverted().transform(
            node_pixel_coord)

        size_in = node_sizes[node]

        ax_pie = inset_axes(
            ax,
            width=size_in,
            height=size_in,
            loc="center",
            bbox_to_anchor=node_inch_coord,
            bbox_transform=fig.dpi_scale_trans,
            axes_kwargs={'zorder': -1},
            borderpad=0,
        )
        ax_pie.set_facecolor("none")
        ax_pie.axis("off")

        col_data = node_data[node]

        ax_pie.pie(
            [1] * len(col_data),
            colors=[
                node_cmap.to_rgba(
                    x, alpha=node_alpha)
                for x in col_data
            ],
            radius=1.0,
            startangle=90,
        )

        outer_circle = mpl.patches.Circle(
            (0, 0),
            radius=1.05,
            edgecolor="grey",
            facecolor="none",
            linewidth=1,
            zorder=0,
        )
        ax_pie.add_patch(outer_circle)
        ax_pie.set_xlim(-1.2, 1.2)
        ax_pie.set_ylim(-1.2, 1.2)
        ax_pie.set_aspect("equal")

    edge_node_sizes = [
        (n, (x * 1.1 * 72 * 0.9) ** 2)
        for n, x in node_sizes.items()
    ]
    nx.draw_networkx_edges(
        G, adjusted_pos, edgelist=edges,
        ax=ax,
        edge_color="gray",
        width=edge_widths,
        alpha=0.7,
        arrowsize=5,
        node_size=[x[1] for x in edge_node_sizes],
        nodelist=[x[0] for x in edge_node_sizes]
    )

    label_dict = {
        n: textwrap.fill(
            l.replace('/', '/ ')
            , 20, break_long_words=False,
            break_on_hyphens=False
        ).replace('/ ', '/')
        for n, l in node_labels.items()
    }
    nx_labels = nx.draw_networkx_labels(
        G, adjusted_pos, label_dict, font_size=10,
        ax=ax,
        font_color="black",
        horizontalalignment="center")
    cbar = fig.colorbar(
        node_cmap, cax=ax_cbar,
        orientation="horizontal",
        alpha=node_alpha)

    ax_cbar.set_title(
        "{}\n{}".format(title, node_cbar_title),
        pad=5, fontsize=11,
        x=0.0, ha='left')

    for t in nx_labels.values():
        t.set_zorder(10)

    np.random.seed(seed)
    random.seed(seed)
    adjust_text(
        list(nx_labels.values()), ax=ax,
        force_text=(0, 0.5),
        force_pull=(0, 0.2),
        force_explode=(0, 0.2),
        force_static=(0, 0.2),
        iter_lim=1e3,
        only_move={
            "text": "y",
            "static": "y",
            "explode": "y",
            "pull": "y",
        }
    )
    ax.axis("off")
    return fig, ax




def plot_pw_neighborhood(rawk, s_pw, filename, node_color_attr,
                         mean_node_color_attr=True,
                         node_color_vmin=None, node_color_vmax=None,
                         node_color_center=None,
                         n_cutoff=10, node_alpha=0.5,
                         node_color_title=None,
                         title=None,
                         figsize=None,
                         draw_labels=True):
    """
    Plot pathway reachable local neighborhood

    Parameters
    ----------
    rawk : Rawk
        A Rawk instance.
    s_pw : str
        A pathway to plot.
    filename : str
        Saved plot filename
    node_color_attr : str
        Node addtribute used for coloring:
        - 'property': node property
        - 'rwv_cpm': random walk visit count per million
        - 'log10p1_rwv_cpm': log10(rwv_cpm + 1)
        - 'rwv_cpm_diff': sample - background rwv_cpm
        - 'rwv_cpm_p1_log2_fc': log2(
          (sample_rwv_cpm + 1) / (background_rwv_cpm + 1))
    mean_node_color_attr : bool
        If True, plot the mean node_color_attr of
        the reactions in each pathway node. If False,
        plot each pathway node as a pie chart, with each
        wedge as the node_color_attr of each reaction in
        the pathway.
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
    node_color_title : str, optional
        Title of the node colorbar.
    title : str, optional
        Title of the graph.
    figsize : (float, float)
        Figure size.

    Returns
    -------
    None
    """
    G = rawk.contrast_pw_graph
    pos = dict([(k, v["pw_pos"]) for k, v in G.nodes.items()])

    # {pw: [rxn1, rxn2, ...], ...}
    pw_rxns_dict = {}
    for i in range(rawk.sample.node_df.shape[0]):
        i_rxn = rawk.sample.node_df.iloc[i, :].loc["rxn"]
        i_pw = rawk.sample.node_df.iloc[i, :].loc["pathway"]

        if i_pw not in pw_rxns_dict:
            pw_rxns_dict[i_pw] = []

        pw_rxns_dict[i_pw].append(i_rxn)

    # random walk visit percentages matrices
    spr = (
        rawk.sample.n2v_num_walks
        * rawk.sample.n2v_walk_length
    )
    bg_spr = (
        rawk.bg_sample.n2v_num_walks
        * rawk.bg_sample.n2v_walk_length
    )
    assert spr == bg_spr
    rw_pw_pw_df = rawk.sample.rw_csp_rtp_spct
    rw_rxn_rxn_df = rawk.sample.rw_csr_rtr_sc / spr * 100
    bg_rw_rxn_rxn_df = rawk.bg_sample.rw_csr_rtr_sc / spr * 100
    rxn_g = rawk.sample.rxn_graph

    # Only implemented one source pw plot
    assert isinstance(s_pw, str)
    s_pw = [s_pw]
    s_pw_rxns_dict = {}
    for k, v in pw_rxns_dict.items():
        s_pw_rxns_dict[k] = [
            x for x in v if x in rw_rxn_rxn_df.columns]

    s_pw_rxns = [s_pw_rxns_dict[x] for x in s_pw]

    def get_topn_a0_inds(s, n, return_set=False):
        r = [k for k, v in s.sort_values(ascending=False)[:n].items()
             if v > 0]
        rset = set(r)
        assert len(r) == len(rset)
        if return_set:
            t = rset
        else:
            t = r
        return t

    s_reached_nodes = get_topn_a0_inds(
        rw_pw_pw_df.loc[:, s_pw[0]], n_cutoff)

    # The first step should always be in the pw.
    assert s_pw[0] in s_reached_nodes

    edge_width_by_edge = {}
    for i in s_reached_nodes:
        i_reached_nodes = get_topn_a0_inds(
            rw_pw_pw_df.loc[:, i], n_cutoff,
            return_set=True)

        for j_rank, j in enumerate(s_reached_nodes):
            if (i, j) in G.edges:
                if j in i_reached_nodes:
                    edge_width_by_edge[(i, j)] = (
                        0.15 + 0.08 * max(0, n_cutoff - j_rank)
                    )

    # Keep only nodes with one or more edges that are not self-loops.
    nodes = set()
    for i in edge_width_by_edge.keys():
        if i[0] != i[1]:
            nodes.add(i[0])
            nodes.add(i[1])

    nodes = sorted(nodes)
    edges = []
    edge_widths = []
    for x, w in edge_width_by_edge.items():
         if (x[0] in nodes) and (x[1] in nodes) and (x[0] != x[1]):
            edges.append(x)
            edge_widths.append(w)

    if len(nodes) == 0:
        nodes = [s_pw[0]]
    assert len(nodes) == len(set(nodes))
    assert len(edges) == len(set(edges))

    # Prepare colors
    node_color_dict = {}
    # convert random walk visit percentage to steps per million
    rw_pw_rxn_s = (
        1e4
        * rw_rxn_rxn_df.loc[:, s_pw_rxns[0]].mean(axis=1)
    )
    bg_rw_pw_rxn_s = (
        1e4
        * bg_rw_rxn_rxn_df.loc[
            rw_rxn_rxn_df.index.tolist(),
            s_pw_rxns[0]].mean(axis=1)
    )
    assert (
        rw_pw_rxn_s.index.tolist()
        == bg_rw_pw_rxn_s.index.tolist()
    )
    # nodes are reached pws from s_pw.
    # - plot placeholder nodes.
    # - plot edges.
    # - plot pie charts
    for i in nodes:
        if node_color_attr == "property":
            # list of all rxn properties in the pw
            i_col = G.nodes[i][node_color_attr]
        else:
            if node_color_attr == "rwv_cpm":
                # rw_pw_rxn_s is pw->rxn percent steps Series
                nc_s = rw_pw_rxn_s
            elif node_color_attr == "log10p1_rwv_cpm":
                # log10 steps percent + 1
                nc_s = np.log10(rw_pw_rxn_s + 1)
            elif node_color_attr == "rwv_cpm_diff":
                nc_s = rw_pw_rxn_s - bg_rw_pw_rxn_s
            elif node_color_attr == "rwv_cpm_p1_log2fc":
                nc_s = np.log2((rw_pw_rxn_s + 1) / (bg_rw_pw_rxn_s + 1))
            else:
                raise ValueError(
                    "Unknown node_color_attr: {}".format(node_color_attr))
            i_col = nc_s.loc[s_pw_rxns_dict[i]].tolist()

        i_col = sorted(i_col)
        if len(i_col) == 0:
            i_col = [0]
        node_color_dict[i] = i_col

    if mean_node_color_attr:
        node_color_dict = {
            n: [np.mean(c)]
            for n, c in node_color_dict.items()
        }

    if node_color_vmin is None:
        node_color_vmin = np.min(
            [np.min(x) for x in node_color_dict.values()])

    if node_color_vmax is None:
        node_color_vmax = np.max(
            [np.max(x) for x in node_color_dict.values()])

    node_color_vmin, node_color_vmax = expand_v0range(
        node_color_vmin, node_color_vmax)

    if node_color_center is None:
        if node_color_attr in ["rwv_cpm_p1_log2fc", "rwv_cpm_diff"]:
            node_color_center = 0
            max_abs = max(
                abs(node_color_vmin),
                abs(node_color_vmax))
            node_color_vmin = -max_abs
            node_color_vmax = max_abs

    node_color_cmap = get_color_mapper(
        node_color_vmin, node_color_vmax,
        node_color_center)

    if draw_labels:
        node_labels = {n: n for n in nodes}
    else:
        node_labels = {n: "" for n in nodes}

    fig, ax = plot_pw_subgraph(
        G, pos, nodes,
        dict((n, 0.4) if n in s_pw else (n, 0.25) for n in nodes),
        node_color_dict, node_color_cmap, node_alpha,
        node_color_title,
        node_labels,
        edges, edge_widths,
        title,
        figsize=figsize)

    plt.savefig(
        filename, dpi=600, format="png",
        bbox_inches='tight',
        pad_inches=0.1)
    plt.close()


def map_strings_to_colors(strings, palette):
    unique_strings = pd.Series(strings).drop_duplicates().tolist()
    string_to_color_map = {
        s: palette[i % len(palette)]
        for i, s in enumerate(unique_strings)
    }
    return [string_to_color_map[s] for s in strings]



def plot_mtx(mtx, row_labels, col_labels, cbar_title,
             figsize, filename, cbar_hm_space,
             vmin=None, vmax=None):

    color_palette = ggsci_pal_jco
    row_colors = map_strings_to_colors(row_labels, color_palette)
    col_colors = map_strings_to_colors(col_labels, color_palette)

    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")
    gs = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=[1, 12],
        hspace=cbar_hm_space
    )
    hm_gs = mpl.gridspec.GridSpecFromSubplotSpec(
        nrows=2,
        ncols=2,
        subplot_spec=gs[1, 0],
        width_ratios=[0.3, 10],
        height_ratios=[0.3, 10],
        wspace=0,
        hspace=0
    )

    ax_cbar = fig.add_subplot(gs[0, 0])
    ax_col  = fig.add_subplot(hm_gs[0, 1])
    ax_row  = fig.add_subplot(hm_gs[1, 0])
    ax_heat = fig.add_subplot(hm_gs[1, 1])

    im = ax_heat.imshow(
        mtx, aspect="auto", cmap="viridis",
        interpolation="none",
        vmin=vmin, vmax=vmax)
    ax_heat.axis("off")

    col_rgb = [cm.colors.to_rgba(c) for c in col_colors]
    ax_col.imshow([col_rgb], aspect="auto")
    ax_col.axis("off")

    row_rgb = [[cm.colors.to_rgba(c)] for c in row_colors]
    ax_row.imshow(row_rgb, aspect="auto")
    ax_row.axis("off")

    cbar = fig.colorbar(
        im, cax=ax_cbar, orientation="horizontal")
    ax_cbar.set_title(
        cbar_title, pad=5, ha="left", x=0)

    plt.savefig(
        filename, dpi=600, format="png",
        bbox_inches='tight',
        pad_inches=0.1)
    plt.close()



def assert_allclose(x, y):
    np.testing.assert_allclose(
        x, y, rtol=0, atol=1e-6)



def pr_sort_key(x):
    if "_____" in x:
        k = x.split("_____")
    else:
        k = x
    return x


def prs_to_ps(prs):
    sep = "_____"
    ps = [x.split(sep)[0] if sep in x else x for x in prs]
    return ps



def plot_rawk_sample_mtx(rawk_sample, save_pfx, base_figsize=2,
                         vmin=None, vmax=None):
    rxn_pw_dict = rawk_sample.pw_by_rxn
    rxns = sorted(rxn_pw_dict.keys())

    spr = (
        rawk_sample.n2v_num_walks
        * rawk_sample.n2v_walk_length
    )

    # metabolic reaction network edge weight
    w = rawk_sample.rxn_to_pw_rxn(
        pd.DataFrame(
            nx.to_numpy_array(
                rawk_sample.rxn_graph,
                nodelist=rxns,
                weight="mn_weight").T,
            index=rxns, columns=rxns))

    assert w.index.is_unique
    assert w.columns.is_unique
    assert w.index.tolist() == w.columns.tolist()

    # dest node property
    dp = rawk_sample.rxn_to_pw_rxn(
        pd.DataFrame(
            nx.to_numpy_array(
                rawk_sample.rxn_graph,
                nodelist=rxns,
                weight="n2v_weight").T,
            index=rxns, columns=rxns))

    tp = rawk_sample.rxn_to_pw_rxn(
        pd.DataFrame(
            nx.to_numpy_array(
                rawk_sample.rxn_graph,
                nodelist=rxns,
                weight="n2v_transition_prob").T,
            index=rxns, columns=rxns))

    # rr, rp, pr, pp. (tgt, src)
    # (tgt rxn, src rxn)
    rr = rawk_sample.rw_cspr_rtpr_sc.copy()
    assert rr.index.is_unique
    assert rr.columns.is_unique
    assert rr.index.tolist() == rr.columns.tolist()
    assert rr.index.tolist() == w.columns.tolist()

    p_by_pr = rawk_sample.p_by_pr

    assert_allclose(rr.sum(), spr)
    rr_cpm = rr / spr * 1e6

    # (tgt pw, src rxn)
    # row sum
    rp = rr.copy().groupby(p_by_pr).sum()
    assert_allclose(rp.sum(), spr)
    assert_frame_equal(
        rp, rawk_sample.rw_cspr_rtp_sc)
    rp_cpm = rp / spr * 1e6

    # (tgt rxn, src pw)
    # col pw mean
    pr = rr.copy().T.groupby(p_by_pr).mean().T
    assert_allclose(pr.sum(), spr)
    pr_cpm = pr / spr * 1e6
    # pp
    pp = pr.copy().groupby(p_by_pr).sum()
    assert_allclose(pp.sum(), spr)
    pp_cpm = pp / spr * 1e6
    assert_frame_equal(
        pp_cpm, rawk_sample.rw_csp_rtp_spct * 1e4)

    plot_config_by_id = {
        "mrn_edge_weight": (
            np.log10(w + 1),
            "reactions", "reactions",
            "log10(edge weight + 1)",
            (base_figsize+1.7, base_figsize+2),
        ),

        "dest_property": (
            dp,
            "reactions", "reactions",
            "destination node property",
            (base_figsize+1.7, base_figsize+2),
        ),

        "n2v_transition_prob": (
            tp,
            "reactions", "reactions",
            "node2vec transition probabity",
            (base_figsize+1.7, base_figsize+2),
        ),

        "rr_rwvcpm": (
            np.log10(rr_cpm + 1),
            "reactions", "reactions",
            "log10(RWVCPM + 1)",
            (base_figsize+1.7, base_figsize+2),
        ),

        "rp_rwvcpm": (
            np.log10(rp_cpm + 1),
            "reactions", "pathways",
            "log10(RWVCPM + 1)",
            (base_figsize+1.7, base_figsize),
        ),

        "pr_rwvcpm": (
            np.log10(pr_cpm + 1),
            "pathways", "reactions",
            "log10(RWVCPM + 1)",
            (base_figsize+0.9, base_figsize+2),
        ),

        "pp_rwvcpm": (
            np.log10(pp_cpm + 1),
            "pathways", "pathways",
            "log10(RWVCPM + 1)",
            (base_figsize+0.9, base_figsize+0.9),
        ),
    }

    for pid, pc in plot_config_by_id.items():
        df = pc[0].copy()
        assert df.index.is_unique
        assert df.columns.is_unique
        df = df.sort_index(key=pr_sort_key).sort_index(
            axis=1, key=pr_sort_key)
        fn = f"{save_pfx}_{pid}_mtx.png"
        cbar_title = (
            f"Columns: source {pc[1]}\n"
            f"Rows: target {pc[2]}\n"
            f"{pc[3]}"
        )
        plot_mtx(
            df.values,
            prs_to_ps(df.index.tolist()),
            prs_to_ps(df.columns.tolist()),
            cbar_title,
            pc[4],
            fn,
            0.8 / pc[4][1],
            vmin=vmin, vmax=vmax,
        )
