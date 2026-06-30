import numpy as np
from pecanpy import pecanpy
from pecanpy.graph import AdjlstGraph
from scipy.sparse import csr_matrix



class RawkNode2Vec:
    """
    Mem efficient chunked Node2Vec
    """
    def __init__(self, node_df, ea_df,
                 p=1.0, q=1.0,
                 num_walks=8000,
                 chunk_num_walks=1000,
                 walk_length=20,
                 workers=1, seed=17,
                 keep_walks=False):

        self.node_df = node_df
        self.ea_df = ea_df

        unique_nodes = sorted(self.node_df["rxn"].tolist())
        assert len(unique_nodes) == len(set(unique_nodes))

        num_nodes = len(unique_nodes)

        node_to_idx = {name: i for i, name in enumerate(unique_nodes)}
        idx_to_node = {i: name for i, name in enumerate(unique_nodes)}

        # use int index and str(int) node names for efficiency
        # node name <> idx <> str(idx). one-to-one mappings
        src_idx = self.ea_df["src"].map(node_to_idx)
        assert src_idx.notna().all()
        src_idx = src_idx.values

        tgt_idx = self.ea_df["dest"].map(node_to_idx)
        assert tgt_idx.notna().all()
        tgt_idx = tgt_idx.values

        weights = self.ea_df["n2v_weight"]
        assert weights.notna().all()
        weights = weights.values

        adj_graph = AdjlstGraph()

        for i in range(len(unique_nodes)):
            adj_graph.add_node(str(i))

        assert len(src_idx) == len(tgt_idx) == len(weights)

        for i in range(len(src_idx)):
            adj_graph.add_edge(
                str(src_idx[i]),
                str(tgt_idx[i]),
                weights[i],
                directed=True)

        self.g = pecanpy.SparseOTF.from_adjlst_graph(
            adj_graph, p=p, q=q,
            workers=workers, random_state=seed)

        self.node_to_idx = node_to_idx
        self.idx_to_node = idx_to_node

        assert len(self.g.nodes) == num_nodes

        self.num_nodes = num_nodes

        self.visit_matrix = csr_matrix(
            (self.num_nodes, self.num_nodes),
            dtype=np.uint32)

        self.num_walks = num_walks
        self.walk_length = walk_length
        self.chunk_num_walks = chunk_num_walks
        self.keep_walks = keep_walks

        if keep_walks:
            self.walks = []
        else:
            self.walks = None

        # only count visits in the first epoch
        self.curr_epoch = 0

    def __iter__(self):
        """
        Generate one chunk of walks
        """
        num_chunks = int(np.ceil(self.num_walks / self.chunk_num_walks))
        self.curr_epoch += 1

        for i in range(num_chunks):
            # Compute walk size allocations
            curr_num_walks = min(
                self.chunk_num_walks,
                # num walks left
                self.num_walks - (i * self.chunk_num_walks))

            walks_chunk = self.g.simulate_walks(
                num_walks=curr_num_walks,
                # walk_length - 1, to be compatible with node2vec
                # package. pecanpy does not count start node as
                # the first step.
                walk_length=self.walk_length - 1)

            if self.keep_walks:
                self.walks += walks_chunk

            assert all([len(x) == self.walk_length for x in walks_chunk])
            assert len(walks_chunk) == curr_num_walks * self.num_nodes

            if self.curr_epoch == 1:
                # walks_chunk are str node indexes
                walks_chunk_arr = np.array(walks_chunk, dtype=np.int32)
                assert (
                    walks_chunk_arr.shape
                    == (curr_num_walks * self.num_nodes,
                        self.walk_length)
                )

                # start node
                # [s1, s1, ..., s1, s2, s2, ...]
                # repeat each element for walk_length times
                # (dest_rxn, src_rxn) count
                visit_cols = np.repeat(
                    walks_chunk_arr[:, 0], self.walk_length)

                # walks_chunk_arr
                # each row is a walk from a start node
                # ravel default order is C.
                #
                # 'C' means to index the elements in row-major, C-style order,
                # with the last axis index changing fastest, back to the first
                # axis index changing slowest.
                # >>> np.array([[1, 2, 3], [11, 12, 13]]).ravel()
                # array([ 1,  2,  3, 11, 12, 13])
                visit_rows = walks_chunk_arr.ravel()

                ones_data = np.ones(len(visit_cols), dtype=np.int32)

                assert len(ones_data) == len(visit_rows) == len(visit_cols)

                chunk_matrix = csr_matrix(
                    (ones_data, (visit_rows, visit_cols)),
                    shape=(self.num_nodes, self.num_nodes),
                    dtype=np.uint32
                )

                self.visit_matrix += chunk_matrix
                del chunk_matrix

            for walk in walks_chunk:
                yield walk

            del walks_chunk

    def get_visit_arr(self):
        # (tgt, src) num steps count matrix
        visit_arr = self.visit_matrix.toarray()
        assert np.all(
            visit_arr.sum(axis=0) == self.walk_length * self.num_walks)
        assert visit_arr.shape == (self.num_nodes, self.num_nodes)
        return visit_arr
