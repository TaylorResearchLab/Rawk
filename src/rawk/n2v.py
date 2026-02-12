# Adapted from GitHub eliorc/node2vec repository.
import os
import random
from collections import defaultdict

import numpy as np
import networkx as nx
from node2vec import Node2Vec



class NoWalkNode2Vec(Node2Vec):
    """
    Create a Node2Vec without generating random walks
    """
    def __init__(self, graph: nx.Graph, dimensions: int = 128,
                 walk_length: int = 80, num_walks: int = 10, p: float = 1,
                 q: float = 1, weight_key: str = 'weight', workers: int = 1,
                 sampling_strategy: dict = None,
                 quiet: bool = False, temp_folder: str = None,
                 seed: int = None):
        self.graph = graph
        self.dimensions = dimensions
        self.walk_length = walk_length
        self.num_walks = num_walks
        self.p = p
        self.q = q
        self.weight_key = weight_key
        self.workers = workers
        self.quiet = quiet
        self.d_graph = defaultdict(dict)

        if sampling_strategy is None:
            self.sampling_strategy = {}
        else:
            self.sampling_strategy = sampling_strategy

        self.temp_folder, self.require = None, None
        if temp_folder:
            if not os.path.isdir(temp_folder):
                raise NotADirectoryError(
                    "temp_folder does not exist or is not "
                    "a directory. ({})".format(temp_folder))

            self.temp_folder = temp_folder
            self.require = "sharedmem"

        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

        self._precompute_probabilities()
        self.walks = None
