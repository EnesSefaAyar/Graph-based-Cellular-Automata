import numpy as np
import networkx as nx
from copy import deepcopy
import random


def apply_single_mutation(Graph, mutation, complex_dict):
    """
    :param Graph: Graph to be updated!
    :param mutation: (HGNC symbol, protein change coordinate) tuple. ex: (TP53, p.S245A)
    :return: None
    """
    effected_comps = []  # List of sets, each set contains the common pathways and complexes of first broken edge.
    effected = []  # First effected targets neighbors of the mutated gene.
    # We determine the first cycle of edges to be inactivated and their common paths and comps
    effected_edges = [] # Keeps the effected edges from mutations
    if mutation in Graph.nodes:
        if random.random() < 0.5:  # If mutation is deleterious
            for edge in Graph.out_edges(mutation):  # all outgoing edges of mutated gene are candidate
                if Graph.edges[edge]["is_directed"] != 0 or Graph.edges[edge]["BIOGRID"] \
                         or Graph.edges[edge]["complex"]:  # If edge is effective
                    # complexes of source node
                    try:
                        source_complexes = set(complex_dict[mutation])
                    except:
                        source_complexes = set([])
                    # complexes of target node
                    try:
                        target_complexes = set(complex_dict[edge[1]])
                    except:
                        target_complexes = set([])
                    # Common complexes
                    common_complexes = source_complexes.intersection(target_complexes)
                    if len(common_complexes) > 0:  # If there are common pathways or complexes
                        effected_comps.append(common_complexes)
                        effected.append(edge[1])
                        nx.set_edge_attributes(Graph, {edge: 0}, "is_directed")
                        nx.set_edge_attributes(Graph, {edge: False}, "BIOGRID")
                        nx.set_edge_attributes(Graph, {edge: False}, "complex")
                        effected_edges.append(edge)
        elif random.random() < 0.9:  # If mutation is not deleterious but possibly hits an interface
            neighs = []
            for edge in Graph.out_edges(mutation):
                if random.random() < 1 / len(Graph.out_edges(mutation)):
                    neighs.append(edge[1])
            for neigh in neighs:
                edge = (mutation, neigh)
                if edge in Graph.edges:
                    if Graph.edges[edge]["is_directed"] != 0 or Graph.edges[edge]["BIOGRID"] \
                            or Graph.edges[edge][
                        "complex"]:  # If edge is effective:
                        # complexes of source node
                        try:
                            source_complexes = set(complex_dict[mutation])
                        except:
                            source_complexes = set([])
                        try:
                            target_complexes = set(complex_dict[neigh])
                        except:
                            target_complexes = set([])
                        # Common complexes
                        common_complexes = source_complexes.intersection(target_complexes)
                        if len(common_complexes) > 0:  # If there are common pathways or complexes
                            effected_comps.append(common_complexes)
                            effected.append(neigh)
                            nx.set_edge_attributes(Graph, {edge: 0}, "is_directed")
                            nx.set_edge_attributes(Graph, {edge: False}, "BIOGRID")
                            nx.set_edge_attributes(Graph, {edge: False}, "complex")
                            reverse_edge = (edge[1], edge[0])
                            if reverse_edge in Graph.edges:
                                nx.set_edge_attributes(Graph, {reverse_edge: 0}, "is_directed")
                                nx.set_edge_attributes(Graph, {reverse_edge: False}, "BIOGRID")
                                nx.set_edge_attributes(Graph, {reverse_edge: False}, "complex")
                                effected_edges.append(edge)
                                effected_edges.append(reverse_edge)
                            else:
                                effected_edges.append(edge)

    #  Downstream flow of mutation effects
    sources = [[i] for i in effected]
    while len(np.array(deepcopy(sources), dtype=object).flatten()) != 0:
        new_sources = [[] for j in range(len(effected))]
        for t, source in enumerate(sources):
            for s in source:
                for edge in Graph.out_edges(s):
                    if Graph.edges[edge]["is_directed"] != 0 or Graph.edges[edge]["BIOGRID"] \
                             or Graph.edges[edge]["complex"]:
                        try:
                            s_comps = set(complex_dict[s])
                        except:
                            s_comps = set([])
                        try:
                            t_comps = set(complex_dict[edge[1]])
                        except:
                            t_comps = set([])
                        c_comps = s_comps.intersection(t_comps)
                        if len(effected_comps[t].intersection(c_comps)) > 0:
                            nx.set_edge_attributes(Graph, {edge: 0}, "is_directed")
                            nx.set_edge_attributes(Graph, {edge: False}, "BIOGRID")
                            nx.set_edge_attributes(Graph, {edge: False}, "complex")
                            reverse_edge = (edge[1], edge[0])
                            if reverse_edge in Graph.edges:
                                nx.set_edge_attributes(Graph, {reverse_edge: 0}, "is_directed")
                                nx.set_edge_attributes(Graph, {reverse_edge: False}, "BIOGRID")
                                nx.set_edge_attributes(Graph, {reverse_edge: False}, "complex")
                                effected_edges.append(edge)
                                effected_edges.append(reverse_edge)
                            else:
                                effected_edges.append(edge)
                            new_sources[t].append(edge[1])
        sources = deepcopy(new_sources)
    print("Effected number of edges:", len(effected_edges), "Node number:", len(Graph.nodes))
    return  effected_edges


def determine_isolated_nodes(G):
    isolated_nodes = []
    switch = True
    for node in G.nodes:
        for edge in list(G.out_edges(node)) + list(G.in_edges(node)):
            if G.edges[edge]["is_directed"] != 0 or G.edges[edge]["BIOGRID"] or \
                    G.edges[edge]["complex"]:
                switch = False
                break
        if switch:
            isolated_nodes.append(node)
        else:
            switch = True

    if len(isolated_nodes) != 0:
        return isolated_nodes
    else:
        return []


def estimator(gene1, gene2, current_exp, type, expressions, window_size=0.05, exp_thr=0.585, pct=0.1, sensitivity=0.01,
              min_samples=50):
    x = np.array(expressions[gene1])
    y = np.array(expressions[gene2])
    t_index = np.array([False] * len(x))
    t_index += (x > exp_thr) & (y > exp_thr)
    if np.sum(t_index) > min_samples:
        x = x[t_index]
        y = y[t_index]
        index = np.array([False] * len(x))
        while np.sum(index) < pct * np.sum(t_index):
            index += (x > (current_exp - window_size)) & (x < (current_exp + window_size))
            window_size += sensitivity
        b = y[index]
        if window_size < 2 * exp_thr:
            return np.mean(b), np.std(b)
        else:
            if current_exp > np.mean(expressions[gene1]):
                if type == 1:
                    return np.mean(b) * (current_exp ** 0.01), np.std(b)
                else:
                    return np.mean(b), np.std(b)
            else:
                if type == 1:
                    return np.min(b), np.std(b)
                else:
                    return np.mean(b) / (current_exp ** 0.01), np.std(b)
