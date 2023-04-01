import pandas as pd
import numpy as np
import networkx as nx
from copy import deepcopy
import sys

# Load data
expressions = pd.read_csv("CCLE_expression.csv", low_memory=False)
mutations = pd.read_csv("CCLE_mutations.csv", low_memory=False)
effective_mutations = pd.read_csv("Effective_mutations.csv", low_memory=False)
complexes_in = pd.read_csv("In_complexes.csv", low_memory=False)
pathways_in = pd.read_csv("In_pathways.csv", low_memory=False)

# Effects of mutations (gene, protein_change)
mutation_effects = list(zip(list(effective_mutations["Hugo_Symbol"]), list(effective_mutations["Protein_Change"])))
# Changed aa position on protein due to mutation
is_deleterious = list(effective_mutations["isDeleterious"])
# Effected interaction partners due to mutation on the hotspot
effected_interaction = list(effective_mutations["effected_interaction"])
# Complexes
in_complexes = dict(zip(list(complexes_in["Gene"]), list(complexes_in["ComplexIDs"])))
# PathwaysComplexIDs
in_pathways = dict(zip(list(pathways_in["Genes"]), list(pathways_in["Pathways"])))


def curate_mutations(Graph):
    """
    :param Graph: Graph to curate mutations.
    :return: list of all mutations as list of tuples. ex: [(TP53, p.S245A), (AR, p.S245fs) etc.]
    """
    # All mutations of given cell line (gene, protein_change)
    all_mutations_tuples = set(zip(list(mutations["Hugo_Symbol"][mutations["DepMap_ID"] == Graph.name]),
                                   list(mutations["Protein_Change"][mutations["DepMap_ID"] == Graph.name])))
    all_mutations_tuples = list(all_mutations_tuples.intersection(mutation_effects))
    return all_mutations_tuples


def apply_single_mutation(Graph, mutation):
    """
    :param Graph: Graph to be updated!
    :param mutation: (HGNC symbol, protein change coordinate) tuple. ex: (TP53, p.S245A)
    :return: None
    """
    effected_paths_comps = []  # List of sets, each set contains the common pathways and complexes of first broken edge.
    effected = []  # First effected targets neighbors of the mutated gene.
    # We determine the first cycle of edges to be inactivated and their common paths and comps
    counter = 0  # Counts the number of effected edges.
    if mutation[0] in Graph.nodes:
        if mutation in mutation_effects:
            if effective_mutations["effected_aa"][mutation_effects.index(mutation)] == -1:  # If mutation is deleterious
                for edge in Graph.out_edges(mutation[0]):  # all outgoing edges of mutated gene are candidate
                    if Graph.edges[edge]["is_directed"] != 0 or Graph.edges[edge]["BIOGRID"] \
                            or Graph.edges[edge]["pathway"] or Graph.edges[edge]["complex"]:  # If edge is effective
                        # complexes of source node
                        if mutation[0] in in_complexes.keys():
                            source_complexes = set(in_complexes[mutation[0]].split(";"))
                        else:
                            source_complexes = set([])
                        # complexes of target node
                        if edge[1] in in_complexes.keys():
                            target_complexes = set(in_complexes[edge[1]].split(";"))
                        else:
                            target_complexes = set([])
                        # Common complexes
                        common_complexes = source_complexes.intersection(target_complexes)
                        # pathways of source node
                        if mutation[0] in in_pathways.keys():
                            source_pathways = set(in_pathways[mutation[0]].split(";"))
                        else:
                            source_pathways = set([])
                        # pathways of target node
                        if edge[1] in in_pathways.keys():
                            target_pathways = set(in_pathways[edge[1]].split(";"))
                        else:
                            target_pathways = set([])
                        # Common pathways
                        common_pathways = source_pathways.intersection(target_pathways)
                        if len(common_pathways.union(common_complexes)) > 0:  # If there are common pathways or complexes
                            effected_paths_comps.append(common_pathways.union(common_complexes))
                            effected.append(edge[1])
                            nx.set_edge_attributes(Graph, {edge: 0}, "is_directed")
                            nx.set_edge_attributes(Graph, {edge: False}, "BIOGRID")
                            nx.set_edge_attributes(Graph, {edge: False}, "pathway")
                            nx.set_edge_attributes(Graph, {edge: False}, "complex")

            else:  # If mutation is not deleterious but possibly hits an interface
                if pd.notna(effected_interaction[mutation_effects.index(mutation)]):  # is there an effected neighbor!
                    neighs = effected_interaction[mutation_effects.index(mutation)].split(
                        ";")  # list of effected neighbors
                    for neigh in neighs:
                        edge = (mutation[0], neigh)
                        if edge in Graph.edges:
                            if Graph.edges[edge]["is_directed"] != 0 or Graph.edges[edge]["BIOGRID"] or Graph.edges[edge]["pathway"] or Graph.edges[edge]["complex"]:  # If edge is effective:
                                # complexes of source node
                                if mutation[0] in in_complexes.keys():
                                    source_complexes = set(in_complexes[mutation[0]].split(";"))
                                else:
                                    source_complexes = set([])
                                # complexes of target node
                                if neigh in in_complexes.keys():
                                    target_complexes = set(in_complexes[neigh].split(";"))
                                else:
                                    target_complexes = set([])
                                # Common complexes
                                common_complexes = source_complexes.intersection(target_complexes)
                                # pathways of source node
                                if mutation[0] in in_pathways.keys():
                                    source_pathways = set(in_pathways[mutation[0]].split(";"))
                                else:
                                    source_pathways = set([])
                                # pathways of target node
                                if neigh in in_pathways.keys():
                                    target_pathways = set(in_pathways[neigh].split(";"))
                                else:
                                    target_pathways = set([])
                                # Common pathways
                                common_pathways = source_pathways.intersection(target_pathways)
                                if len(common_pathways.union(common_complexes)) > 0:  # If there are common pathways or complexes
                                    effected_paths_comps.append(common_pathways.union(common_complexes))
                                    effected.append(neigh)
                                    nx.set_edge_attributes(Graph, {edge: 0}, "is_directed")
                                    nx.set_edge_attributes(Graph, {edge: False}, "BIOGRID")
                                    nx.set_edge_attributes(Graph, {edge: False}, "pathway")
                                    nx.set_edge_attributes(Graph, {edge: False}, "complex")
                                    reverse_edge = (edge[1], edge[0])
                                    if reverse_edge in Graph.edges:
                                        nx.set_edge_attributes(Graph, {reverse_edge: 0}, "is_directed")
                                        nx.set_edge_attributes(Graph, {reverse_edge: False}, "BIOGRID")
                                        nx.set_edge_attributes(Graph, {reverse_edge: False}, "pathway")
                                        nx.set_edge_attributes(Graph, {reverse_edge: False}, "complex")

            #  Downstream flow of mutation effects
            sources = [[i] for i in effected]
            counter = len(effected)
            while len(np.array(deepcopy(sources), dtype=object).flatten()) != 0:
                new_sources = [[] for j in range(len(effected))]
                for t, source in enumerate(sources):
                    for s in source:
                        for edge in Graph.out_edges(s):
                            if Graph.edges[edge]["is_directed"] != 0 or Graph.edges[edge]["BIOGRID"] or Graph.edges[edge]["pathway"] or Graph.edges[edge]["complex"]:
                                if s in in_complexes.keys():
                                    s_comps = set(in_complexes[s].split(";"))
                                else:
                                    s_comps = set([])
                                if edge[1] in in_complexes.keys():
                                    t_comps = set(in_complexes[edge[1]].split(";"))
                                else:
                                    t_comps = set([])
                                if s in in_pathways.keys():
                                    s_paths = set(in_pathways[s].split(";"))
                                else:
                                    s_paths = set([])
                                if edge[1] in in_pathways.keys():
                                    t_paths = set(in_pathways[edge[1]].split(";"))
                                else:
                                    t_paths = set([])
                                c_comps = s_comps.intersection(t_comps)
                                c_paths = s_paths.intersection(t_paths)
                                u_paths_comps = c_comps.union(c_paths)
                                if len(effected_paths_comps[t].intersection(u_paths_comps)) > 0:
                                    nx.set_edge_attributes(Graph, {edge: 0}, "is_directed")
                                    nx.set_edge_attributes(Graph, {edge: False}, "BIOGRID")
                                    nx.set_edge_attributes(Graph, {edge: False}, "pathway")
                                    nx.set_edge_attributes(Graph, {edge: False}, "complex")
                                    reverse_edge = (edge[1], edge[0])
                                    if reverse_edge in Graph.edges:
                                        nx.set_edge_attributes(Graph, {reverse_edge: 0}, "is_directed")
                                        nx.set_edge_attributes(Graph, {reverse_edge: False}, "BIOGRID")
                                        nx.set_edge_attributes(Graph, {reverse_edge: False}, "pathway")
                                        nx.set_edge_attributes(Graph, {reverse_edge: False}, "complex")
                                    counter += 1
                                    new_sources[t].append(edge[1])
                sources = deepcopy(new_sources)
    print(f"{mutation} --> Effected number of edges:", counter)
    return counter


def determine_isolated_nodes(G):
    isolated_nodes = []
    switch = True
    for node in G.nodes:
        for edge in list(G.out_edges(node)) + list(G.in_edges(node)):
            if G.edges[edge]["is_directed"] != 0 or G.edges[edge]["BIOGRID"] or G.edges[edge]["pathway"] or \
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
        return None


def estimator(gene1, gene2, current_exp, type, window_size=0.05, exp_thr=0.585, pct=0.1, sensitivity=0.01, min_samples=50):
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


# Simple progress bar
def progres(current, total, status='', bar_len=50):
    filled_len = int(round(bar_len * current / float(total)))

    percents = round(100.0 * current / float(total), 1)
    bar = '#' * filled_len + ' ' * (bar_len - filled_len)

    fmt = '[%s] %s%s ...%s' % (bar, percents, '%', status)
    print('\b' * len(fmt), end='')  # clears the line
    sys.stdout.write(fmt)
    sys.stdout.flush()
