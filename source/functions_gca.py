import pandas as pd
import numpy as np
import networkx as nx
from copy import deepcopy
import random
import sys
import json

# Load data
mutations = pd.read_csv("raw_data/CCLE_mutations.csv", low_memory=False)
effective_mutations = pd.read_csv("Effective_mutations.csv", low_memory=False)
complexes_in = pd.read_csv("In_complexes.csv", low_memory=False)
pathways_in = pd.read_csv("In_pathways.csv", low_memory=False)

with open("estimator.json", "r") as f:
    estimation_dict = json.load(f)
estimated_genes = set(list(estimation_dict.keys()))

# Effects of mutations (gene, protein_change)
mutation_effects = list(zip(list(effective_mutations["Hugo_Symbol"]), list(effective_mutations["Protein_Change"])))
# Changed aa position on protein due to mutation
is_deleterious = list(effective_mutations["isDeleterious"])
# Effected interaction partners due to mutation on the hotspot
effected_interaction = list(effective_mutations["effected_interaction"])
# Complexes
in_complexes = dict(zip(list(complexes_in["Genes"]), list(complexes_in["ComplexIDs"])))
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


def apply_single_mutation(Graph, mutation, randomize=False):
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
                if randomize:
                    mutation = (random.choice(list(Graph.nodes)), mutation[1])
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
                        if len(common_pathways.union(
                                common_complexes)) > 0:  # If there are common pathways or complexes
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
                    if randomize:
                        mutation = (random.choice(list(Graph.nodes)), mutation[1])
                        new_neighs = []
                        count = 0
                        for edge in Graph.out_edges(mutation[0]):
                            if count < len(neighs):
                                if random.random() < 1 / len(Graph.out_edges(mutation[0])):
                                    new_neighs.append(edge[1])
                                    count += 1
                        neighs = new_neighs
                    for neigh in neighs:
                        edge = (mutation[0], neigh)
                        if edge in Graph.edges:
                            if Graph.edges[edge]["is_directed"] != 0 or Graph.edges[edge]["BIOGRID"] or \
                                    Graph.edges[edge]["pathway"] or Graph.edges[edge][
                                "complex"]:  # If edge is effective:
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
                                if len(common_pathways.union(
                                        common_complexes)) > 0:  # If there are common pathways or complexes
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
                            if Graph.edges[edge]["is_directed"] != 0 or Graph.edges[edge]["BIOGRID"] or \
                                    Graph.edges[edge]["pathway"] or Graph.edges[edge]["complex"]:
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


def estimator(gene1, gene2, current_exp):
    if gene1 + "_" + gene2 in estimated_genes:
        exact_exp = current_exp
        current_exp = int(current_exp)
        try:
            lower = estimation_dict[gene1 + "_" + gene2][str(current_exp)]
            upper = estimation_dict[gene1 + "_" + gene2][str(current_exp + 1)]
            return (exact_exp - current_exp) * (upper - lower) + lower
        except KeyError:
            limit = max(np.array(list(estimation_dict[gene1 + "_" + gene2].keys()), dtype=int))
            return estimation_dict[gene1 + "_" + gene2][str(limit)]
    else:
        return None


def directed_edge_swap(G, *, nswap=1, max_tries=100, seed=None):
    if nswap > max_tries:
        raise nx.NetworkXError("Number of swaps > number of tries allowed.")
    if len(G) < 4:
        raise nx.NetworkXError("DiGraph has fewer than four nodes.")
    if len(G.edges) < 3:
        raise nx.NetworkXError("DiGraph has fewer than 3 edges")

    # Instead of choosing uniformly at random from a generated edge list,
    # this algorithm chooses nonuniformly from the set of nodes with
    # probability weighted by degree.
    tries = 0
    swapcount = 0
    keys, degrees = zip(*G.degree())  # keys, degree
    cdf = nx.utils.cumulative_distribution(degrees)  # cdf of degree
    discrete_sequence = nx.utils.discrete_sequence

    while swapcount < nswap:
        # choose source node index from discrete distribution
        start_index = discrete_sequence(1, cdistribution=cdf, seed=seed)[0]
        start = keys[start_index]
        tries += 1

        if tries > max_tries:
            msg = f"Maximum number of swap attempts ({tries}) exceeded before desired swaps achieved ({nswap})."
            raise nx.NetworkXAlgorithmError(msg)

        # If the given node doesn't have any out edges, then there isn't anything to swap
        if G.out_degree(start) == 0:
            continue
        second = random.choice(list(G.succ[start]))
        if start == second:
            continue

        if G.out_degree(second) == 0:
            continue
        third = random.choice(list(G.succ[second]))
        if second == third:
            continue

        if G.out_degree(third) == 0:
            continue
        fourth = random.choice(list(G.succ[third]))
        if third == fourth:
            continue

        if (
                third not in G.succ[start]
                and fourth not in G.succ[second]
                and second not in G.succ[third]
        ):
            # Swap nodes
            G.add_edge(start, third)
            G.add_edge(third, second)
            G.add_edge(second, fourth)
            # Move attributes
            for attr in G.edges[(start, second)].keys():
                nx.set_edge_attributes(G, {(start, third): G.edges[(start, second)][attr]}, attr)
                nx.set_edge_attributes(G, {(third, second): G.edges[(third, fourth)][attr]}, attr)
                nx.set_edge_attributes(G, {(second, fourth): G.edges[(second, third)][attr]}, attr)
            G.remove_edge(start, second)
            G.remove_edge(second, third)
            G.remove_edge(third, fourth)
            swapcount += 1
    return G


# Simple progress bar
def progres(current, total, status='', bar_len=50):
    filled_len = int(round(bar_len * current / float(total)))

    percents = round(100.0 * current / float(total), 1)
    bar = '#' * filled_len + ' ' * (bar_len - filled_len)

    fmt = '[%s] %s%s ...%s' % (bar, percents, '%', status)
    print('\b' * len(fmt), end='')  # clears the line
    sys.stdout.write(fmt)
    sys.stdout.flush()
