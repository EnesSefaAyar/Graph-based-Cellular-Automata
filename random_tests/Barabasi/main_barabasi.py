import functions_barabasi as f
import pandas as pd
import networkx as nx
import random
import numpy as np
from copy import deepcopy
import json

# parameters
def main(r, n):
    CL = "ACH-000001"
    num_iter = 50
    mutations_at = 20
    m = 7
    exp_thr = 0.585
    random.seed(r)
    # Random expression dataframe
    expressions = pd.DataFrame(np.random.random(size=(300, n)), columns=list(range(0, n)))
    # Reference
    reference = pd.read_csv(r"Reference.csv", low_memory=False)

    # Create a random Barabasi graph with n nodes.
    G1 = nx.barabasi_albert_graph(n, m, seed=r)
    G1 = nx.convert_node_labels_to_integers(G1, first_label=0)
    G1 = nx.DiGraph(G1)
    G1.remove_edges_from([edge for edge in G1.edges if random.random() < 0.5])
    # Random activities
    activity_dict = dict(zip(list(G1.nodes), [random.randint(0, 10) for i in range(n)]))
    # Set node states
    nx.set_node_attributes(G1, activity_dict, "activity")
    # Random attributes
    biogrid_freq = np.sum(np.array(reference["BIOGRID"]) == True) / len(reference["BIOGRID"])
    complex_freq = np.sum(np.array(reference["complex"]) == True) / len(reference["complex"])
    pathway_freq = np.sum(np.array(reference["pathway"]) == True) / len(reference["pathway"])
    directed_freq = 1 - np.sum(np.array(reference["is_directed"]) == 0) / len(reference["is_directed"])
    positive_freq = np.sum(np.array(reference["is_directed"]) == 1) / len(reference["is_directed"])

    for edge in G1.edges:
        nx.set_edge_attributes(G1, {edge: 0.5 + random.random()}, "Confidence")
        if random.random() < biogrid_freq:
            nx.set_edge_attributes(G1, {edge: True}, "BIOGRID")
        else:
            nx.set_edge_attributes(G1, {edge: False}, "BIOGRID")
        if random.random() < directed_freq:
            if random.random() < positive_freq:
                nx.set_edge_attributes(G1, {edge: 1}, "is_directed")
            else:
                nx.set_edge_attributes(G1, {edge: -1}, "is_directed")
        else:
            nx.set_edge_attributes(G1, {edge: 0}, "is_directed")
    nx.set_edge_attributes(G1, False, "is_effected")
    nx.set_node_attributes(G1, False, "is_mutated")
    # Determine random complexes
    complexes = list(range(0, n//10))
    complex_dict = {}
    for node in G1.nodes:
        for complex in complexes:
            if random.random() < 0.09:
                if node in complex_dict.keys():
                    complex_dict[node].append(random.choice(complexes))
                else:
                    complex_dict[node] = [random.choice(complexes)]
    for edge in G1.edges:
      try:
        if len(set(complex_dict[edge[0]]).intersection(set(complex_dict[edge[1]]))) > 0:
            nx.set_edge_attributes(G1, {edge: True}, "complex")
        else:
            nx.set_edge_attributes(G1, {edge: False}, "complex")
      except:
        nx.set_edge_attributes(G1, {edge: False}, "complex")
    ################################################################################################################
    list_of_graphs = [G1, deepcopy(G1)]
    states = []
    node_states = nx.get_node_attributes(G1, "activity")
    states.append(node_states)
    mutated_nodes = []
    isolated_nodes = []
    for i in range(num_iter):
        # GCA starts
        for node in list_of_graphs[-1].nodes:
            estimated_change = 0
            weights = 0
            for edge in list_of_graphs[-2].in_edges(node):
                if list_of_graphs[-2].nodes[edge[0]]["activity"] > exp_thr:
                    if list_of_graphs[-2].edges[edge]["is_directed"] == 1:  # If stimulation, estimate is over one std. of mean
                        estimate_parameters = f.estimator(edge[0], edge[1], list_of_graphs[-2].nodes[edge[0]]["activity"],
                                                          1, expressions)
                        if estimate_parameters is None:
                            continue
                        estimated_change += (estimate_parameters[0] + 2 * estimate_parameters[1]) * \
                                            list_of_graphs[-2].edges[edge]["Confidence"] * (
                                                        list_of_graphs[-2].nodes[edge[0]]["activity"] / max(
                                                    list(expressions[edge[0]])))
                        weights += list_of_graphs[-2].edges[edge]["Confidence"] * (
                                list_of_graphs[-2].nodes[edge[0]]["activity"] / max(list(expressions[edge[0]])))
                    elif list_of_graphs[-2].edges[edge]["is_directed"] == -1:  # If repression, estimate is under one std. of mean
                        estimate_parameters = f.estimator(edge[0], edge[1], list_of_graphs[-2].nodes[edge[0]]["activity"],
                                                          -1, expressions)
                        if estimate_parameters is None:
                            continue
                        estimated_change += (estimate_parameters[0] - 2 * estimate_parameters[1]) * \
                                            list_of_graphs[-2].edges[edge]["Confidence"] * (
                                                        list_of_graphs[-2].nodes[edge[0]]["activity"] / max(
                                                    list(expressions[edge[0]])))
                        weights += list_of_graphs[-2].edges[edge]["Confidence"] * (
                                    list_of_graphs[-2].nodes[edge[0]]["activity"] / max(list(expressions[edge[0]])))
            if weights > 0:
                list_of_graphs[-1].nodes[node]["activity"] = max(
                    (estimated_change / weights + list_of_graphs[-2].nodes[node]["activity"]) / 2, 0)
        list_of_graphs.append(deepcopy(list_of_graphs[-1]))
        # GCA ends
        node_states = nx.get_node_attributes(list_of_graphs[-2], "activity")
        states.append(node_states)
        # Apply random mutations
        if i == mutations_at:
            mutation_list = [random.choice(list(G1.nodes)) for i in range(20)]
            for mut in mutation_list:
                f.apply_single_mutation(list_of_graphs[-2], mut, complex_dict)
        print(f"Iteration {CL}, {i + 1}")

    # Output
    with open(f'Barabasi_n{n}_r{r}_states.json', 'w') as fout:
        json.dump(states, fout)
    with open(f'Barabasi_n{n}_r{r}_mutated_nodes.json', 'w') as f1:
        json.dump(mutated_nodes, f1)
    with open(f'Barabasi_n{n}_r{r}_isolated_nodes.json', 'w') as f1:
        json.dump(isolated_nodes, f1)