import functions as f
import pandas as pd
import networkx as nx
from copy import deepcopy
import json
import random


def main(cell_line, r):
    # parameters
    CL = cell_line
    random.seed(r)
    num_iter = 50
    mutations_at = 20
    exp_thr = 0.585
    # Reference Network Construct
    expressions = pd.read_csv(r"CCLE_expression.csv", low_memory=False)
    info = pd.read_csv(r"sample_info.csv", low_memory=False)
    reference = pd.read_csv(r"Reference.csv", low_memory=False)
    G = nx.from_pandas_edgelist(reference, source="Interactor 1", target="Interactor 2",
                                edge_attr=["Confidence", "is_directed", "BIOGRID", "pathway", "complex"],
                                create_using=nx.DiGraph())

    # Remove nodes do not having a state or isolated.
    G.remove_nodes_from([node for node in G.nodes if node not in expressions.columns])
    G.remove_nodes_from(list(nx.isolates(G)))
    # Linage specific node states (mean of lineage samples for each gene)
    lineage = info["lineage"][info["DepMap_ID"] == CL].values[0]
    cell_lines = list(info["DepMap_ID"][info["lineage"] == lineage])
    activity_dict = expressions[expressions["DepMap_ID"].isin(cell_lines)].drop("DepMap_ID",
                                                                                axis=1).mean().to_dict()
    # Set node states
    nx.set_node_attributes(G, activity_dict, "activity")

    # Curate significant mutations for cell lines
    G.name = CL
    # List of mutations
    mutation_list = f.curate_mutations(G)
    ################################################################################################################
    list_of_graphs = [G, deepcopy(G)]
    states = []
    node_states = nx.get_node_attributes(G, "activity")
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
                    if list_of_graphs[-2].edges[edge][
                        "is_directed"] == 1:  # If stimulation, estimate is over one std. of mean
                        estimate_parameters = f.estimator(edge[0], edge[1],
                                                          list_of_graphs[-2].nodes[edge[0]]["activity"], 1)
                        if estimate_parameters is None:
                            continue
                        estimated_change += (estimate_parameters[0] + 2 * estimate_parameters[1]) * \
                                            list_of_graphs[-2].edges[edge]["Confidence"] * (
                                                        list_of_graphs[-2].nodes[edge[0]]["activity"] / max(
                                                    list(expressions[edge[0]])))
                        weights += list_of_graphs[-2].edges[edge]["Confidence"] * (
                                list_of_graphs[-2].nodes[edge[0]]["activity"] / max(list(expressions[edge[0]])))
                    elif list_of_graphs[-2].edges[edge][
                        "is_directed"] == -1:  # If repression, estimate is under one std. of mean
                        estimate_parameters = f.estimator(edge[0], edge[1],
                                                          list_of_graphs[-2].nodes[edge[0]]["activity"], -1)
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
        # Apply mutations
        if i == mutations_at:
            for mut in mutation_list:
                effected_edges = f.apply_single_mutation(list_of_graphs[-2], mut, randomize=True)
                if effected_edges != 0:
                    mutated_nodes.append(mut)
            isolated_nodes = f.determine_isolated_nodes(list_of_graphs[-2])
        print(f"Iteration {CL}, {i + 1}")

    # Output
    with open(f'{cell_line}_seed_random{r}_states.json', 'w') as fout:
        json.dump(states, fout)
    with open(f'{cell_line}_seed_random{r}_isolated_nodes.json', 'w') as fout:
        json.dump(isolated_nodes, fout)
    with open(f'{cell_line}_seed_random{r}_mutated_nodes.json', 'w') as fout:
        json.dump(mutated_nodes, fout)
