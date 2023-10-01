import functions_gca as f
import pandas as pd
import networkx as nx
from copy import deepcopy
from joblib import Parallel, delayed
import json


def main(cell_line):
    # parameters
    CL = cell_line
    num_iter = 30
    exp_thr = 0.585
    # Reference Network Construct
    expressions = pd.read_csv(r"raw_data/CCLE_expression.csv", low_memory=False)
    info = pd.read_csv(r"raw_data/sample_info.csv", low_memory=False)
    reference = pd.read_csv(r"Reference.csv", low_memory=False)
    G = nx.from_pandas_edgelist(reference, source="Interactor 1", target="Interactor 2",
                                edge_attr=["Confidence", "BIOGRID", "is_directed", "pathway", "complex"],
                                create_using=nx.DiGraph())

    # Remove nodes do not having a state or isolated.
    G.remove_nodes_from([node for node in G.nodes if node not in expressions.columns])
    G.remove_nodes_from(list(nx.isolates(G)))

    # Linage specific node states (mean of lineage samples for each gene)
    lineage = info["lineage"][info["DepMap_ID"] == CL].values[0]
    cell_lines = list(info["DepMap_ID"][info["lineage"] == lineage])
    activity_dict = expressions[expressions["DepMap_ID"].isin(cell_lines)].drop("DepMap_ID", axis=1).mean().to_dict()
    # Set node states
    nx.set_node_attributes(G, activity_dict, "activity")
    G.name = CL
    ################################################################################################################
    list_of_graphs = [G, deepcopy(G)]
    states = []
    node_states = nx.get_node_attributes(G, "activity")
    states.append(node_states)
    for i in range(num_iter):
        # GCA starts
        for node in list_of_graphs[-1].nodes:
            estimated_change = 0
            weights = 0
            for edge in list_of_graphs[-2].in_edges(node):
                if list_of_graphs[-2].nodes[edge[0]]["activity"] > exp_thr:
                    if list_of_graphs[-2].edges[edge]["is_directed"] != 0:
                        estimate = f.estimator(edge[0], edge[1], list_of_graphs[-2].nodes[edge[0]]["activity"])
                        if estimate is None:
                            continue
                        estimated_change += estimate * list_of_graphs[-2].edges[edge]["Confidence"] * (
                                list_of_graphs[-2].nodes[edge[0]]["activity"] / max(list(expressions[edge[0]])))
                        weights += list_of_graphs[-2].edges[edge]["Confidence"] * (
                                list_of_graphs[-2].nodes[edge[0]]["activity"] / max(list(expressions[edge[0]])))
            if weights > 0:
                list_of_graphs[-1].nodes[node]["activity"] = max(
                    (estimated_change / weights + list_of_graphs[-2].nodes[node]["activity"]) / 2, 0)
        list_of_graphs.append(deepcopy(list_of_graphs[-1]))
        # GCA ends
        node_states = nx.get_node_attributes(list_of_graphs[-2], "activity")
        states.append(node_states)
        print(f"Iteration {CL}, {i + 1}")
    # Output
    with open(f'wo_states/{lineage}_wo_states.json', 'w') as fout:
        json.dump(states, fout)


samples = pd.read_csv("raw_data/sample_info.csv", low_memory=False)
cell_lines = set([])
linages = set([])
for i, row in samples.iterrows():
    if row[17] not in linages:
        linages.add(row[17])
        cell_lines.add(row[0])
res = Parallel(n_jobs=-1)(delayed(main)(cell_line) for cell_line in cell_lines)
