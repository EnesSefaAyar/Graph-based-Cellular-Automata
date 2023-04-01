import functions as f
import pandas as pd
import numpy as np
import networkx as nx
from copy import deepcopy
import json

# parameters
CL = "ACH-000001"
num_iter = 50
mutations_at = 20
exp_thr = 0.585

# Load data
expressions = pd.read_csv(r"CCLE_expression.csv", low_memory=False)
info = pd.read_csv(r"sample_info.csv", low_memory=False)
reference = pd.read_csv(r"Reference.csv", low_memory=False)

# Reference Network Construct
G = nx.from_pandas_edgelist(reference, source="Interactor 1", target="Interactor 2",
                            edge_attr=["Confidence", "is_directed", "BIOGRID", "pathway", "complex"],
                            create_using=nx.DiGraph())

# Remove nodes do not having a state or isolated.
G.remove_nodes_from([node for node in G.nodes if node not in expressions.columns])
G.remove_nodes_from(list(nx.isolates(G)))

# Linage specific node states (mean of lineage samples for each gene)
lineage = info["lineage"][info["DepMap_ID"] == CL].values[0]
cell_lines = list(info["DepMap_ID"][info["lineage"] == lineage])
activity_dict = expressions[expressions["DepMap_ID"].isin(cell_lines)].drop("DepMap_ID", axis=1).mean().to_dict()

# Set initial node states
nx.set_node_attributes(G, activity_dict, "activity")

# Curate significant mutations for cell lines
G.name = CL
mutation_list = f.curate_mutations(G)

################################################################################################################
list_of_graphs = [G, deepcopy(G)]  # Keep track of updated graphs
states = []
node_states = nx.get_node_attributes(G, "activity")
states.append(node_states)  # Ordered list of dictionaries keep state of each node during iterations
mutated_nodes = []  # Ordered list of applied mutations
isolated_nodes = []  # Isolated node (gene) names

# Iterations start
f.progres(0, num_iter, status='Started', bar_len=num_iter)
for i in range(num_iter):
    # GCA starts
    for node in list_of_graphs[-1].nodes:
        estimated_change = 0
        weights = 0
        for edge in list_of_graphs[-2].in_edges(node):
            if list_of_graphs[-2].nodes[edge[0]]["activity"] > exp_thr:
                if list_of_graphs[-2].edges[edge]["is_directed"] == 1:  # If stimulation, estimate is over one std. of mean
                    estimate_parameters = f.estimator(edge[0], edge[1], list_of_graphs[-2].nodes[edge[0]]["activity"], 1)
                    if estimate_parameters is None:
                        continue
                    estimated_change += (estimate_parameters[0] + 2 * estimate_parameters[1]) * list_of_graphs[-2].edges[edge]["Confidence"] * (list_of_graphs[-2].nodes[edge[0]]["activity"] / max(list(expressions[edge[0]])))
                    weights += list_of_graphs[-2].edges[edge]["Confidence"] * (
                            list_of_graphs[-2].nodes[edge[0]]["activity"] / max(list(expressions[edge[0]])))
                elif list_of_graphs[-2].edges[edge]["is_directed"] == -1:  # If repression, estimate is under one std. of mean
                    estimate_parameters = f.estimator(edge[0], edge[1], list_of_graphs[-2].nodes[edge[0]]["activity"], -1)
                    if estimate_parameters is None:
                        continue
                    estimated_change += (estimate_parameters[0] - 2 * estimate_parameters[1]) * list_of_graphs[-2].edges[edge]["Confidence"] * (list_of_graphs[-2].nodes[edge[0]]["activity"] / max(list(expressions[edge[0]])))
                    weights += list_of_graphs[-2].edges[edge]["Confidence"] * (list_of_graphs[-2].nodes[edge[0]]["activity"] / max(list(expressions[edge[0]])))
        if weights > 0:
            list_of_graphs[-1].nodes[node]["activity"] = max((estimated_change / weights + list_of_graphs[-2].nodes[node]["activity"]) / 2, 0)
    list_of_graphs.append(deepcopy(list_of_graphs[-1]))
    # GCA ends
    node_states = nx.get_node_attributes(list_of_graphs[-2], "activity")
    states.append(node_states)
    # Apply mutations
    if i == mutations_at:
        effected_edges = 0
        for mutation in mutation_list:
            effected_edges = f.apply_single_mutation(list_of_graphs[-2], mutation)
            if effected_edges != 0:
                mutated_nodes.append(mutation)
                effected_edges = 0
        isolated_nodes = f.determine_isolated_nodes(list_of_graphs[-2])
    list_of_graphs = list_of_graphs[1:]  # Release unused Graph from memory.
    f.progres(i, num_iter, status='Running', bar_len=num_iter)
f.progres(num_iter, num_iter, status='Iterations Finished', bar_len=num_iter)

# GCA nodes
# Load states in without mutations simulation
with open(f'{lineage}_wo_states.json', 'r') as fin:
    wo_states = json.load(fin)
gca_nodes = set([])
for key in list(states[0].keys()):
    mut = []
    not_mut = []
    for i, iter in enumerate(wo_states):
        not_mut.append(iter[key])
        mut.append(states[i][key])
    difference = np.max(np.abs(np.array(mut) - np.array(not_mut)))
    if difference > exp_thr:
        gca_nodes.add(key)

# List of mutated, isolated and GCA nodes
with open(f'{CL}_mutated_nodes.json', 'w') as fout:
    json.dump(mutated_nodes, fout)
with open(f'{CL}_isolated_nodes.json', 'w') as fout:
    json.dump(isolated_nodes, fout)
with open(f'{CL}_gca_nodes.json', 'w') as fout:
    json.dump(isolated_nodes, fout)

# Calculate the giant component
all_nodes = set(isolated_nodes).union(set(mutated_nodes).union(gca_nodes))
sub = nx.Graph(nx.subgraph(G, all_nodes))
sub.remove_nodes_from(list(nx.isolates(sub)))
giant_sub = nx.DiGraph(nx.subgraph(G, max(nx.connected_components(sub), key=len)))

# Add specificity attributes
with open("random_seeds_specificities.json", "r") as f:
    seed_spec = json.load(f)
with open("random_edge_specificities.json", "r") as f:
    edge_spec = json.load(f)
for node in giant_sub.nodes():
    nx.set_node_attributes(giant_sub, {node: seed_spec[node]}, "seed_spec")
    nx.set_node_attributes(giant_sub, {node: edge_spec[node]}, "edge_spec")

# Output subnetwork in .graphml format
nx.write_graphml(giant_sub, f"{CL}_giant.graphml")
