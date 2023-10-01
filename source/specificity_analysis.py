import json
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

ref = pd.read_csv("Reference.csv", low_memory=False)
samples = pd.read_csv("raw_data/sample_info.csv", low_memory=False)

exp_thr = 0.585
CL_name = "ACH-000001"
# list of all nodes and zero count list of each node
list_of_nodes = list(set(ref["Interactor 1"]).union(set(ref["Interactor 2"])))
counts = [0] * len(list_of_nodes)

for i in range(100):
    # Load random simulation data
    linage = list(samples["lineage"][samples["DepMap_ID"] == CL_name])[0]
    with open(fr"wo_states/{linage}_wo_states.json", "r") as read_file:
        without_mutations = json.load(read_file)
    with open(fr"random_seeds/{CL_name}_seed_random{i}_isolated_nodes.json",
              "r") as read_file:  # Isolated nodes in random simulation
        isolated_nodes = json.load(read_file)
    with open(fr"random_seeds/{CL_name}_seed_random{i}_states.json",
              "r") as read_file:  # States of nodes through iterations in random simulation
        states = json.load(read_file)
    # Determine gca nodes for the randomization
    gca_nodes = set([])
    for key in list(states[0].keys()):
        mut = []
        not_mut = []
        for j, iteration in enumerate(without_mutations):
            not_mut.append(iteration[key])
            mut.append(states[j][key])
        difference = np.max(np.abs(np.array(mut) - np.array(not_mut)))
        if difference > exp_thr:
            gca_nodes.add(key)
    # Determine isolated nodes if any
    isolated_nodes = set()
    try:  # In case there is no isolated nodes
        isolated_nodes = set(isolated_nodes)
    except:
        pass
    all_nodes = isolated_nodes.union(gca_nodes)  # Resulting number of nodes (output)
    # Increase the occurance count of each node
    for k, node in enumerate(list_of_nodes):
        if node in all_nodes:
            counts[k] += 1

# Specificities
specificities = dict(zip(list_of_nodes, list(1 - np.array(counts) / 100)))
with open(fr"random_seeds/random_seeds_specificities.json", "w") as read_file:
    json.dump(specificities, read_file)

# Histogram of counts
counts = np.array(counts)
plt.axis([0, 101, 0, 800])
plt.hist(counts[counts > 0], bins=100, range=(0.5, 100 + 0.5))
plt.xlabel("Number of Occurences")
plt.ylabel("Number of Nodes")
plt.title("Random Seeds Occurance Distribution")
plt.savefig(r"random_seeds/specificities_histogram.svg")
