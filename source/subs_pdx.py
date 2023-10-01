import os
import json
import networkx as nx
import pandas as pd
import functions_sub as fun


ref = pd.read_csv("Reference.csv", low_memory=False)
G = nx.from_pandas_edgelist(ref, source="Interactor 1", target="Interactor 2",
                            edge_attr=["Confidence", "is_directed", "BIOGRID", "pathway", "complex"],
                            create_using=nx.DiGraph())

files = os.listdir("pdx_simulations/")
CLs = []
for file in files:
    if file.startswith("X"):
        CLs.append(file.split("_")[0])
CLs = set(CLs)

for CL in CLs:
    result = fun.find_all_sub_nodes("pdx_simulations/" + CL, pdx=True)
    if result is None:
        continue
    mutated, isolated, gca_nodes = result

    all_nodes = (mutated.union(isolated)).union(gca_nodes)
    with open("random_seeds/random_seeds_specificities.json", "r") as f:
        seed_spec = json.load(f)
    
    if len(all_nodes) != 0:
        with open(f"pdx_modules/{CL}.txt", "w") as f1:
           f1.write("Gene\t" + "seed_specificity" + "\n")
           for gene in all_nodes:
               try:
                   f1.write(gene + "\t" + str(seed_spec[gene]) + "\n")
               except KeyError:
                   f1.write(gene + "\t" + str(0) + "\n")
           f1.close()
        sub_net = G.subgraph(all_nodes)
        nx.write_edgelist(sub_net, f"pdx_networks/{CL}.txt")

