import os
import json
import pandas as pd
import functions_sub as fun
import networkx as nx
from joblib import Parallel, delayed

reference = pd.read_csv(r"Reference.csv", low_memory=False)
G = nx.from_pandas_edgelist(reference, source="Interactor 1", target="Interactor 2",
                            edge_attr=["Confidence", "is_directed", "BIOGRID", "pathway", "complex"],
                            create_using=nx.DiGraph())


def main(CL):
    result = fun.find_all_sub_nodes("cell_line_simulations/" + CL)
    if result is None:
        return
    mutated, isolated, gca_nodes = result
    all_nodes = (mutated.union(isolated)).union(gca_nodes)
    with open("random_seeds/random_seeds_specificities.json", "r") as f:
        seed_spec = json.load(f)

    if len(all_nodes) != 0:
        with open(f"cell_line_modules/{CL}.txt", "w") as f1:
            f1.write("Gene\t" + "seed_specificity" + "\n")
            for gene in all_nodes:
                try:
                    f1.write(gene + "\t" + str(seed_spec[gene]) + "\n")
                except KeyError:
                    f1.write(gene + "\t" + str(0) + "\n")
            f1.close()
        sub = G.subgraph(all_nodes)
        nx.write_edgelist(sub, f"cell_line_networks/{CL}.txt")

files = os.listdir("cell_line_simulations/")
CLs = set()
for file in files:
    if file.startswith("ACH"):
        CLs.add(file.split("_")[0])

res = Parallel(n_jobs=2)(delayed(main)(cell_line) for cell_line in CLs)
