import json
import pandas as pd
import numpy as np
import networkx as nx

drugs = pd.read_csv("sorted_drugs.csv", low_memory=False)
ref = pd.read_csv("Reference.csv", low_memory=False)
drugs = drugs.groupby("DRUG_ID")
G = nx.from_pandas_edgelist(ref, source="Interactor 1", target="Interactor 2",
                            edge_attr=["Confidence", "is_directed", "BIOGRID", "pathway", "complex"],
                            create_using=nx.DiGraph())

for sub in drugs:
    z_scores = np.array(sub[1]["Z_SCORE_PUBLISHED"])
    CLs = np.array(sub[1]["DepMap_ID"])
    sensitives = CLs[z_scores < 0]
    resistants = CLs[z_scores > 0]
    sensitive_nodes = []
    for CL in sensitives:
        try:
            df = pd.read_csv("cell_line_modules/" + CL + ".txt", delimiter="\t", low_memory=False)
            sensitive_nodes.append(set(df[df.columns[0]]))
        except FileNotFoundError:
            pass
    resistant_nodes = []
    for CL in resistants:
        try:
            df = pd.read_csv("cell_line_modules/" + CL + ".txt", delimiter="\t", low_memory=False)
            resistant_nodes.append(set(df[df.columns[0]]))
        except FileNotFoundError:
            pass
    if len(sensitive_nodes) == 0:
        continue
    sens_union = sensitive_nodes[0]
    for sens in sensitive_nodes:
        sens_union = sens_union.union(sens)
    if len(resistant_nodes) == 0:
        continue
    res_union = resistant_nodes[0]
    for res in resistant_nodes:
        res_union = res_union.union(res)

    sens_majority = set([])
    for node in sens_union:
        counter = 0
        for sens in sensitive_nodes:
            if node in sens:
                counter += 1
        if counter > len(sensitive_nodes) // 2:
            sens_majority.add(node)

    res_majority = set([])
    for node in res_union:
        counter = 0
        for res in resistant_nodes:
            if node in res:
                counter += 1
        if counter > len(resistant_nodes) // 2:
            res_majority.add(node)

    sens_unique = sens_majority.difference(res_majority)
    res_unique = res_majority.difference(sens_majority)

    with open("random_seeds/random_seeds_specificities.json", "r") as f:
        seed_spec = json.load(f)

    if len(sens_unique) != 0 and len(res_unique) != 0:
        with open(f"drug_modules/{sub[0]}sensitive.txt", "w") as f1:
            f1.write("Gene\t" + "seed_specificity" + "\n")
            for gene in sens_unique:
                if gene in seed_spec:
                    f1.write(gene + "\t" + str(seed_spec[gene]) + "\n")
                else:
                    f1.write(gene + "\t" + str(0) + "\n")
            f1.close()

        sub_net = G.subgraph(sens_unique)
        nx.write_edgelist(sub_net, f"drug_networks/{sub[0]}sensitive.txt")

        with open(f"drug_modules/{sub[0]}resistant.txt", "w") as f2:
            f2.write("Gene\t" + "seed_specificity" + "\n")
            for gene in res_unique:
                if gene in seed_spec:
                    f2.write(gene + "\t" + str(seed_spec[gene]) + "\n")
                else:
                    f2.write(gene + "\t" + str(0) + "\n")
            f2.close()
        sub_net = G.subgraph(res_unique)
        nx.write_edgelist(sub_net, f"drug_networks/{sub[0]}resistant.txt")
