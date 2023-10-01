import numpy as np
import pandas as pd
import os
import json
import networkx as nx
from matplotlib import pyplot as plt

"""
drugidresistant or drugidsensitive is the key.
genes set are the values. ex: ({158resistant: (EGFR, TP53 ...)}
"""
folder = "drug_modules/"
files = os.listdir(folder)
drug_modules = {}
for file in files:
    df = pd.read_csv(folder+file, delimiter="\t", low_memory=False)
    genes = set(df[list(df.columns)[0]])
    key = file.split(".")[0]
    drug_modules[key] = genes

"""
pdx names is the key.
genes set are the values. ex: ({X-1095: (EGFR, TP53 ...)}
"""
folder = "pdx_modules/"
files = os.listdir(folder)
pdx_modules = {}
for file in files:
    df = pd.read_csv(folder+file, delimiter="\t", low_memory=False)
    genes = set(df[list(df.columns)[0]])
    key = file.split(".")[0]
    pdx_modules[key] = genes

"""
Mapping drugs from names to drug ids.
"""
pdx_responses = pd.read_csv("pdx_drugs.csv", low_memory=False)
pdx_responses = pdx_responses[(pdx_responses["MATCH_TYPE"] == "name")]
pdx_responses.reset_index(inplace=True)
pdx_responses.drop("index", axis=1, inplace=True)

PDXs = []
responses = []
for i, row in pdx_responses.iterrows():
    if "PD" in row[10]:
        responses.append(str(int(row[-2])) + "resistant")
        PDXs.append(row[0])
    else:
        responses.append(str(int(row[-2])) + "sensitive")
        PDXs.append(row[0])

"""
pdxs are ordered list of pdxs. [X-1095, ... ]
responses are ordered list of drug responses [1095resistant ...]
"""

# Ready to use network proximity analysis.
with open("pdx_all_proximities.json", "r") as f2:
    pdx_proximities = json.load(f2)

"""
pdx_proximities are dict of dict {pdx: {drugidresistant: (d, z, p) ...}, ...}
"""

accuracies = []
number_of_predictions = []
precisions = []
recalls = []
mccs = []
tp_rate = []
fp_rate = []
z_scores = list(np.arange(-12, -20, -0.25))

correct_drugs = {'1199': [], '1598': [], '1507': [], '1372': [], '1873': []}

for i in z_scores:
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    z_thr = i
    pdx_predicted = set()
    drugs_predicted = set()
    for pdx, all_scores in pdx_proximities.items():
        actual_responses = np.array(responses)[(np.array(PDXs) == pdx)]
        for response_name in actual_responses:
            if response_name in all_scores.keys():
                if "res" in response_name:
                    drug = response_name.split("res")[0]
                else:
                    drug = response_name.split("sens")[0]
                if "res" in response_name:
                    d_r, z_r, p_r = all_scores[response_name]
                    d_s, z_s, p_s = all_scores[response_name.replace("resistant", "sensitive")]
                else:
                    d_s, z_s, p_s = all_scores[response_name]
                    d_r, z_r, p_r = all_scores[response_name.replace("sensitive", "resistant")]
                if d_r < d_s and z_r < i:
                    prediction = "resistant"
                elif d_r > d_s and z_s < i:
                    prediction = "sensitive"
                else:
                    continue

                if prediction in response_name:
                    if prediction == "resistant":
                        tp += 1
                        pdx_predicted.add(pdx)
                        drugs_predicted.add(drug)
                        if i == -14.75:
                            correct_drugs[str(drug)].append(pdx)
                    else:
                        tn += 1
                        pdx_predicted.add(pdx)
                        drugs_predicted.add(drug)
                else:
                    if prediction == "sensitive":
                        fn += 1
                        pdx_predicted.add(pdx)
                        drugs_predicted.add(drug)
                    else:
                        fp += 1
                        pdx_predicted.add(pdx)
                        drugs_predicted.add(drug)
    print("TP:", tp)
    print("FP:", fp)
    print("TN:", tn)
    print("FN:", fn)
    print("Z:", i)
    print("Number of PDXs:", len(pdx_predicted))
    print("Number of Drugs:", drugs_predicted)
    try:
        mccs.append(((tp * tn) - (fp * fn)) / ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5)
    except ZeroDivisionError as e:
        mccs.append(0)
    accuracies.append(100 * (tp + tn) / (fp + tp + fn + tn))
    number_of_predictions.append(tp + fp + tn + fn)
    precisions.append(100 * tp/(tp + fp))
    recalls.append(100 * tp/(tp + fn))
    tp_rate.append(tp / (tp + fn))
    fp_rate.append(fp / (fp + tn))


color = 'tab:blue'
plt.ylabel('Accuracy', color=color)
plt.xlabel("Z-score")
plt.plot(z_scores, accuracies, color=color)
plt.ylim([0, 100])
plt.savefig("PDX_Accuracies.svg")
plt.clf()
plt.gcf()

####################################################

color = 'tab:orange'
plt.ylabel('Precision', color=color)
plt.xlabel("Z-score")
plt.plot(z_scores, precisions, color=color)
plt.ylim([0, 100])
plt.savefig("PDX_Precisions.svg")
plt.gcf()
plt.clf()

####################################################

color = 'tab:green'
plt.ylabel('Recall', color=color)
plt.xlabel("Z-score")
plt.plot(z_scores, recalls, color=color)
plt.ylim([0, 100])
plt.savefig("PDX_Recalls.svg")
plt.gcf()
plt.clf()


####################################################

color = 'tab:purple'
plt.ylabel('MCC', color=color)
plt.xlabel("Z-score")
plt.plot(z_scores, mccs, color=color)
plt.ylim([-1, 1])
plt.savefig("PDX_MCC.svg")
plt.gcf()
plt.clf()

##########################

plt.clf()
plt.title("PDX")
plt.scatter(tp_rate, fp_rate, color='tab:purple')
plt.plot(np.arange(0, 1, 0.1), np.arange(0, 1, 0.1), color="blue", linestyle='dashed')
plt.ylim([0, 1])
plt.xlim([0, 1])
plt.xlabel("FP Rate")
plt.ylabel("TP Rate")
plt.savefig("PDX_ROC.svg")

# Majortiy
pdx_all = set([])
for pdx in correct_drugs["1873"]:
    pdx_mod = set(pd.read_csv(f"pdx_modules/{pdx}.txt", delimiter="\t")["Gene"])
    for gene in pdx_mod:
        pdx_all.add(gene)
pdx_all = list(pdx_all)
pdx_count = [0]*len(pdx_all)
print("Number of PDX samples", len(correct_drugs["1873"]))
for pdx in correct_drugs["1873"]:
    pdx_mod = set(pd.read_csv(f"pdx_modules/{pdx}.txt", delimiter="\t")["Gene"])
    for gene in pdx_mod:
        pdx_count[pdx_all.index(gene)] += 1
pdx_count = np.array(pdx_count)
pdx_all = np.array(pdx_all)
pdx_majority = pdx_all[pdx_count >= len(correct_drugs["1873"])//2]
pdx_majority = set(pdx_majority)
# Add drug targets
pdx_majority.add("MTOR")
pdx_majority.add("PIK3CA")
pdx_majority.add("PIK3CB")
pdx_majority.add("PIK3CD")
pdx_majority.add("PIK3CG")
pdx_majority.add("PIK3C3")

# Drug-PDX module Relations
reference = pd.read_csv(r"Reference.csv", low_memory=False)
G = nx.from_pandas_edgelist(reference, source="Interactor 1", target="Interactor 2",
                            edge_attr=["Confidence", "is_directed", "BIOGRID", "pathway", "complex"],
                            create_using=nx.DiGraph())
# Remove nodes isolated.
G.remove_nodes_from(list(nx.isolates(G)))

# For sens true predicted PDX
sens_drug_module = set(pd.read_csv(f"drug_modules/1873sensitive.txt", delimiter="\t")["Gene"])
res_drug_module = set(pd.read_csv(f"drug_modules/1873resistant.txt", delimiter="\t")["Gene"])

all_nodes = pdx_majority.union(sens_drug_module).union(res_drug_module)

print("All nodes:",len(all_nodes))

sub = G.subgraph(all_nodes)
sub_nodes = [c for c in sorted(nx.weakly_connected_components(sub), key=len, reverse=True)]
sub = G.subgraph(sub_nodes[0])
for edge in sub.edges:
    if edge[0] in sens_drug_module and edge[0] not in res_drug_module:
        if edge[0] in pdx_majority:
            nx.set_node_attributes(sub, {edge[0]: "pdx-sens"}, "type")
        else:
            nx.set_node_attributes(sub, {edge[0]: "sens"}, "type")
    elif edge[0] in res_drug_module and edge[0] not in sens_drug_module:
        if edge[0] in pdx_majority:
            nx.set_node_attributes(sub, {edge[0]: "pdx-res"}, "type")
        else:
            nx.set_node_attributes(sub, {edge[0]: "res"}, "type")
    elif edge[0] in res_drug_module and edge[0] in sens_drug_module:
        if edge[0] in pdx_majority:
            nx.set_node_attributes(sub, {edge[0]: "pdx-res-sens"}, "type")
        else:
            nx.set_node_attributes(sub, {edge[0]: "res-sens"}, "type")
    else:
        nx.set_node_attributes(sub, {edge[0]: "pdx"}, "type")

    if edge[1] in sens_drug_module and edge[1] not in res_drug_module:
        if edge[1] in pdx_majority:
            nx.set_node_attributes(sub, {edge[1]: "pdx-sens"}, "type")
        else:
            nx.set_node_attributes(sub, {edge[1]: "sens"}, "type")
    elif edge[1] in res_drug_module and edge[1] not in sens_drug_module:
        if edge[1] in pdx_majority:
            nx.set_node_attributes(sub, {edge[1]: "pdx-res"}, "type")
        else:
            nx.set_node_attributes(sub, {edge[1]: "res"}, "type")
    elif edge[1] in res_drug_module and edge[1] in sens_drug_module:
        if edge[1] in pdx_majority:
            nx.set_node_attributes(sub, {edge[1]: "pdx-res-sens"}, "type")
        else:
            nx.set_node_attributes(sub, {edge[1]: "res-sens"}, "type")
    else:
        nx.set_node_attributes(sub, {edge[1]: "pdx"}, "type")

print(len(sub.nodes),len(sub.edges))
nx.write_graphml(sub, "res_test.graphml")

# Analyze certain genes
genes_to_analyze = ["SOS2", "FGF1", "PGF", "TFDP2", "VEGFC", "IL3", "GNGT1"]
genes_to_analyze = ["TRAF5"]
with open(f"wo_states/esophagus_wo_states.json") as f:
    wostates = json.load(f)

for gene in genes_to_analyze:
    pdx_w = []
    pdx_wo = []
    for pdx in correct_drugs["1873"]:
        with open(f"pdx_simulations/{pdx}_states.json") as f:
            wstates = json.load(f)
        pdx_mod = set(pd.read_csv(f"pdx_modules/{pdx}.txt", delimiter="\t")["Gene"])
        if gene in pdx_mod:
            with_states = []
            without_states = []
            for i, wstate in enumerate(wstates):
                with_states.append(wstate[gene])
                without_states.append(wostates[i][gene])
            pdx_w.append(with_states)
            pdx_wo.append(without_states)
    plt.clf()
    for j in range(len(pdx_w)):
        plt.plot(range(0, 31), pdx_w[j], color='tab:red', linestyle="dashed")
    plt.xlabel("Iterations")
    plt.ylabel("Activity State")
    plt.title(gene)
    plt.plot(range(0, 31), pdx_wo[j], color='tab:blue')
    plt.show()

