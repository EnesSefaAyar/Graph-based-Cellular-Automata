import pandas as pd
import numpy as np

"""Process mutations"""
mutations = pd.read_csv("raw_data/CCLE_mutations.csv", usecols=["Hugo_Symbol", "Protein_Change", "isDeleterious",
                                                                "isTCGAhotspot", "isCOSMIChotspot",
                                                                "Variant_Classification"],
                        low_memory=False)
# Remove silent mutations and mutations that has no effect on protein level
mutations = mutations[(mutations["Variant_Classification"] != "Silent") & (pd.notna(mutations["Protein_Change"]))]
mutations.reset_index(drop=True, inplace=True)
mutations = mutations[[False if change[2] == change[-1] else True for change in list(mutations["Protein_Change"])]]
mutations.reset_index(drop=True, inplace=True)
mutations.drop("Variant_Classification", axis=1, inplace=True)
effected_aa = []
for i, row in mutations.iterrows():
    if "*" in row[1] or "fs" in row[1] or "del" in row[1]:
        effected_aa.append(-1)
    else:
        try:
            effected_aa.append(int(row[1][3:-1]))
        except:
            effected_aa.append(None)
mutations["effected_aa"] = effected_aa
mutations.dropna(inplace=True)
mutations.reset_index(drop=True, inplace=True)

# Protein interfaces
ref = pd.read_csv("Reference.csv", low_memory=False)
interfaces = pd.read_csv("raw_data/H_sapiens_interfacesHQ.txt", delimiter="\t", usecols=["P1", "P2", "P1_IRES",
                                                                                         "P2_IRES"], low_memory=False)
mapper = pd.read_csv("raw_data/interface_id_map.tsv", delimiter="\t", low_memory=False)
map_dict = dict(zip(mapper["From"], mapper["To"]))
interfaces[interfaces.columns[0]] = interfaces[interfaces.columns[0]].map(map_dict)
interfaces[interfaces.columns[1]] = interfaces[interfaces.columns[1]].map(map_dict)
interfaces.dropna(inplace=True)
interfaces.reset_index(drop=True, inplace=True)
interfaces["P1_IRES"] = interfaces["P1_IRES"].apply(lambda x: x.replace("[", "").replace("]", ""))
interfaces["P2_IRES"] = interfaces["P2_IRES"].apply(lambda x: x.replace("[", "").replace("]", ""))

# Optimization before loop.
interfaces = interfaces[(interfaces["P1_IRES"] != "") | (interfaces["P2_IRES"] != "")]
interfaces["P1_IRES"] = interfaces["P1_IRES"].apply(lambda x: x.split(","))
interfaces["P2_IRES"] = interfaces["P2_IRES"].apply(lambda x: x.split(","))
# Fix the - in between residues 70-74 --> 70, 71, 72, 73, 74
P1_res = list(interfaces["P1_IRES"])
P2_res = list(interfaces["P2_IRES"])
for i in range(len(P1_res)):
    new = P1_res[i].copy()
    for j in range(len(P1_res[i])):
        if "-" in P1_res[i][j]:
            for res in range(int(P1_res[i][j].split("-")[0]), int(P1_res[i][j].split("-")[1]) + 1):
                new.append(str(res))
            new.remove(P1_res[i][j])
    P1_res[i] = set(new)
    new = P2_res[i].copy()
    for j in range(len(P2_res[i])):
        if "-" in P2_res[i][j]:
            for res in range(int(P2_res[i][j].split("-")[0]), int(P2_res[i][j].split("-")[1]) + 1):
                new.append(str(res))
            new.remove(P2_res[i][j])
    P2_res[i] = set(new)
interfaces["P1_IRES"] = P1_res
interfaces["P2_IRES"] = P2_res

effected_interaction = []
l = len(mutations)

for i, row in mutations.iterrows():
    effected_neighbors = ""
    if row[-1] != -1:
        coordinates = list(interfaces[interfaces["P1"] == row[0]]["P1_IRES"])
        negihbors = list(interfaces[interfaces["P1"] == row[0]]["P2"])
        for j, coordinate in enumerate(coordinates):
            if str(int(row[-1])) in coordinate:
                effected_neighbors += negihbors[j] + ";"
        coordinates = list(interfaces[interfaces["P2"] == row[0]]["P2_IRES"])
        negihbors = list(interfaces[interfaces["P2"] == row[0]]["P1"])
        for j, coordinate in enumerate(coordinates):
            if str(int(row[-1])) in coordinate:
                effected_neighbors += negihbors[j] + ";"
                print(100*i/l)
        if len(effected_neighbors) != 0:
            effected_neighbors = effected_neighbors[:-1]
    effected_interaction.append(effected_neighbors)
mutations["effected_interaction"] = effected_interaction
mutations = mutations[(mutations["effected_interaction"] != "") | (mutations["effected_aa"] == -1)]
mutations.to_csv("mutations_test.csv", index=False)

######################################################################################################################
"""Process Pathways, Complexes and Drugs"""

# Pathways (Only KEGG and Wikipathways used)
pathways = pd.read_csv("raw_data/CPDB_pathways_genes.tab", delimiter="\t", usecols=["external_id", "source",
                                                                                    "hgnc_symbol_ids"], low_memory=False)
pathways = pathways[(pathways["source"] == "KEGG") | (pathways["source"] == "Wikipathways")]
pathways_dict = {list(pathways["external_id"])[i]: item.split(",") for i, item in enumerate(list(pathways["hgnc_symbol_ids"]))}

genes_set = set([])
for id, p_genes in pathways_dict.items():
    for gene in p_genes:
        genes_set.add(gene)
ids = []
in_pathways = {}
for gene in genes_set:
    paths = set([])
    for id, p_genes in pathways_dict.items():
        if gene in p_genes:
            paths.add(id)
    in_pathways[gene] = paths

with open("In_pathways.csv", "w") as f:
    f.write("Genes,Pathways\n")
    for key, val in in_pathways.items():
        f.write(key+",")
        line = ""
        for v in val:
            line += v+";"
        f.write(line[:-1]+"\n")

# Complexes (CORUM database used)
complexes = pd.read_csv("raw_data/Complexes_has.csv", low_memory=False)
complexes_dict = {list(complexes["Complex_id"])[i]: set(item.split(";")) for i, item in enumerate(list(complexes["Genes"]))}

genes_set = set([])
for id, c_genes in complexes_dict.items():
    for gene in c_genes:
        genes_set.add(gene)
ids = []
in_complexes = {}
for gene in genes_set:
    comps = set([])
    for id, c_genes in complexes_dict.items():
        if gene in c_genes:
            comps.add(id)
    in_complexes[gene] = comps

with open("In_complexes.csv", "w") as f:
    f.write("Genes,ComplexIDs\n")
    for key, val in in_complexes.items():
        if len(val) > 0:
            f.write(key+",")
            line = ""
            for v in val:
                line += str(v) + ";"
            f.write(line[:-1]+"\n")

# Process drugs-responses
drugs = pd.read_csv("raw_data/sanger-dose-response.csv", usecols=["DATASET", "DRUG_ID", "Z_SCORE_PUBLISHED", "ARXSPAN_ID",
                                                                  "DRUG_NAME"], low_memory=False)
drugs.columns = usecols = ["DATASET", "DRUG_ID", "Z_SCORE_PUBLISHED", "DepMap_ID", "DRUG_NAME"]
drugs = drugs[np.abs(drugs["Z_SCORE_PUBLISHED"]) > 2]
groups = drugs.groupby("DRUG_ID")
new_df = pd.DataFrame(columns=drugs.columns)
for group in groups:
    df = group[1].sort_values("Z_SCORE_PUBLISHED")
    if np.sum(np.array(df["Z_SCORE_PUBLISHED"]) > 0) > 0 and np.sum(np.array(df["Z_SCORE_PUBLISHED"]) < 0) > 0:
        new_df = new_df.merge(df, how="outer")
new_df.dropna(inplace=True)
new_df.to_csv("sorted_drugs.csv", index=False)


patient_drugs = pd.read_csv("raw_data/Xenograft_drug_responses.csv", low_memory=False)
