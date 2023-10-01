import numpy as np
import pandas as pd

# iref interactome (downloaded from the omics integrator github page)
iref = pd.read_csv("raw_data/iref_mitab_miscore_2013_08_12_interactome.txt", delimiter="\t",
                   names=["Interactor 1", "Interactor 2", "iref_confidence"], header=None, low_memory=False)

# string interactome (downloaded only human interactome and names mapped, low confidence filtered out)
string = pd.read_csv("raw_data/9606.protein.links.v11.5.txt", delimiter=" ", low_memory=False)
string_info = pd.read_csv("raw_data/9606.protein.info.v11.5.txt", delimiter="\t", low_memory=False)
string_mapper = dict(zip(list(string_info["#string_protein_id"]), list(string_info["preferred_name"])))
string[string.columns[0]] = string[string.columns[0]].map(string_mapper)
string[string.columns[1]] = string[string.columns[1]].map(string_mapper)
string.columns = ["Interactor 1", "Interactor 2", "string_confidence"]
string = string[string["string_confidence"] > 700]
string["string_confidence"] = np.array(string["string_confidence"]) / 1000
string.dropna(inplace=True)
string.drop_duplicates(["Interactor 1", "Interactor 2"], inplace=True)
string.reset_index(drop=True, inplace=True)

# Intact interactome (all human interactions used.)
intact = pd.read_csv("raw_data/intact.txt", usecols=["Alias(es) interactor A", "Alias(es) interactor B",
                                                     "Taxid interactor A", "Taxid interactor B",
                                                     "Confidence value(s)"], delimiter="\t",  low_memory=False)
intact = intact[intact["Taxid interactor A"] == "taxid:9606(human)|taxid:9606(Homo sapiens)"]
intact = intact[intact["Taxid interactor B"] == "taxid:9606(human)|taxid:9606(Homo sapiens)"]
intact.drop(["Taxid interactor A", "Taxid interactor B"], axis=1, inplace=True)
intact.columns = ["Interactor 1", "Interactor 2", "intact_confidence"]
intact["intact_confidence"] = np.array(intact["intact_confidence"].apply(lambda x: x[-4:]), dtype=float)
intact = intact[intact["Interactor 1"].str.contains("(gene name)")]
intact = intact[intact["Interactor 2"].str.contains("(gene name)")]
intact.reset_index(drop=True, inplace=True)
intact["Interactor 1"] = intact["Interactor 1"].apply(lambda x: x.split("(gene name)")[0].split("uniprotkb:")[-1])
intact["Interactor 2"] = intact["Interactor 2"].apply(lambda x: x.split("(gene name)")[0].split("uniprotkb:")[-1])
intact.dropna(inplace=True)
intact.drop_duplicates(["Interactor 1", "Interactor 2"], inplace=True)
intact.reset_index(drop=True, inplace=True)

# Hippie interactome (All human interactions used.)
hippie = pd.read_csv("raw_data/hippie_current.txt", delimiter="\t", header=None, usecols=[0, 2, 4],
                     names=["Interactor 1", "Interactor 2", "iref_confidence"], low_memory=False)
hippie.columns = ["Interactor 1", "Interactor 2", "hippie_confidence"]
hippie.dropna(inplace=True)
hippie["Interactor 1"] = hippie["Interactor 1"].apply(lambda x: x.split("_")[0])
hippie["Interactor 2"] = hippie["Interactor 2"].apply(lambda x: x.split("_")[0])
hippie.dropna(inplace=True)
hippie.drop_duplicates(["Interactor 1", "Interactor 2"], inplace=True)
hippie.reset_index(drop=True, inplace=True)

# OmniPath (All interactions are used.)
omnipath = pd.read_csv("raw_data/AllInteractions.csv", usecols=["source", "target", "is_directed", "is_stimulation",
                                                                "is_inhibition", "consensus_direction",
                                                                "consensus_stimulation", "consensus_inhibition"],
                       low_memory=False)
omnipath_info = pd.read_csv("raw_data/idmapping_2023_07_11.tsv", delimiter="\t", low_memory=False)
omnipath_mapper = dict(zip(list(omnipath_info["From"]), list(omnipath_info["To"])))
omnipath[omnipath.columns[0]] = omnipath[omnipath.columns[0]].map(omnipath_mapper)
omnipath[omnipath.columns[1]] = omnipath[omnipath.columns[1]].map(omnipath_mapper)
omnipath.dropna(inplace=True)
omnipath.drop_duplicates(["source", "target"], inplace=True)
omnipath.reset_index(drop=True, inplace=True)
omni_directed = np.zeros(len(omnipath))
for i, row in omnipath.iterrows():
    if row[7] and not row[6]:  # consensus inhibition.
        omni_directed[i] = int(-1)
    elif row[6] and not row[7]:  # consensus stimulation.
        omni_directed[i] = int(1)
    elif row[6] and row[7]:  # consensus conflict.
        if row[3] and not row[4]:
            omni_directed[i] = int(1)  # stimulation
        elif row[4] and not row[3]:
            omni_directed[i] = int(-1)  # repression

omnipath["omni_directed"] = omni_directed
omnipath = omnipath[omnipath["omni_directed"] != 0]
omnipath.drop(["is_directed", "is_stimulation", "is_inhibition", "consensus_direction", "consensus_stimulation",
               "consensus_inhibition"], axis=1, inplace=True)
omnipath.columns = ["Interactor 1", "Interactor 2", "omni_directed"]

# TRRUST (All directed interactions used.)
trrust = pd.read_csv("raw_data/trrust_rawdata.human.tsv", delimiter="\t", header=None, usecols=[0, 1, 2],
                     low_memory=False)
trrust = trrust[trrust[2] != "Unknown"]
trrust.columns = ["Interactor 1", "Interactor 2", "tr_directed"]
trrust.dropna(inplace=True)
trrust.drop_duplicates(["Interactor 1", "Interactor 2"], inplace=True)
trrust.reset_index(drop=True, inplace=True)
trrust[trrust.columns[2]] = trrust[trrust.columns[2]].map({"Repression": -1, "Activation": 1})

# Biogrid (Human interactions limited to (direct interaction), (physical association), (association).)
biogrid = pd.read_csv("raw_data/BIOGRID-ALL-4.4.223.mitab.txt", delimiter="\t", usecols=["Alt IDs Interactor A",
                                                                                         "Alt IDs Interactor B",
                                                                                         "Taxid Interactor A",
                                                                                         "Taxid Interactor B",
                                                                                         "Interaction Types"], low_memory=False)

biogrid = biogrid[biogrid["Taxid Interactor A"] == "taxid:9606"]
biogrid = biogrid[biogrid["Taxid Interactor B"] == "taxid:9606"]
biogrid = biogrid[(biogrid["Interaction Types"] == 'psi-mi:"MI:0407"(direct interaction)') |
                  (biogrid["Interaction Types"] == 'psi-mi:"MI:0915"(physical association)') |
                  (biogrid["Interaction Types"] == 'psi-mi:"MI:0914"(association)')]

biogrid["Alt IDs Interactor B"] = biogrid["Alt IDs Interactor B"].apply(lambda x: x.split("locuslink:")[1].split("|")[0])
biogrid["Alt IDs Interactor A"] = biogrid["Alt IDs Interactor A"].apply(lambda x: x.split("locuslink:")[1].split("|")[0])
biogrid.drop(["Taxid Interactor A", "Taxid Interactor B", "Interaction Types"], axis=1, inplace=True)
biogrid.dropna(inplace=True)
biogrid.reset_index(drop=True, inplace=True)
biogrid.columns = ["Interactor 1", "Interactor 2"]
biogrid["BIOGRID"] = np.array([True]*len(biogrid))

# Pathways (Only KEGG and Wikipathways used)
pathways = pd.read_csv("raw_data/CPDB_pathways_genes.tab", delimiter="\t", usecols=["external_id", "source",
                                                                                    "hgnc_symbol_ids"], low_memory=False)
pathways = pathways[(pathways["source"] == "KEGG") | (pathways["source"] == "Wikipathways")]
pathways_dict = {list(pathways["external_id"])[i]: item.split(",") for i, item in enumerate(list(pathways["hgnc_symbol_ids"]))}

# Complexes (CORUM database used)
complexes = pd.read_csv("raw_data/Complexes_has.csv", low_memory=False)
complexes_dict = {list(complexes["Complex_id"])[i]: set(item.split(";")) for i, item in enumerate(list(complexes["Genes"]))}

# Merge All (All merged into one dataframe.)
all_merged = iref.merge(string, how="outer")
all_merged = all_merged.merge(hippie, how="outer")
all_merged = all_merged.merge(intact, how="outer")
all_merged = all_merged.merge(omnipath, how="outer")
all_merged = all_merged.merge(trrust, how="outer")
all_merged = all_merged.merge(biogrid, how="outer")

all_merged.drop_duplicates(inplace=True)
all_merged.reset_index(drop=True, inplace=True)

# Pathway & Complex Annotation
all_merged["pathway"] = np.array([False]*len(all_merged))
all_merged["complex"] = np.array([False]*len(all_merged))

for genes in pathways_dict.values():
    all_merged.loc[(all_merged["Interactor 1"].isin(genes)) & (all_merged["Interactor 2"].isin(genes)), "pathway"] = True

for genes in complexes_dict.values():
    all_merged.loc[(all_merged["Interactor 1"].isin(genes)) & (all_merged["Interactor 2"].isin(genes)), "complex"] = True

all_merged["BIOGRID"] = all_merged["BIOGRID"].fillna(False)
all_merged = all_merged.fillna(0)

# Advanced filtering: all directed kept, all transition rules applicable hippi high conf. kept)
all_merged = all_merged[(all_merged["omni_directed"] != 0) | (all_merged["tr_directed"] != 0) |
                        (((all_merged["pathway"] == True) | (all_merged["complex"] == True))
                         & (all_merged["BIOGRID"] == True) & (all_merged["hippie_confidence"] > 0.83))]

# Confidence values determined.
all_merged["Confidence"] = np.zeros(len(all_merged))
for i, row in all_merged.iterrows():
    if row[6] != 0 or row[7] != 0:
        if row[4] != 0:
            all_merged["Confidence"].at[i] = row[4]
        elif row[2] != 0 or row[3] != 0 or row[5] != 0:
            all_merged["Confidence"].at[i] = max(row[2], row[3], row[5])
        else:
            all_merged["Confidence"].at[i] = 1
    elif row[4] != 0:
        all_merged["Confidence"].at[i] = row[4]
    else:
        all_merged["Confidence"].at[i] = max(row[2], row[3], row[5])

# OmniPath prioritized in conflicting directionality. 
all_merged["is_directed"] = list(all_merged["omni_directed"])
all_merged["is_directed"] = [list(all_merged["tr_directed"])[i] if element == 0 else element for i, element in enumerate(list(all_merged["omni_directed"]))]

# Extra columns dropped.
all_merged.drop(["iref_confidence",	"string_confidence", "hippie_confidence", "intact_confidence", "omni_directed",
                 "tr_directed"], axis=1, inplace=True)


# Remove self-edges and add the bidirected edges for undirected
all_merged = all_merged[np.array(all_merged["Interactor 1"]) != np.array(all_merged["Interactor 2"])]
intact.reset_index(drop=True, inplace=True)
edges = set(zip(list(all_merged["Interactor 1"]), list(all_merged["Interactor 2"])))
new_df = pd.DataFrame(columns=all_merged.columns)
for i, row in all_merged.iterrows():
    if row[-1] == 0:
        if (row[1], row[0]) not in edges:
            new_row = row.copy()
            new_row["Interactor 1"] = row[1]
            new_row["Interactor 2"] = row[0]
            new_df = pd.concat([new_df, new_row.to_frame().T], ignore_index=True)
all_merged = all_merged.merge(new_df, how="outer")
all_merged.to_csv("Reference.csv", index=False)
