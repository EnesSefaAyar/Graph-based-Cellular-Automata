import pandas as pd
import os
import network_proximity as proximity
from joblib import Parallel, delayed
import json


def main(pdx):
    """
    drugidresistant or drugidsensitive is the key.
    genes set are the values. ex: ({158resistant: (EGFR, TP53 ...)}
    """
    folder = "drug_modules/"
    files = os.listdir(folder)
    drug_modules = {}
    for file in files:
        df = pd.read_csv(folder + file, delimiter="\t", low_memory=False)
        df = df[df[df.columns[1]] > 0.9]
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
        df = pd.read_csv(folder + file, delimiter="\t", low_memory=False)
        df = df[df[df.columns[1]] > 0.9]
        genes = set(df[list(df.columns)[0]])
        key = file.split(".")[0]
        pdx_modules[key] = genes

    pdx_proximities = {}
    c = 0
    l2 = len(drug_modules.keys())

    scores = {}
    genes1 = pdx_modules[pdx]
    for sens_res, genes2 in drug_modules.items():
        c += 1
        if len(genes1) > 0 and len(genes2) > 0:
            d, z, p = proximity.proximity(genes1, genes2, 100, 1)
            scores[sens_res] = (d, z, p)
        else:
            print("Empty gene set!")
            scores[sens_res] = (100, 0, 1)
        print(100*c/l2)
    pdx_proximities[pdx] = scores
    # Save the analysis
    with open(f"pdx_temp_proximity/{pdx}_proximities.json", "w") as f2:
        json.dump(pdx_proximities, f2)
    """
    pdx_proximities are dict of dict {pdx: {drugidresistant: (d, z, p) ...}, ...}
    """


folder = "pdx_modules/"
files = os.listdir(folder)
new_files = set([])
for file in files:
    new_files.add(file.split(".")[0])

res = Parallel(n_jobs=-1)(delayed(main)(cell_line) for cell_line in new_files)

folder = "pdx_temp_proximity/"
files = os.listdir(folder)
pdxs_proximities = {}
for file in files:
    with open(folder+file, "r") as f:
        pdx_proximity = json.load(f)
        for key, val in pdx_proximity.items():
            pdxs_proximities[key] = val
with open("pdx_all_proximities.json", "w") as f:
    json.dump(pdxs_proximities, f)

