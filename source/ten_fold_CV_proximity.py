import pandas as pd
import os
import network_proximity as proximity
from joblib import Parallel, delayed
import json


def main(CL, k):
    """
    drugidresistant or drugidsensitive is the key.
    genes set are the values. ex: ({158resistant: (EGFR, TP53 ...)}
    """
    folder = f"ten_fold_CV/{k}/"
    files = os.listdir(folder)
    drug_modules = {}
    for file in files:
        if ".txt" in file:
            df = pd.read_csv(folder + file, delimiter="\t", low_memory=False)
            df = df[df[df.columns[1]] > 0.9]
            genes = set(df[list(df.columns)[0]])
            key = file.split(".")[0]
            drug_modules[key] = genes

    """
    cell line names is the key.
    genes set are the values. ex: ({ACH-00001: (EGFR, TP53 ...)}
    """
    folder = f"ten_fold_CV/{k}/test_modules/"
    files = os.listdir(folder)
    test_modules = {}
    for file in files:
        df = pd.read_csv(folder + file, delimiter="\t", low_memory=False)
        df = df[df[df.columns[1]] > 0.9]
        genes = set(df[list(df.columns)[0]])
        key = file.split(".")[0]
        test_modules[key] = genes

    cell_line_proximities = {}

    c = 0
    l2 = len(drug_modules.keys())

    scores = {}
    genes1 = test_modules[CL]
    for sens_res, genes2 in drug_modules.items():
        c += 1
        # scores[sens_res] = len(genes.intersection(second_genes)) / (len(genes) + len(second_genes))
        if len(genes1) > 0 and len(genes2) > 0:
            d, z, p = proximity.proximity(genes1, genes2, 100, 1)
            scores[sens_res] = (d, z, p)
        else:
            print("Empty gene set!")
            scores[sens_res] = (100, 0, 1)
        print(100*c/l2)
    cell_line_proximities[CL] = scores
    # Save the analysis
    with open(f"ten_fold_CV/{k}/temp_proximity/{CL}_proximities.json", "w") as f2:
        json.dump(cell_line_proximities, f2)
    """
    cell_line_proximities are dict of dict {cell line: {drugidresistant: (d, z, p) ...}, ...}
    """


for i in range(10):
    folder = f"ten_fold_CV/{i}/test_modules/"
    files = os.listdir(folder)
    new_files = set([])
    for file in files:
        new_files.add(file.split(".")[0])

    res = Parallel(n_jobs=-1)(delayed(main)(cell_line, i) for cell_line in new_files)

    folder = f"ten_fold_CV/{i}/temp_proximity/"
    files = os.listdir(folder)
    cell_line_proximities = {}
    for file in files:
        with open(folder+file, "r") as f:
            CL_proximity = json.load(f)
            for key, val in CL_proximity.items():
                cell_line_proximities[key] = val
    with open(f"ten_fold_CV/cell_line_proximities{i}.json", "w") as f:
        json.dump(cell_line_proximities, f)
