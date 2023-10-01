import json
import os
import pandas as pd
import functions_sub as fun
import random
import numpy as np


random.seed(1)
drugs = pd.read_csv("sorted_drugs.csv", low_memory=False)
CLs = list(set(drugs["DepMap_ID"]))

ratio = len(CLs)//10

for j in range(10):
    R1 = list(random.sample(CLs, ratio))
    drugs[drugs["DepMap_ID"].isin(R1)].to_csv(f"ten_fold_CV/{j}/test.csv", index=False)
    drugs[np.invert(drugs["DepMap_ID"].isin(R1))].to_csv(f"ten_fold_CV/{j}/train.csv", index=False)
    CLs = [CL for CL in CLs if CL not in R1]

for i in range(10):
    drugs = pd.read_csv(f"ten_fold_CV/{i}/train.csv", low_memory=False)
    ref = pd.read_csv("Reference.csv", low_memory=False)
    drugs = drugs.groupby("DRUG_ID")
    for sub in drugs:
        z_scores = np.array(sub[1]["Z_SCORE_PUBLISHED"])
        CLs = np.array(sub[1]["DepMap_ID"])
        sensitives = CLs[z_scores < 0]
        resistants = CLs[z_scores > 0]
        sensitive_nodes = []
        for CL in sensitives:
            try:
                df = pd.read_csv(f"cell_line_modules/{CL}.txt", delimiter="\t", low_memory=False)
            except FileNotFoundError as e:
                continue
            sensitive_nodes.append(set(df[df.columns[0]]))
        resistant_nodes = []
        for CL in resistants:
            try:
                df = pd.read_csv(f"cell_line_modules/{CL}.txt", delimiter="\t", low_memory=False)
            except FileNotFoundError as e:
                continue
            resistant_nodes.append(set(df[df.columns[0]]))
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
            if counter > len(sensitive_nodes)//2:
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
            with open(f"ten_fold_CV/{i}/{sub[0]}sensitive.txt", "w") as f1:
                f1.write("Gene\t" + "seed_specificity" + "\n")
                for gene in sens_unique:
                    try:
                        f1.write(gene + "\t" + str(seed_spec[gene]) + "\n")
                    except KeyError:
                        f1.write(gene + "\t" + str(0) + "\n")
                f1.close()

            with open(f"ten_fold_CV/{i}/{sub[0]}resistant.txt", "w") as f2:
                f2.write("Gene\t" + "seed_specificity" + "\n")
                for gene in res_unique:
                    try:
                        f2.write(gene + "\t" + str(seed_spec[gene]) + "\n")
                    except KeyError:
                        f2.write(gene + "\t" + str(0) + "\n")
                f2.close()

for i in range(10):
    df = pd.read_csv(f"ten_fold_CV/{i}/test.csv")
    CLs = set(df["DepMap_ID"])
    for CL in CLs:
        result = fun.find_all_sub_nodes("cell_line_simulations/" + CL)
        if result is None:
            continue
        mutated, isolated, gca_nodes = result
        all_nodes = (mutated.union(isolated)).union(gca_nodes)

        if len(all_nodes) != 0:
            with open(f"ten_fold_CV/{i}/test_modules/{CL}.txt", "w") as f1:
                f1.write("Gene\t" + "seed_specificity" + "\n")
                for gene in all_nodes:
                    try:
                        f1.write(gene + "\t" + str(seed_spec[gene]) + "\n")
                    except KeyError:
                        f1.write(gene + "\t" + str(0) + "\n")
                f1.close()
