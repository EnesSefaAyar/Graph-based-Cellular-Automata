import json
import pandas as pd

drugs = pd.read_csv("sorted_drugs.csv", low_memory=False)
ref = pd.read_csv("Reference.csv", low_memory=False)
samples = pd.read_csv("raw_data/sample_info.csv", low_memory=False)
exp_thr = 0.585


def find_all_sub_nodes(CL_name, pdx=False, tcga=False):
    try:
        # Load mutated nodes
        with open(fr"{CL_name}_mutated_nodes.json", "r") as read_file:
            mutated = json.load(read_file)
        # Load without mutations for lineage
        if pdx:
            linage = "esophagus"
        elif tcga:
            linage = "lung"
        else:
            linage = list(samples["lineage"][samples["DepMap_ID"] == CL_name.split("/")[1]])[0]
        with open(fr"wo_states/{linage}_wo_states.json", "r") as read_file:
            without_mutations = json.load(read_file)
        # Load mutated CL states
        with open(fr"{CL_name}_states.json", "r") as read_file:
            with_mutations = json.load(read_file)

        # Load isolated nodes
        with open(fr"{CL_name}_isolated_nodes.json", "r") as read_file:
            isolated = json.load(read_file)
        if isolated is None:
            isolated = {}
        gca_nodes = set()
        for key in list(with_mutations[0].keys()):
            for i, iteration in enumerate(without_mutations):
                if abs(iteration[key] - with_mutations[i][key]) > exp_thr:
                    gca_nodes.add(key)
        mutated_nodes = set([])
        for mut in mutated:
            mutated_nodes.add(mut[1][0])
        return mutated_nodes, set(isolated), gca_nodes
    except (IndexError, FileNotFoundError):
        return None
