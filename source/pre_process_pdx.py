import pandas as pd

pdx_drugs = pd.read_csv("raw_data/Xenograft_drug_responses.csv", low_memory=False)
pdx_drugs = pdx_drugs[(pdx_drugs["Treatment type"] == "single") & (pdx_drugs["Treatment"] != "untreated")]

GDSC = pd.read_csv("raw_data/GDSC2_fitted_dose_response.csv", usecols=["DRUG_ID", "DRUG_NAME", "PUTATIVE_TARGET"], low_memory=False)

GDSC = GDSC[(pd.notna(GDSC["DRUG_NAME"])) & (pd.notna(GDSC["DRUG_ID"])) & (pd.notna(GDSC["PUTATIVE_TARGET"]))]

pdx_drug_names = list(pdx_drugs["Treatment"])
pdx_drug_targets = list(pdx_drugs["Treatment target"])

GDSC["DRUG_NAME"] = GDSC["DRUG_NAME"].apply(lambda x: x.lower())
GDSC["PUTATIVE_TARGET"] = GDSC["PUTATIVE_TARGET"].apply(lambda x: x.lower())
GDSC_drug_ids = list(GDSC["DRUG_ID"])
GDSC_drug_names = list(GDSC["DRUG_NAME"])
GDSC_drug_targets = list(GDSC["PUTATIVE_TARGET"])

drug_mapper = pd.read_csv("raw_data/pdx_drug_mapping.csv", low_memory=False)
drug_map_dict = dict(zip(list(drug_mapper["pdx_drug_names"]), list(drug_mapper["alternatives"])))


pdx_drug_ids = []
match_type = []
for i, name in enumerate(pdx_drug_names):
    alternatives = set()
    try:
        alternatives = drug_map_dict[name.lower()].split(";")
    except:
        pass
    if name.lower() in GDSC_drug_names:
        pdx_drug_ids.append(GDSC_drug_ids[GDSC_drug_names.index(name)])
        match_type.append("name")
    elif len(set(alternatives).intersection(set(GDSC_drug_names))) > 0:
        pdx_drug_ids.append(GDSC_drug_ids[GDSC_drug_names.index(list(set(alternatives).intersection(set(GDSC_drug_names)))[0])])
        match_type.append("name")
    else:
        flag = True
        for target in pdx_drug_targets[i].lower().split(","):
            if "/" in target:
                target = target.split("/")[0]
            for j, targets in enumerate(GDSC_drug_targets):
                if target.lower() in targets.lower() or targets.lower() in target.lower():
                    pdx_drug_ids.append(GDSC_drug_ids[j])
                    match_type.append("target")
                    flag = False
                    break
            if not flag:
                break
        if flag:
            pdx_drug_ids.append(None)
            match_type.append(None)
pdx_drugs["DRUG_ID"] = pdx_drug_ids
pdx_drugs["MATCH_TYPE"] = match_type
pdx_drugs.to_csv("pdx_drugs.csv", index=False)
