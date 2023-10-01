import json
import pandas as pd
import numpy as np

df = pd.read_csv("raw_data/CCLE_expression.csv", low_memory=False)
ref = pd.read_csv("Reference.csv", low_memory=False)


def estimator(gene1, gene2, target, type=-1, window_size=0.1, exp_thr=0.585, pct=0.1, sensitivity=0.01, min_points=50):
    try:
        x = np.array(df[gene1])
        y = np.array(df[gene2])
        t_index = np.array([False] * len(x))
        t_index += (x > exp_thr) & (y > exp_thr)
        x = x[t_index]
        y = y[t_index]
        o_index = np.array([False] * len(x))
        o_index += (x > (np.mean(x) - np.std(x) * 2)) & (x < (np.mean(x) + np.std(x) * 2)) \
                   & (y > (np.mean(y) - np.std(y) * 2)) & (y < (np.mean(y) + np.std(y) * 2))
        x = x[o_index]
        y = y[o_index]
    except:
        return None
    if np.sum(t_index) > min_points:
        index = np.array([False] * len(x))
        while np.sum(index) < pct * np.sum(t_index):
            index += (x > (target - window_size)) & (x < (target + window_size))
            window_size += sensitivity
        b = y[index]
        if window_size < 2 * exp_thr:
            return np.mean(b), np.std(b)
        else:
            if target > np.mean(df[gene1]):
                if type == 1:
                    return np.mean(b) * (target ** 0.01), np.std(b)
                else:
                    return np.mean(b), np.std(b)
            else:
                if type == 1:
                    return np.min(b), np.std(b)
                else:
                    return np.mean(b), np.std(b)
    else:
        return None


exp_thr = 0.585
estimations_dict = dict()
for j, row in ref.iterrows():
    try:
        if pd.notna(row[2]):
            if row[2] == 1:
                x = np.array(df[row[0]])
                y = np.array(df[row[1]])
                t_index = np.array([False] * len(x))
                t_index += (x > exp_thr) & (y > exp_thr)
                x = x[t_index]
                y = y[t_index]
                inside_dict = {}
                for i in np.arange(0, max(list(df[row[0]])) + 3, 1):
                    if estimator(row[0], row[1], i, 1) is None:
                        continue
                    mean, std = estimator(row[0], row[1], i, 1)
                    if mean + 2 * std > max(y):
                        inside_dict[int(i)] = max(y)
                    else:
                        inside_dict[int(i)] = mean + 2 * std
            elif row[2] == -1:
                try:
                    x = np.array(df[row[0]])
                    y = np.array(df[row[1]])
                    t_index = np.array([False] * len(x))
                    t_index += (x > exp_thr) & (y > exp_thr)
                    x = x[t_index]
                    y = y[t_index]

                    inside_dict = {}
                    for i in np.arange(0, max(list(df[row[0]])) + 3, 1):
                        if estimator(row[0], row[1], i, 1) is None:
                            continue
                        mean, std = estimator(row[0], row[1], i, -1)

                        if mean - 2 * std < min(y):
                            inside_dict[int(i)] = min(y)
                        else:
                            inside_dict[int(i)] = mean - 2 * std
                except:
                    continue
        else:
            continue
    except:
        continue
    if len(inside_dict) != 0:
        estimations_dict[row[0] + "_" + row[1]] = inside_dict
    else:
        pass

with open('estimator.json', 'w') as f:
    json.dump(estimations_dict, f)
