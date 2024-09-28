import pandas as pd
import os
from tqdm import tqdm
import collections
import numpy as np
from openTSNE import TSNE

root = os.path.dirname(os.path.abspath(__file__))

demb = pd.read_csv(os.path.join(root, "..", "results", "cheese_eos39co.csv"))
x = np.array(demb[list(demb.columns)[2:]])

p = TSNE().fit(x)
dp = pd.DataFrame(p, columns=["proj_x", "proj_y"])

dh = pd.read_csv(os.path.join(root, "..", "results", "cheese_search_flat_ranked.csv"))
dm = pd.read_csv(os.path.join(root, "..", "data", "all_molecules_with_cheese_results.csv"))
ds = pd.read_csv(os.path.join(root, "..", "results", "cheese_search.csv"))
ds = ds[ds["score"] > 0.5]

iks_counts_05 = collections.defaultdict(list)
for ik in tqdm(dm["inchikey"].tolist()):
    ds_ = ds[ds["query_inchikey"] == ik]
    for ik_ in ds_["inchikey"].tolist():
        iks_counts_05[ik_] += [ik]
iks_counts_05 = {k: len(set(v)) for k, v in iks_counts_05.items()}
print(iks_counts_05)

ds = ds[ds["score"] > 0.7]
iks_counts_07 = collections.defaultdict(list)
for ik in tqdm(dm["inchikey"].tolist()):
    ds_ = ds[ds["query_inchikey"] == ik]
    for ik_ in ds_["inchikey"].tolist():
        iks_counts_07[ik_] += [ik]
iks_counts_07 = {k: len(set(v)) for k, v in iks_counts_07.items()}

R = []
for ik in dh["inchikey"].tolist():
    if ik not in iks_counts_05:
        iks_counts_05[ik] = 0
    if ik not in iks_counts_07:
        iks_counts_07[ik] = 0
    R += [[iks_counts_07[ik], iks_counts_05[ik]]]

df = pd.DataFrame(R, columns=["num_manual_at_07", "num_manual_at_05"])

dh = pd.concat([dh, df, dp], axis=1)

dh.to_csv(os.path.join(root, "..", "results", "cheese_hits_publish.csv"), index=False)


unique_iks = set(dh["inchikey"].tolist())

ds = pd.read_csv(os.path.join(root, "..", "results", "cheese_search.csv"))

ds = ds[ds["score"] > 0.5]
iks_counts_05 = collections.defaultdict(int)
for ik in dm["inchikey"].tolist():
    ds_ = ds[ds["query_inchikey"] == ik]
    for ik_ in ds_["inchikey"].tolist():
        if ik_ in unique_iks:
            iks_counts_05[ik_] += 1

ds = ds[ds["score"] > 0.7]
iks_counts_07 = collections.defaultdict(int)
for ik in dm["inchikey"].tolist():
    ds_ = ds[ds["query_inchikey"] == ik]
    for ik_ in ds_["inchikey"].tolist():
        if ik_ in unique_iks:
            iks_counts_07[ik_] += 1

R = []
for ik in dm["inchikey"].tolist():
    if ik not in iks_counts_05:
        iks_counts_05[ik] = 0
    if ik not in iks_counts_07:
        iks_counts_07[ik] = 0
    R += [[iks_counts_07[ik], iks_counts_05[ik]]]

df = pd.DataFrame(R, columns=["num_manual_at_07", "num_manual_at_05"])

dm = pd.concat([dm, df], axis=1)
dm.to_csv(os.path.join(root, "..", "results", "queries_publish.csv"), index=False)