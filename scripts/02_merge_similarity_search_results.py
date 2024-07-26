import os
import pandas as pd
from tqdm import tqdm
import numpy as np
import collections

root = os.path.dirname(os.path.abspath(__file__))

df = pd.read_csv(os.path.join(root, "..", "data", "all_molecules.csv"))
df = df[df["inchikey"] != ""]
df = df[df["inchikey"].notnull()]

inchikeys = sorted(set(df["inchikey"].tolist()))
print("Unique molecules in the dataset", len(inchikeys))

ik2names = collections.defaultdict(list)
ik2cat = collections.defaultdict(list)

for r in tqdm(df.values):
    ik = r[1]
    ik2names[ik] += [r[0]]
    ik2cat[ik] += [r[-1]]

ik2name = dict((k, "; ".join(sorted(set(v)))) for k, v in ik2names.items())
ik2cat = dict((k, sorted(set(v))) for k, v in ik2cat.items())

R = []
for r in tqdm(df.values):
    ik = r[1]
    name = ik2name[ik]
    cats = ik2cat[ik]
    if len(cats) == 2:
        cat = "both"
    else:
        cat = cats[0]
    if cat == "both":
        catbin = [1, 1]
    elif cat == "natural":
        catbin = [0, 1]
    elif cat == "synthetic":
        catbin = [1, 0]
    else:
        raise ValueError(cat)
    fn = os.path.join(root, "..", "results", "cheese", f"{ik}.csv")
    if os.path.exists(fn):
        ds = pd.read_csv(fn)
        for s in ds.values:
            R += [[name, cat] + catbin + list(s)]
    else:
        print("Missing", ik, len(ik))

columns = [
    "query_name",
    "query_category",
    "query_is_synthetic",
    "query_is_natural",
    "query_smiles",
    "query_inchikey",
    "smiles",
    "identifier",
    "search_type",
    "score",
    "database",
]
df = pd.DataFrame(R, columns=columns)

df.to_csv(os.path.join(root, "..", "results", "cheese_search.csv"), index=False)

df = pd.read_csv(os.path.join(root, "..", "results", "cheese_search.csv"))

databases = ["zinc15", "enamine-real"]
categories = ["natural", "synthetic", "both"]
search_types = ["consensus", "morgan", "espsim_shape", "espsim_electrostatic"]

R = []
for database in databases:
    for category in categories:
        for search_type in search_types:
            df_ = df[df["database"] == database]
            df_ = df_[df_["query_category"] == category]
            df_ = df_[df_["search_type"] == search_type]
            r = [
                database,
                category,
                search_type,
                df_.shape[0],
                df_["score"].mean(),
                df_["score"].std(),
                df_["score"].median(),
                np.percentile(df_["score"], 25),
                np.percentile(df_["score"], 75),
                df_["score"].min(),
                df_["score"].max(),
            ]
            R.append(r)

dr = pd.DataFrame(
    R,
    columns=[
        "database",
        "category",
        "search_type",
        "counts",
        "mean",
        "std",
        "median",
        "perc_25",
        "perc_75",
        "min",
        "max",
    ],
)
dr.to_csv(
    os.path.join(root, "..", "results", "cheese_search_aggregate.csv"), index=False
)
