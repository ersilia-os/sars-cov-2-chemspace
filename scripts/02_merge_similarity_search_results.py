import os
import pandas as pd
from tqdm import tqdm
import numpy as np
import collections
from rdkit import Chem

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
        inchikeys = []
        for s in ds["smiles"].tolist():
            try:
                inchikeys += [Chem.MolToInchiKey(Chem.MolFromSmiles(s))]
            except:
                inchikeys += [None]
        ds["inchikey"] = inchikeys
        ds = ds[ds["inchikey"].notnull()]
        ds = ds[
            [
                "query_smiles",
                "query_inchikey",
                "smiles",
                "inchikey",
                "identifier",
                "search_type",
                "score",
                "database",
            ]
        ]
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
    "inchikey",
    "identifier",
    "search_type",
    "score",
    "database",
]
df = pd.DataFrame(R, columns=columns)
df = df.drop_duplicates(inplace=False).reset_index(drop=True)

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

# Drugbank results

db_smiles = pd.read_csv(os.path.join(root, "..", "data", "drugbank_smiles.csv"))[
    "Smiles"
].tolist()
db_inchikeys = []
for smiles in db_smiles:
    try:
        inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
    except:
        inchikey = ""
    if len(inchikey) != 27:
        print("Invalid inchikey for", smiles)
        inchikey = None
    db_inchikeys += [inchikey]

db = pd.DataFrame({"smiles": db_smiles, "inchikey": db_inchikeys})
db = db[db["inchikey"].notnull()]
db.to_csv(os.path.join(root, "..", "data", "drugbank_inchikeys.csv"), index=False)

key2smi = {}
key2consenus = collections.defaultdict(list)
key2morgan = collections.defaultdict(list)
key2esp_shape = collections.defaultdict(list)
key2esp_elec = collections.defaultdict(list)

drugbank_counts = 0
for ik in tqdm(db["inchikey"].tolist()):
    fn = os.path.join(root, "..", "results", "cheese", f"{ik}.csv")
    if not os.path.exists(fn):
        continue
    df = pd.read_csv(fn)
    for r in df.iterrows():
        r = r[1]
        smi = r["smiles"]
        ik = Chem.MolToInchiKey(Chem.MolFromSmiles(smi))
        src = r["database"]
        search_type = r["search_type"]
        key = (ik, src)
        key2smi[key] = smi
        if search_type == "consensus":
            key2consenus[key] += [r["score"]]
        elif search_type == "morgan":
            key2morgan[key] += [r["score"]]
        elif search_type == "espsim_shape":
            key2esp_shape[key] += [r["score"]]
        elif search_type == "espsim_electrostatic":
            key2esp_elec[key] += [r["score"]]
        else:
            raise ValueError(search_type)
    drugbank_counts += 1

keys = list(key2smi.keys())
R = []
for key in keys:
    ik, src = key
    smi = key2smi[key]
    for search_type in search_types:
        r = [ik, smi, src, search_type]
        if search_type == "consensus":
            scores = key2consenus[key]
        elif search_type == "morgan":
            scores = key2morgan[key]
        elif search_type == "espsim_shape":
            scores = key2esp_shape[key]
        elif search_type == "espsim_electrostatic":
            scores = key2esp_elec[key]
        if not scores:
            continue
        r += [
            len(scores),
            len(scores) / drugbank_counts,
            np.sum(scores),
            np.sum(scores) / drugbank_counts,
            np.mean(scores),
            np.std(scores),
            np.median(scores),
            np.percentile(scores, 25),
            np.percentile(scores, 75),
            np.min(scores),
            np.max(scores),
        ]
        R.append(r)
ds = pd.DataFrame(
    R,
    columns=[
        "inchikey",
        "smiles",
        "database",
        "search_type",
        "counts",
        "fraction",
        "sum",
        "sum_fraction",
        "mean",
        "std",
        "median",
        "perc_25",
        "perc_75",
        "min",
        "max",
    ],
)
ds = ds.drop_duplicates(inplace=False).reset_index(drop=True)
ds.to_csv(
    os.path.join(root, "..", "results", "cheese_drugbank_search.csv"), index=False
)
