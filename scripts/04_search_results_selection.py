import pandas as pd
import collections
import os
from tqdm import tqdm
from rdkit import Chem
import csv
import numpy as np

root = os.path.dirname(os.path.abspath(__file__))

smiles_and_ik_file = os.path.join(
    root, "..", "results", "cheese_search_smiles_and_inchikey.csv"
)
if not os.path.exists(smiles_and_ik_file):
    df = pd.read_csv(os.path.join(root, "..", "results", "cheese_search.csv"))
    smiles_list = list(set(df["smiles"].tolist()))
    inchikeys = [
        Chem.MolToInchiKey(Chem.MolFromSmiles(smi)) for smi in tqdm(smiles_list)
    ]

    with open(smiles_and_ik_file, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["smiles", "inchikey"])
        for smi, ik in zip(smiles_list, inchikeys):
            if len(ik) != 27:
                continue
            writer.writerow([smi, ik])

smi2ik = {}
for r in pd.read_csv(smiles_and_ik_file).values:
    smi2ik[r[0]] = r[1]

databases = {}
consensus_scores = collections.defaultdict(list)
consensus_syn_scores = collections.defaultdict(list)
consensus_nat_scores = collections.defaultdict(list)
consensus_hits = collections.defaultdict(int)
consensus_syn_hits = collections.defaultdict(int)
consensus_nat_hits = collections.defaultdict(int)

df = pd.read_csv(os.path.join(root, "..", "results", "cheese_search.csv"))
df = df[df["search_type"] == "consensus"][
    ["query_category", "smiles", "score", "database"]
]

for r in tqdm(df.values):
    category = r[0]
    smi = r[1]
    ik = smi2ik[smi]
    score = r[2]
    database = r[3]
    databases[ik] = database
    consensus_scores[ik] += [score]
    consensus_hits[ik] += 1
    if category == "synthetic":
        consensus_syn_scores[ik] += [score]
        consensus_syn_hits[ik] += 1
    elif category == "natural":
        consensus_nat_scores[ik] += [score]
        consensus_nat_hits[ik] += 1
    else:
        continue

smiles = [r[0] for r in pd.read_csv(smiles_and_ik_file).values]

def max_score(scores):
    if len(scores) == 0:
        return 0
    return np.max(scores)

def sum_score(scores):
    if len(scores) == 0:
        return 0
    return np.sum(scores)

R = []
for smi in smiles:
    ik = smi2ik[smi]
    if ik not in consensus_scores:
        continue
    R.append(
        [
            smi,
            ik,
            consensus_hits[ik],
            consensus_syn_hits[ik],
            consensus_nat_hits[ik],
            sum_score(consensus_scores[ik]),
            sum_score(consensus_syn_scores[ik]),
            sum_score(consensus_nat_scores[ik]),
            max_score(consensus_scores[ik]),
            max_score(consensus_syn_scores[ik]),
            max_score(consensus_nat_scores[ik]),
            databases[ik],
        ]
    )

df = pd.DataFrame(
    R,
    columns=[
        "smiles",
        "inchikey",
        "consensus_hits",
        "consensus_syn_hits",
        "consensus_nat_hits",
        "consensus_score",
        "consensus_syn_score",
        "consensus_nat_score",
        "consensus_max_score",
        "consensus_syn_max_score",
        "consensus_nat_max_score",
        "database",
    ],
)

df = df.sort_values(
    by=["consensus_hits", "consensus_score"], ascending=[False, False]
).reset_index(drop=True)
df.to_csv(os.path.join(root, "..", "results", "cheese_popular_hits.csv"), index=False)

# analyze hits per query

df = pd.read_csv(os.path.join(root, "..", "results", "cheese_search.csv"))
print(len(set(df["query_inchikey"])))
print(df.shape)
df = df[df["search_type"] == "consensus"].reset_index(drop=True)

da = pd.read_csv(os.path.join(root, "..", "data", "all_molecules.csv"))
da = da[["inchikey", "smiles", "category"]]
da = da[da["inchikey"].notnull()]
print(da[da["inchikey"].isnull()].shape)
print(len(set(da["inchikey"])))
print(len(set(da["smiles"])))
print(len(smi2ik))
for r in da.values:
    ik = r[0]
    smi = r[1]
    cat = r[2]
    dq = df[df["query_inchikey"] == ik]
    dq = dq.drop_duplicates(inplace=False)
    dq = dq.reset_index(drop=True)
    if dq.shape[0] != 100:
        print(ik, dq)
