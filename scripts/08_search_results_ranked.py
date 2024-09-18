import os
from scipy.stats import rankdata
import pandas as pd
import collections

root = os.path.dirname(os.path.abspath(__file__))

print("Loading databases for inchikey")

ik2db_ = collections.defaultdict(list)

df = pd.read_csv(os.path.join(root, "..", "results", "cheese_search.csv"))
for r in df.iterrows():
    r = r[1]
    ik2db_[r["inchikey"]].append(r["database"])
ik2db_ = dict((k, list(set(v))) for k,v in ik2db_.items())
ik2db = {}
for k,v in ik2db_.items():
    if len(v) == 1:
        ik2db[k] = v[0]
    else:
        ik2db[k] = "both"

print("Loading manually annotated molecules inchikeys")
ma_iks = set(pd.read_csv(os.path.join(root, "..", "data", "all_molecules_with_cheese_results.csv"))['inchikey'])

print("Loading classifier data")

df = pd.read_csv(os.path.join(root, "..", "results", "classifiers", "classifier_output_macs.csv"))

data = {
    "inchikey": df["inchikey"],
    "smiles": df["smiles"],
    "database": [ik2db[ik] for ik in df["inchikey"]],
    "is_manually_annotated": [1 if ik in ma_iks else 0 for ik in df["inchikey"].tolist()],
    "ma_vs_db": rankdata(df["ma_vs_db"], method="ordinal"),
    "ma_vs_ch0": rankdata(df["ma_vs_ch0"], method="ordinal"),
    "ma_vs_ch1": rankdata(df["ma_vs_ch1"], method="ordinal"),
    "ma_vs_ch2": rankdata(df["ma_vs_ch2"], method="ordinal"),
    "cdc_vs_cdg": rankdata(df["cdc_vs_cdg"], method="ordinal"),
}

df = pd.DataFrame(data)

df["ma_vs_ch"] = (df["ma_vs_ch0"] + df["ma_vs_ch1"] + df["ma_vs_ch2"]) / 3

df["clf_rank"] = (df["ma_vs_db"] + df["ma_vs_ch"] + df["cdc_vs_cdg"]) / 3

df = df[["inchikey", "smiles", "database", "is_manually_annotated", "clf_rank"]]

print("Loading similarity data")

dc = pd.read_csv(os.path.join(root, "..", "results", "cheese_search.csv"))
dc = dc[["inchikey", "score"]]
dc = dc.drop_duplicates(subset="inchikey").reset_index(drop=True)

ik2score_sum = collections.defaultdict(float)
for r in dc.iterrows():
    r = r[1]
    ik2score_sum[r["inchikey"]] += r["score"]

iks = []
scores = []
for k,v in ik2score_sum.items():
    iks.append(k)
    scores.append(v)

ds = pd.DataFrame({"inchikey": iks, "score": scores})
ds["sim_rank"] = rankdata(ds["score"], method="ordinal")

df = df.merge(ds[["inchikey", "sim_rank"]], on="inchikey")
df = df[df["sim_rank"].notnull()]
df["tot_rank"] = (df["clf_rank"] + df["sim_rank"]) / 2

print("Re-scaling ranks")
df["clf_rank"] = rankdata(df["clf_rank"], method="ordinal")/len(df)*100
df["sim_rank"] = rankdata(df["sim_rank"], method="ordinal")/len(df)*100
df["tot_rank"] = rankdata(df["tot_rank"], method="ordinal")/len(df)*100

print("Saving ranked search results")
df = df.sort_values("tot_rank", ascending=False).reset_index(drop=True)

df.to_csv(os.path.join(root, "..", "results", "cheese_search_ranked.csv"), index=False)

print("Removing stereoisomer information from inchikeys")
df["inchikey_flat"] = [ik.split("-")[0] for ik in df["inchikey"]]
iks_flat_to_dedupe = set(df[df["is_manually_annotated"] == 1]["inchikey_flat"])
df_1 = df[(df["inchikey_flat"].isin(iks_flat_to_dedupe)) & (df["is_manually_annotated"] == 1)]
df_0 = df[~df["inchikey_flat"].isin(iks_flat_to_dedupe)]
df = pd.concat([df_1, df_0]).sort_values("tot_rank", ascending=False).reset_index(drop=True)
max_rank_idx = df.groupby('inchikey_flat')['tot_rank'].idxmax()

df = df.loc[max_rank_idx].reset_index(drop=True)
df = df.sort_values("tot_rank", ascending=False).reset_index(drop=True)
df = df.drop(columns=["inchikey_flat"])
print("Re-scaling ranks")
df["clf_rank"] = rankdata(df["clf_rank"], method="ordinal")/len(df)*100
df["sim_rank"] = rankdata(df["sim_rank"], method="ordinal")/len(df)*100
df["tot_rank"] = rankdata(df["tot_rank"], method="ordinal")/len(df)*100

df.to_csv(os.path.join(root, "..", "results", "cheese_search_flat_ranked.csv"), index=False)
