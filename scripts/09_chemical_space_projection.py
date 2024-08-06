import os
import pandas as pd
import numpy as np
from openTSNE import TSNE

root = os.path.dirname(os.path.abspath(__file__))

df_ma = pd.read_csv(os.path.join(root, "..", "data", "all_molecules.csv"))[["inchikey", "smiles"]]
#df_db = pd.read_csv(os.path.join(root, "..", "data", "drugbank_inchikeys.csv"))[["inchikey", "smiles"]]
df_cd = pd.read_csv(os.path.join(root, "..", "data", "chemdiv", "chemdiv_molecules.csv"))[["inchikey", "smiles"]]
df_cs = pd.read_csv(os.path.join(root, "..", "results", "cheese_search_flat_ranked.csv"))[["inchikey", "smiles", "tot_rank"]]

dm_ma = pd.read_csv(os.path.join(root, "..", "data", "all_molecules_eos39co.csv"))
dm_cd = pd.read_csv(os.path.join(root, "..", "data", "chemdiv", "chemdiv_eos39co.csv"))
dm_cs = pd.read_csv(os.path.join(root, "..", "results", "cheese_eos39co.csv"))

df_ma = df_ma[df_ma["inchikey"].isin(dm_ma["key"])].drop_duplicates(subset="inchikey").sort_values(by="inchikey").reset_index(drop=True)
df_cd = df_cd[df_cd["inchikey"].isin(dm_cd["key"])].drop_duplicates(subset="inchikey").sort_values(by="inchikey").reset_index(drop=True)
df_cs = df_cs[df_cs["inchikey"].isin(dm_cs["key"])].drop_duplicates(subset="inchikey").sort_values(by="inchikey").reset_index(drop=True)

dm_ma = dm_ma[dm_ma["key"].isin(df_ma["inchikey"])].drop_duplicates(subset="key").sort_values(by="key").reset_index(drop=True)
dm_cd = dm_cd[dm_cd["key"].isin(df_cd["inchikey"])].drop_duplicates(subset="key").sort_values(by="key").reset_index(drop=True)
dm_cs = dm_cs[dm_cs["key"].isin(df_cs["inchikey"])].drop_duplicates(subset="key").sort_values(by="key").reset_index(drop=True)

x_ma = np.array(dm_ma[list(dm_ma.columns)[2:]])
x_cd = np.array(dm_cd[list(dm_cd.columns)[2:]])
x_cs = np.array(dm_cs[list(dm_cs.columns)[2:]])

print("Projecting chemical space")

x = np.concatenate([x_ma, x_cd, x_cs], axis=0)

inchikeys = list(df_ma["inchikey"]) + list(df_cd["inchikey"]) + list(df_cs["inchikey"])
smiles = list(df_ma["smiles"]) + list(df_cd["smiles"]) + list(df_cs["smiles"])
tot_rank = [-1]*df_ma.shape[0] + [-1]*df_cd.shape[0] + list(df_cs["tot_rank"])
category = ["Manually annotated"]*df_ma.shape[0] + ["ChemDiv Coronavirus"]*df_cd.shape[0] + ["CHEESE search"]*df_cs.shape[0]

data = {
    "inchikey": inchikeys,
    "smiles": smiles,
    "tot_rank": tot_rank,
    "category": category,
}
df = pd.DataFrame(data)
print(df.shape)

p = TSNE().fit(x)

df["x"] = p[:,0]
df["y"] = p[:,1]

df.to_csv(os.path.join(root, "..", "results", "chemical_space_projection.csv"), index=False)