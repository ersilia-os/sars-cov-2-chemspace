from rdkit import Chem
import pandas as pd
import os

DATAPATH = "../data"


def add_relative_freq(df1, df2):
    total = len(df1)
    df2["rel_freq"] = [int(x) / total for x in df2["abs_freq"].tolist()]
    return df2


def sdf_to_smiles_table(sdf_path, df1):
    smiles_list = []
    names = []
    freqs = []
    supplier = Chem.SDMolSupplier(sdf_path)
    for mol in supplier:
        if mol is None:
            continue  # Skip any invalid molecules
        smiles = Chem.MolToSmiles(mol)
        smiles_list.append(smiles)
        name = mol.GetProp("name") if mol.HasProp("name") else None
        names.append(name)
        freq = mol.GetProp("freq") if mol.HasProp("freq") else None
        freqs.append(freq)
    df = pd.DataFrame({"smiles": smiles_list, "name": names, "abs_freq": freqs})
    df = add_relative_freq(df1, df)
    return df


npsdf = os.path.join(DATAPATH,"original", "scaffolds", "scaffold_np_fin3_recap.sdf")
df = pd.read_csv(os.path.join(DATAPATH, "all_molecules.csv"))
npdf = df[df["category"] == "natural"]
sdsdf = os.path.join(DATAPATH,"original", "scaffolds", "scaffold_sd_fin3_recap.sdf")
sddf = df[df["category"] == "synthetic"]
dbsdf = os.path.join(DATAPATH, "drugbank_scaffolds.sdf")
dbdf = pd.read_csv(os.path.join(DATAPATH, "drugbank_smiles.csv"))
np = sdf_to_smiles_table(npsdf, npdf)
sd = sdf_to_smiles_table(sdsdf, sddf)
db = sdf_to_smiles_table(dbsdf, dbdf)

print(np.shape, sd.shape, db.shape)
np.to_csv(os.path.join(DATAPATH, "np_scaffolds.csv"), index=False)
sd.to_csv(os.path.join(DATAPATH, "sd_scaffolds.csv"), index=False)
db.to_csv(os.path.join(DATAPATH, "db_scaffolds.csv"), index=False)

np = np.rename(columns={"abs_freq": "np_abs_freq", "rel_freq": "np_rel_freq"})
sd = sd.rename(columns={"abs_freq": "sd_abs_freq", "rel_freq": "sd_rel_freq"})
db = db.rename(columns={"abs_freq": "db_abs_freq", "rel_freq": "db_rel_freq"})
print(len(list(set(np["name"]) - set(sd["name"]))))
print(len(list(set(np["name"]) - set(db["name"]))))
print(len(list(set(db["name"]) - set(sd["name"]))))

merged_df = pd.merge(np, sd, on=["smiles", "name"], how="outer")
merged_df = pd.merge(merged_df, db, on=["smiles", "name"], how="outer")
merged_df = merged_df.fillna(0)
print(merged_df.shape, len(set(merged_df["name"])))
print(merged_df.head())

merged_df.to_csv(os.path.join(DATAPATH, "all_scaffolds.csv"), index=False)
