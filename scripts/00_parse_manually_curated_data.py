import os
import pandas as pd
import rdkit
from standardiser import standardise
from rdkit import RDLogger
from rdkit import Chem
import csv
import shutil
from tqdm import tqdm

root = os.path.dirname(os.path.abspath(__file__))

RDLogger.DisableLog("rdApp.*")

data_dir = os.path.join(root, "..", "data")
np_dir = os.path.join(data_dir, "original", "NP_curated")
sd_dir = os.path.join(data_dir, "original", "SD_curated")


def molecule_loader(subfolder):
    sdf_paths = []
    mol_paths = []
    mol2_paths = []
    non_opened = []

    for fn in os.listdir(subfolder):
        if fn.endswith(".sdf"):
            sdf_paths.append(os.path.join(subfolder, fn))
        if fn.endswith(".mol"):
            mol_paths.append(os.path.join(subfolder, fn))
        if fn.endswith(".mol2"):
            mol2_paths.append(os.path.join(subfolder, fn))

    mols = []
    paths = []
    names = []
    for sdf_path in sdf_paths:
        name = sdf_path.split("/")[-1][:-4]
        suppl = rdkit.Chem.SDMolSupplier(sdf_path)
        mols_ = [mol for mol in suppl if mol is not None]
        if len(mols_) == 0:
            non_opened += [name]
            continue
        if len(mols_) > 1:
            mols_ = [mols_[0]]
        mols += mols_
        paths += [sdf_path]
        names += [name]

    for mol_path in mol_paths:
        name = mol_path.split("/")[-1][:-4]
        mol = rdkit.Chem.MolFromMolFile(mol_path)
        if mol is None:
            non_opened += [name]
            continue
        mols += [mol]
        paths += [mol_path]
        names += [name]

    for mol2_path in mol2_paths:
        name = mol2_path.split("/")[-1][:-5]
        mol = rdkit.Chem.MolFromMol2File(mol2_path)
        if mol is None:
            non_opened += [name]
            continue
        mols += [mol]
        paths += [mol2_path]
        names += [name]
    print(len(sdf_paths), len(mol_paths), len(mol2_paths))
    print(non_opened)
    assert len(mols) == len(names)
    print("TOTAL MOLS", len(mols))
    mols_ = []
    non_parsed_mols = []
    c = 0
    for i, mol in enumerate(mols):
        try:
            mol = standardise.run(mol)
            if mol is not None:
                mols_ += [(names[i], mol)]
        except:
            c += 1
            non_parsed_mols += [names[i]]
            continue
    print(
        "Number of non-standardized molecules (skipped) {0}. File: {1}".format(
            c, subfolder
        )
    )
    return mols_, non_parsed_mols

np_mols, np_non_parsed_mols = molecule_loader(np_dir)
sd_mols, sd_non_parsed_mols = molecule_loader(sd_dir)


np_non_parsed = pd.DataFrame({"file_name": np_non_parsed_mols, "category": "natural"})
sd_non_parsed = pd.DataFrame({"file_name": sd_non_parsed_mols, "category": "synthetic"})
non_parsed = pd.concat([np_non_parsed, sd_non_parsed])
non_parsed.to_csv(os.path.join(data_dir, "original", "non_parsed_mols.csv"), index=False)

def mols_to_table(mols, category):
    mols_ = []
    for name, mol in mols:
        smiles = rdkit.Chem.MolToSmiles(mol)
        mols_ += [
            (
                name,
                rdkit.Chem.MolToInchiKey(Chem.MolFromSmiles(smiles)),
                smiles,
                category,
            )
        ]
    df = pd.DataFrame(mols_, columns=["file_name", "inchikey", "smiles", "category"])
    return df


np_df = mols_to_table(np_mols, "natural")
sd_df = mols_to_table(sd_mols, "synthetic")
np_df_dup = np_df.drop_duplicates(subset=["inchikey"])
print(len(np_df), len(np_df_dup))
sd_df_dup = sd_df.drop_duplicates(subset=["inchikey"])
print(len(sd_df), len(sd_df_dup))

np_duplicates = np_df[np_df.duplicated(subset=["inchikey"], keep=False)]
sd_duplicates = sd_df[sd_df.duplicated(subset=["inchikey"], keep=False)]
duplicates = pd.concat([np_duplicates, sd_duplicates])
duplicates.to_csv(os.path.join(data_dir, "original", "duplicated_mols.csv"), index=False)

df = pd.concat([np_df_dup, sd_df_dup]).reset_index(drop=True)
df_ = df.drop_duplicates(subset=["inchikey"])
assert len(df) == len(df_)
print(len(np_df_dup), len(sd_df_dup), len(df))

df.to_csv(os.path.join(data_dir, "all_molecules.csv"), index=False)


"""
# parse chemdiv

chemdiv_sdf = os.path.join(
    data_dir, "chemdiv", "ChemDiv_SDF_CORONAVIRUS_Library_20750.sdf"
)

suppl = Chem.SDMolSupplier(chemdiv_sdf)
mols = [mol for mol in suppl if mol is not None]

with open(os.path.join(data_dir, "chemdiv", "chemdiv_molecules.csv"), "w") as f:
    writer = csv.writer(f)
    writer.writerow(["name", "inchikey", "smiles"])
    for mol in tqdm(mols):
        name = mol.GetProp("IDNUMBER")
        inchikey = rdkit.Chem.MolToInchiKey(mol)
        smiles = rdkit.Chem.MolToSmiles(mol)
        writer.writerow([name, inchikey, smiles])

# parse chemdiv generalistic molecules

chemdiv_sdf = os.path.join(data_dir, "Preplated-set-100K.sdf")

suppl = Chem.SDMolSupplier(chemdiv_sdf)
mols = [mol for mol in suppl if mol is not None]

with open(os.path.join(data_dir, "chemdiv_100k_generalistic.csv"), "w") as f:
    writer = csv.writer(f)
    writer.writerow(["name", "inchikey", "smiles"])
    for mol in tqdm(mols):
        name = mol.GetProp("IDNUMBER")
        inchikey = rdkit.Chem.MolToInchiKey(mol)
        smiles = rdkit.Chem.MolToSmiles(mol)
        writer.writerow([name, inchikey, smiles])


# parse reference library molecules from chembl

reference_library = os.path.join(
    data_dir, "reference_library.txt"
)  # this file was downloaded as specified in the ersilia-os/compound-embedding repository.

if not os.path.exists(os.path.join(data_dir, "reference_library_inchikeys.csv")):
    with open(reference_library, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        smiles_list = [row[0] for row in reader]

    with open(os.path.join(data_dir, "reference_library_inchikeys.csv"), "w") as f:
        writer = csv.writer(f)
        writer.writerow(["inchikey", "smiles"])
        for smiles in tqdm(smiles_list):
            mol = rdkit.Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            inchikey = rdkit.Chem.MolToInchiKey(mol)
            writer.writerow([inchikey, smiles])

"""