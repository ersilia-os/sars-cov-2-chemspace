import os
import pandas as pd
import rdkit
from standardiser import standardise
from rdkit import RDLogger

root = os.path.dirname(os.path.abspath(__file__))

RDLogger.DisableLog('rdApp.*')                                                                                                                                                           

data_dir = os.path.join(root, '..', 'data')
np_dir = os.path.join(data_dir, 'all', 'NP2')
sd_dir = os.path.join(data_dir, 'all', 'SD2')

def molecule_loader(subfolder):
    sdf_paths = []
    mol_paths = []
    mol2_paths = []
    names = []
    for fn in os.listdir(subfolder):
        if fn.endswith(".sdf"):
            sdf_paths.append(os.path.join(subfolder, fn))
            names += [fn[:-4]]
        if fn.endswith(".mol"):
            mol_paths.append(os.path.join(subfolder, fn))
            names += [fn[:-4]]
        if fn.endswith(".mol2"):
            mol2_paths.append(os.path.join(subfolder, fn))
            names += [fn[:-5]]
    mols = []
    for sdf_path in sdf_paths:
        mols.append(rdkit.Chem.SDMolSupplier(sdf_path))
    for mol_path in mol_paths:
        mols.append(rdkit.Chem.MolFromMolFile(mol_path))
    for mol2_path in mol2_paths:
        mols.append(rdkit.Chem.MolFromMol2File(mol2_path))
    mols_ = []
    c = 0
    for i, mol in enumerate(mols):
        try:
            mol = standardise.run(mol)
            if mol is not None:
                mols_ += [(names[i], mol)]
        except:
            c += 1
            continue
    print("Number of non-standardized molecules (skipped) {0}".format(c))
    return mols_

np_mols = molecule_loader(np_dir)
sd_mols = molecule_loader(sd_dir)

def mols_to_table(mols, category):
    mols_ = []
    for name, mol in mols:
        mols_ += [(name, rdkit.Chem.MolToInchiKey(mol), rdkit.Chem.MolToSmiles(mol), category)]
    df = pd.DataFrame(mols_, columns=['file_name', "inchikey", 'smiles', "category"])
    return df

np_df = mols_to_table(np_mols, "natural")
sd_df = mols_to_table(sd_mols, "synthetic")

df = pd.concat([np_df, sd_df]).drop_duplicates().reset_index(drop=True)

df.to_csv(os.path.join(data_dir, 'all_molecules.csv'), index=False)