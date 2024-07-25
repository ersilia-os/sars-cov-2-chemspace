import pandas as pd
import os

#prepare files for Ersilia Models

root = os.path.dirname(os.path.abspath(__file__))


cheese_smiles = pd.read_csv(os.path.join(root, "../results/cheese_search.csv"))["smiles"].tolist()
chemdiv_smiles = pd.read_csv(os.path.join(root, "../data/chemdiv_molecules.csv"))["smiles"].tolist()

print(len(cheese_smiles), len(set(cheese_smiles)))
print(len(chemdiv_smiles), len(set(chemdiv_smiles)))
print(len(set(cheese_smiles) & set(chemdiv_smiles)))

cheese = pd.DataFrame({"smiles": list(set(cheese_smiles))})
print(cheese.shape)
cheese.to_csv(os.path.join(root, "../results/cheese_search_smiles.csv"), index=False)
chemdiv = pd.DataFrame({"smiles": list(set(chemdiv_smiles))})
chemdiv.to_csv(os.path.join(root, "../data/chemdiv_smiles.csv"), index=False)
print(chemdiv.shape)