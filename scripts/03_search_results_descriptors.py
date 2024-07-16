# from eosce import ErsiliaCompoundEmbeddings
import pandas as pd
import os

root = os.path.dirname(os.path.abspath(__file__))

query_smiles_list = pd.read_csv(os.path.join(root, "../data/all_molecules.csv"))["smiles"].tolist()

emb = ErsiliaCompoundEmbeddings()

smiles_list = pd.read_csv(os.path.join(root, "../results/cheese_search.csv"))["smiles"].tolist()