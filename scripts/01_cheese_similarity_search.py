import os
import pandas as pd
from dotenv import load_dotenv
import time
from rdkit import Chem
import requests
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

root = os.path.dirname(os.path.abspath(__file__))

load_dotenv(os.path.join(root, "..", ".env"))
CHEESE_API_KEY = os.getenv("CHEESE_API_KEY")


def _query_molecule(smiles, search_type, search_quality, n_neighbors):
    inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
    api_key = CHEESE_API_KEY
    data = requests.get(
        "https://api.cheese.themama.ai/molsearch",
        {
            "search_input": smiles,
            "search_type": search_type,
            "n_neighbors": n_neighbors,
            "search_quality": search_quality,
            "descriptors": False,
            "properties": False,
            "filter_molecules": True,
        },
        headers={"Authorization": f"Bearer {api_key}"},
        verify=False,
    ).json()
    R = []
    for r in data["neighbors"]:
        identifier = r["zinc_id"]
        if identifier.startswith("ZINC"):
            database = "zinc15"
        elif identifier.startswith("ENAMINE"):
            database = "enamine-real"
        else:
            database = None
        R += [
            [
                smiles,
                inchikey,
                r["smiles"],
                identifier,
                search_type,
                r["Morgan Tanimoto"],
                database,
            ]
        ]
    df = (
        pd.DataFrame(
            R,
            columns=[
                "query_smiles",
                "query_inchikey",
                "smiles",
                "identifier",
                "search_type",
                "score",
                "database",
            ],
        )
        .sort_values("score", ascending=False)
        .reset_index(drop=True)
    )
    return df


def query_molecule(
    smiles, search_type="consensus", search_quality="very accurate", n_neighbors=100
):
    for _ in range(3):
        try:
            df = _query_molecule(
                smiles,
                search_type=search_type,
                search_quality=search_quality,
                n_neighbors=n_neighbors,
            )
            if df.shape[0] == n_neighbors:
                return df
        except:
            print("Error, retrying in 5 seconds")
            time.sleep(5)
    return None


def query_molecule_all_similarities(
    smiles, search_quality="very accurate", n_neighbors=100
):
    search_types = ["consensus", "morgan", "espsim_electrostatic", "espsim_shape"]
    dfs = []
    for search_type in search_types:
        print(smiles, search_type)
        dfs += [
            query_molecule(
                smiles,
                search_type=search_type,
                search_quality=search_quality,
                n_neighbors=n_neighbors,
            )
        ]
    dfs = [df for df in dfs if df is not None]
    if len(dfs) == 0:
        return None
    df = pd.concat(dfs).reset_index(drop=True)
    return df


def query_and_write(smiles):
    try:
        inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
    except:
        return
    if len(inchikey) != 27:
        print("Invalid inchikey for", smiles)
        return
    file_name = os.path.join("..", "results", "cheese", f"{inchikey}.csv")
    if os.path.exists(file_name):
        print("Already done for", inchikey)
        return
    df = query_molecule_all_similarities(smiles)
    if df is not None:
        df.to_csv(file_name, index=False)


smiles_list = pd.read_csv(os.path.join(root, "../data/all_molecules.csv"))[
    "smiles"
].tolist()

for smiles in smiles_list:
    query_and_write(smiles)

smiles_list_drugbank = pd.read_csv(os.path.join(root, "../data/drugbank_smiles.csv"))[
    "Smiles"
].tolist()

for smiles in smiles_list_drugbank:
    query_and_write(smiles)
