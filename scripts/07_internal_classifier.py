import os
import pandas as pd
import random
import collections
import lazyqsar as lq
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from tqdm import tqdm
import json
import numpy as np
from rdkit import Chem
import random

TIME_BUDGET_SEC = 60
ESTIMATOR_LIST = ["rf", "lgbm"]
REDUCED = False
N_FOLDS = 5

root = os.path.dirname(os.path.abspath(__file__))

if not os.path.exists(os.path.join(root, "..", "results", "classifiers")):
    os.makedirs(os.path.join(root, "..", "results", "classifiers"))


def smiles_to_inchikey(smiles_list, forbidden_inchikeys=[]):
    smiles_list_ = []
    inchikeys = []
    for smi in tqdm(smiles_list):
        ik = Chem.MolToInchiKey(Chem.MolFromSmiles(smi))
        if ik is None:
            continue
        if ik in forbidden_inchikeys:
            continue
        inchikeys += [ik]
        smiles_list_ += [smi]
    return smiles_list_, inchikeys


class ClassifierPipeline(object):
    def __init__(self, name, dp, du, max_neg_to_pos=10):
        self.name = name
        smiles_p = dp["smiles"].tolist()
        if "inchikey" not in list(dp.columns):
            smiles_p, iks_p = smiles_to_inchikey(smiles_p)
        else:
            iks_p = dp["inchikey"].tolist()
        smiles_u = du["smiles"].tolist()
        if "inchikey" not in list(du.columns):
            smiles_u, iks_u = smiles_to_inchikey(smiles_u)
        else:
            iks_u = du["inchikey"].tolist()
        if len(smiles_u) > max_neg_to_pos * len(smiles_p):
            idxs = list(range(len(smiles_u)))
            idxs = random.sample(idxs, max_neg_to_pos * len(smiles_p))
            smiles_u = [smiles_u[i] for i in idxs]
            iks_u = [iks_u[i] for i in idxs]
        self.smiles_1 = smiles_p
        self.iks_1 = iks_p
        self.smiles_0 = smiles_u
        self.iks_0 = iks_u
        smiles_list = self.smiles_0 + self.smiles_1
        y = [0] * len(self.smiles_0) + [1] * len(self.smiles_1)
        idxs = list(range(len(smiles_list)))
        random.shuffle(idxs)
        self.smiles_list = [smiles_list[i] for i in idxs]
        self.y = [y[i] for i in idxs]
        self.model = None

    def cross_validate(self):

        cross_validation_data = collections.defaultdict(list)

        for _ in tqdm(range(N_FOLDS)):
            smiles_train, smiles_test, y_train, y_test = train_test_split(
                self.smiles_list, self.y, test_size=0.2, stratify=self.y
            )

            model = lq.MorganBinaryClassifier(
                reduced=REDUCED, time_budget_sec=TIME_BUDGET_SEC, estimator_list=ESTIMATOR_LIST
            )
            model.fit(smiles_train, y_train)

            y_pred = model.predict_proba(smiles_test)[:, 1]

            fpr, tpr, thr = roc_curve(y_test, y_pred)
            roc_auc = auc(fpr, tpr)
            J = tpr - fpr

            best_thr_index = np.argmax(J)
            best_thr = thr[best_thr_index]

            cross_validation_data["roc_auc"] += [roc_auc]
            cross_validation_data["thr"] += [best_thr]
            cross_validation_data["y_hat"] += [list(y_pred)]
            cross_validation_data["y"] += [list(y_test)]

        with open(
            os.path.join(root, "..", "results", "classifiers", f"{self.name}.json"),
            "w",
        ) as f:
            json.dump(cross_validation_data, f, indent=4)

    def train(self):
        self.model = lq.MorganBinaryClassifier(
            reduced=REDUCED, time_budget_sec=TIME_BUDGET_SEC, estimator_list=ESTIMATOR_LIST
        )
        self.model.fit(self.smiles_list, self.y)

    def predict(self, name, df):
        print(f"Predicting {name}", df.shape)
        smiles_list = df["smiles"].tolist()
        inchikeys = df["inchikey"].tolist()
        y_hat = self.model.predict_proba(smiles_list)[:, 1]
        file_name = os.path.join(root, "..", "results", "classifiers", f"classifier_output_{name}.csv")
        if os.path.exists(file_name):
            df = pd.read_csv(file_name)
        else:
            df = pd.DataFrame({"inchikey": inchikeys, "smiles": smiles_list})
        df[self.name] = y_hat
        df.to_csv(file_name, index=False)

print("Gathering prediction datasets")
df_manually_annotated = pd.read_csv(os.path.join(root, "..", "data", "all_molecules_with_cheese_results.csv"))[['smiles', 'inchikey']]
print("... manually annotated", df_manually_annotated.shape)
df_cheese_search = pd.read_csv(os.path.join(root, "..", "results", "cheese_search_smiles_and_inchikey.csv"))[['smiles', 'inchikey']]
print("... CHEESE search SMILES and InChIKeys", df_cheese_search.shape)
df_drugbank = pd.read_csv(os.path.join(root, "..", "data", "drugbank_inchikeys_with_cheese_results.csv"))[['smiles', 'inchikey']]
print("... Drugbank", df_drugbank.shape)
df = pd.read_csv(os.path.join(root, "..", "results", "cheese_drugbank_search.csv"))
smiles = df["smiles"].tolist()
inchikeys = df["inchikey"].tolist()
pairs = [(smi, ik) for smi, ik in zip(smiles, inchikeys)]
pairs = list(set(pairs))
smiles = [p[0] for p in pairs]
inchikeys = [p[1] for p in pairs]
df_drugbank_cheese_search = pd.DataFrame({'smiles': smiles, 'inchikey': inchikeys}).sample(n=df_cheese_search.shape[0])
print("... Drugbank CHEESE search", len(inchikeys), "subsampled", df_drugbank_cheese_search.shape)

print("Manually annotated vs Drugbank")
df_p = pd.read_csv(os.path.join(root, "..", "data", "all_molecules_with_cheese_results.csv"))
df_u = pd.read_csv(os.path.join(root, "..", "data", "drugbank_inchikeys_with_cheese_results.csv"))
classifier = ClassifierPipeline("ma_vs_db", df_p, df_u)
classifier.cross_validate()
classifier.train()
classifier.predict("ma", df_manually_annotated)
classifier.predict("macs", df_cheese_search)
classifier.predict("db", df_drugbank)
classifier.predict("dbcs", df_drugbank_cheese_search)

print("Manually annotated vs ChEMBL (0)")
df_p = pd.read_csv(os.path.join(root, "..", "data", "all_molecules_with_cheese_results.csv"))
df_u = pd.read_csv(os.path.join(root, "..", "data", "reference_library_inchikeys.csv"))
classifier = ClassifierPipeline("ma_vs_ch0", df_p, df_u)
classifier.cross_validate()
classifier.train()
classifier.predict("ma", df_manually_annotated)
classifier.predict("macs", df_cheese_search)
classifier.predict("db", df_drugbank)
classifier.predict("dbcs", df_drugbank_cheese_search)

print("Manually annotated vs ChEMBL (1)")
df_p = pd.read_csv(os.path.join(root, "..", "data", "all_molecules_with_cheese_results.csv"))
df_u = pd.read_csv(os.path.join(root, "..", "data", "reference_library_inchikeys.csv"))
classifier = ClassifierPipeline("ma_vs_ch1", df_p, df_u)
classifier.cross_validate()
classifier.train()
classifier.predict("ma", df_manually_annotated)
classifier.predict("macs", df_cheese_search)
classifier.predict("db", df_drugbank)
classifier.predict("dbcs", df_drugbank_cheese_search)

print("Manually annotated vs ChEMBL (2)")
df_p = pd.read_csv(os.path.join(root, "..", "data", "all_molecules_with_cheese_results.csv"))
df_u = pd.read_csv(os.path.join(root, "..", "data", "reference_library_inchikeys.csv"))
classifier = ClassifierPipeline("ma_vs_ch2", df_p, df_u)
classifier.cross_validate()
classifier.train()
classifier.predict("ma", df_manually_annotated)
classifier.predict("macs", df_cheese_search)
classifier.predict("db", df_drugbank)
classifier.predict("dbcs", df_drugbank_cheese_search)

print("Chemdiv Covid vs Chemdiv")
df_p = pd.read_csv(os.path.join(root, "..", "data", "chemdiv", "chemdiv_molecules.csv"))
df_u = pd.read_csv(os.path.join(root, "..", "data", "chemdiv", "chemdiv_100k_generalistic.csv"))
classifier = ClassifierPipeline("cdc_vs_cdg", df_p, df_u)
classifier.cross_validate()
classifier.train()
classifier.predict("ma", df_manually_annotated)
classifier.predict("macs", df_cheese_search)
classifier.predict("db", df_drugbank)
classifier.predict("dbcs", df_drugbank_cheese_search)
