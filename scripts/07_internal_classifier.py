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

TIME_BUDGET_SEC = 60
ESTIMATOR_LIST = ["rf", "lgbm"]
REDUCED = False
N_FOLDS = 1

root = os.path.dirname(os.path.abspath(__file__))

df = pd.read_csv(os.path.join(root, "..", "data", "chemdiv_molecules.csv"))

smiles_1 = df["smiles"].tolist()
iks_1 = df["inchikey"].tolist()

df = pd.read_csv(os.path.join(root, "..", "data", "all_molecules.csv"))
smiles_1_ = df["smiles"].tolist()
iks_1_ = df["inchikey"].tolist()
for smi, ik in zip(smiles_1_, iks_1_):
    if ik not in iks_1:
        smiles_1 += [smi]
        iks_1 += [ik]

df = pd.read_csv(os.path.join(root, "..", "data", "chemdiv_100k_generalistic.csv"))
df = df[~df["inchikey"].isin(iks_1)]

smiles_0 = df["smiles"].tolist()

smiles_list = smiles_0 + smiles_1
y = [0] * len(smiles_0) + [1] * len(smiles_1)

idxs = list(range(len(smiles_list)))
random.shuffle(idxs)

smiles_list = [smiles_list[i] for i in idxs]
y = [y[i] for i in idxs]

cross_validation_data = collections.defaultdict(list)

for _ in tqdm(range(N_FOLDS)):
    smiles_train, smiles_test, y_train, y_test = train_test_split(
        smiles_list, y, test_size=0.2, random_state=42, stratify=y
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
    os.path.join(root, "..", "results", "internal_classifier_cross_validation.json"),
    "w",
) as f:
    json.dump(cross_validation_data, f, indent=4)

print("Training a full model on the entire dataset")

model = lq.MorganBinaryClassifier(
    reduced=REDUCED, time_budget_sec=TIME_BUDGET_SEC, estimator_list=ESTIMATOR_LIST
)

model.fit(smiles_list, y)

thr = np.mean(cross_validation_data["thr"])


def get_prediction_dataframe(smiles, inchikeys):
    y_pred = model.predict_proba(smiles)[:, 1]
    y_bin = []
    for yp in y_pred:
        if yp > thr:
            y_bin += [1]
        else:
            y_bin += [0]
    df = pd.DataFrame(
        {"inchikey": inchikeys, "smiles": smiles, "y_hat": y_pred, "y_bin": y_bin}
    )
    return df


print("Predicting on manually annotated molecules")

df = pd.read_csv(os.path.join(root, "..", "data", "all_molecules.csv"))
smiles = df["smiles"].tolist()
inchikeys = df["inchikey"].tolist()
df = get_prediction_dataframe(smiles, inchikeys)
df.to_csv(
    os.path.join(root, "..", "results", "internal_classifier_all_molecules.csv"),
    index=False,
)

print("Predicting on CHEESE molecules")

df = pd.read_csv(
    os.path.join(root, "..", "results", "cheese_search_smiles_and_inchikey.csv")
)
smiles = df["smiles"].tolist()
inchikeys = df["inchikey"].tolist()
df = get_prediction_dataframe(smiles, inchikeys)
y_pred = model.predict_proba(smiles)[:, 1]
df.to_csv(
    os.path.join(root, "..", "results", "internal_classifier_cheese_search.csv"),
    index=False,
)

print("Predicting on DrugBank molecules")

df = pd.read_csv(os.path.join(root, "..", "data", "drugbank_inchikeys.csv"))
smiles = df["smiles"].tolist()
inchikeys = df["inchikey"].tolist()
df = get_prediction_dataframe(smiles, inchikeys)
df.to_csv(
    os.path.join(root, "..", "results", "internal_classifier_drugbank.csv"), index=False
)

print("Predicting on Drugbank CHEESE search molecules")

df = pd.read_csv(os.path.join(root, "..", "results", "cheese_drugbank_search.csv"))
smiles = df["smiles"].tolist()
inchikeys = df["inchikey"].tolist()
pairs = [(smi, ik) for smi, ik in zip(smiles, inchikeys)]
pairs = list(set(pairs))
smiles = [p[0] for p in pairs]
inchikeys = [p[1] for p in pairs]
df = get_prediction_dataframe(smiles, inchikeys)
df.to_csv(
    os.path.join(
        root, "..", "results", "internal_classifier_drugbank_cheese_search.csv"
    ),
    index=False,
)
