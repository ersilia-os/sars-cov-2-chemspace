# REDIAL
ersilia serve eos8fth
ersilia -v run -i ../data/drugbank_smiles.csv -o ../data/drugbank_eos8fth.csv
ersilia close

#ImageMol
ersilia serve eos4cxk
ersilia -v run -i ../data/drugbank_smiles.csv -o ../data/drugbank_eos4cxk.csv
ersilia close

#SarsCov1 Chemprop
ersilia serve eos9f6t
ersilia -v run -i ../data/drugbank_smiles.csv -o ../data/drugbank_eos9f6t.csv
ersilia close