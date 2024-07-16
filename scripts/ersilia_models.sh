## assuming all models are fetched

#ersilia compound embeddings
ersilia serve eos2gw4 
ersilia -v run -i ../data/all_molecules.csv -o ../data/all_molecules_eos2gw4.csv
ersilia close

#RA Score
ersilia serve eos2r5a
ersilia -v run -i ../data/all_molecules.csv -o ../data/all_molecules_eos2r5a.csv
ersilia close

#Morgan FPS
ersilia serve eos4wt0 
ersilia -v run -i ../data/all_molecules.csv -o ../data/all_molecules_eos4wt0.csv
ersilia close

#AdmetAI
ersilia serve eos7d58
ersilia -v run -i ../data/all_molecules.csv -o ../data/all_molecules_eos7d58.csv
ersilia close

#SYBA
ersilia serve eos7pw8
ersilia -v run -i ../data/all_molecules.csv -o ../data/all_molecules_eos7pw8.csv
ersilia close

#NP Score
ersilia serve eos8ioa
ersilia -v run -i ../data/all_molecules.csv -o ../data/all_molecules_eos8ioa.csv
ersilia close

#SA Score
ersilia serve eos9ei3
ersilia -v run -i ../data/all_molecules.csv -o ../data/all_molecules_eos9ei3.csv
ersilia close

#NP Likeness
ersilia serve eos9yui
ersilia -v run -i ../data/all_molecules.csv -o ../data/all_molecules_eos9yui.csv
ersilia close

#Whales
ersilia serve eos24ur
ersilia -v run -i ../data/all_molecules.csv -o ../data/all_molecules_eos24ur.csv
ersilia close