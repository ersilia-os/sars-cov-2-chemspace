#REDIAL
ersilia serve eos8fth
ersilia -v run -i ../results/cheese_search_flat_ranked.csv -o ../results/cheese_eos8fth.csv
ersilia close

#SYBA
ersilia serve eos7pw8
ersilia -v run -i ../results/cheese_search_flat_ranked.csv -o ../results/cheese_eos7pw8.csv
ersilia close

#NP Score
ersilia serve eos8ioa
ersilia -v run -i ../results/cheese_search_flat_ranked.csv -o ../results/cheese_eos8ioa.csv
ersilia close

#ersilia compound embeddings
ersilia serve eos2gw4 
ersilia -v run -i ../results/cheese_search_flat_ranked.csv -o ../results/cheese_eos2gw4.csv
ersilia close

#Whales
ersilia serve eos24ur
ersilia -v run -i ../results/cheese_search_flat_ranked.csv -o ../results/cheese_eos24ur.csv
ersilia close

#UniMol
ersilia serve eos39co
ersilia -v run -i ../results/cheese_search_flat_ranked.csv -o ../results/cheese_eos39co.csv
ersilia close

#Morgan FPS
ersilia serve eos4wt0 
ersilia -v run -i ../results/cheese_search_flat_ranked.csv -o ../results/cheese_eos4wt0.csv
ersilia close



