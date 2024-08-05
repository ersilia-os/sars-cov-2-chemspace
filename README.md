# SARS-CoV-2 Chemical Space
This repository explores the chemical space associated with SARS-CoV-2 according to manually curated data at [UB-CeDD](https://www.ub-cedd.org/) (Buea, Cameroon).
The project is developed in the context of a research visit of [Prof. Fidele Ntie-Kang](https://scholar.google.de/citations?user=XvORr_kAAAAJ&hl=en).

## Data
All data used in the project is publicly available. The manually curated molecules can be found under data/original and a compilation of both Natural Products and Synthetic Derivatives is available in all_molecules.csv. Additional datasets downloaded from their public sources (DrugBank, ChemDiv Coronavirus) are also in the /data folder.
In addition to manually curated and internet-retrieved data, the folder contains the calculations of several descriptors and molecular properties for each of the datasets, obtained via the Ersilia Model Hub. Each dataset is referenced with an EOS ID. 

The /results folder contains the molecules obtained during the Ultra Large Library Screening using [Cheese](https://cheese-docs.deepmedchem.com/)

## Scripts
The /scripts folder contains the data curation pipelines and the necessary scripts to precalculate all descriptors and molecular properties using Ersilia.

## Notebooks
The notebooks folder showcases the data analysis and figure generation for the publication.

## Citation
Please cite this work as [TBC]

## Usage
You will need a .env file with the CHEESE_API_KEY **TO DO**

## License
The code in this repository is available under a GPLv3 license and the data and figures under a CC-BY-4 License.

## About us
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit research organisation aimed at developing data science capacity in the global south.
