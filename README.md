# SARS-CoV-2 Chemical Space
This repository explores the chemical space associated with SARS-CoV-2 according to manually curated data at [UB-CeDD](https://www.ub-cedd.org/) (Buea, Cameroon).
The project is developed in the context of a research visit of [Prof. Fidele Ntie-Kang](https://scholar.google.de/citations?user=XvORr_kAAAAJ&hl=en).

## Overview
An important component of UB-CeDD's work is the manual curation of synthetic and natural product compounds targetting a specific pathogen. UB-CeDD has assembled a team of curators who have catalogued compounds of potential SARS-CoV-2 activity. In this project, we set up a chemoinformatics pipeline to expand the chemical space around these manually-curated compounds using tools that can work in low-resource settings. Importantly, while SARS-CoV-2 may not be a top priority globally, the pipeline is conceived such that in can be applied to other pathogens and disease areas. Therefore, this repository should be viewed as a proof-of-concept. Central to the pipeline is the [Ersilia Model Hub](https://ersilia.io), our resource of open source AI/ML models for drug discovery.

The pipeline consists of the following steps:
1. Manual curation of SARS-CoV-2 related compounds at UB-CeDD. Both synthetic and natural products are annotated.
1. Automated characterisation of the manually-curated compounds using the Ersilia Model Hub. This includes calculation of interesting features such as synthetic-accessibility and natural-product-likeness, as well as ADME properties.
1. Ultra-large scale similarity search against [Zinc](https://zinc15.docking.org/) and [Enamine REAL](https://enamine.net/compound-collections/real-compounds/real-database) databases. For this, [CHEESE](https://cheese-docs.deepmedchem.com/) is used to query using 2D and 3D similarities.
1. Post-processing and aggregation of similarity search results.
1. Characterisation of the resulting chemical space based on SARS-CoV-2 predictors as provided by the [REDIAL-2020](https://github.com/sirimullalab/redial-2020) suite of models, as available from the Ersilia Model Hub.

## Data
All data used in the project is publicly available. The manually curated molecules can be found under `data/original` and a compilation of both Natural Products and Synthetic Derivatives is available in `all_molecules.csv`. Additional datasets downloaded from their public sources ([DrugBank](https://go.drugbank.com/), [ChemDiv Coronavirus](https://www.chemdiv.com/)) are also in the `/data` folder. In addition to manually curated and publicly available data, the folder contains the calculations of several descriptors and molecular properties for each of the datasets, obtained via the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia). Each dataset is referenced with an `eos` identifier. 

## Scripts
The `/scripts` folder contains the data processing pipelines and the necessary scripts to precalculate all descriptors and molecular properties using Ersilia's models. The Python scripts (`.py`) can be executed one-by-one following their numerical order. The `ersilia_models_*` scripts can be executed assuming all the necessary models have been fetched using the `ersilia` [command-line interface](https://github.com/ersilia-os/ersilia). The `/results` folder contains the molecules obtained during the ultra-large scale library screening using [CHEESE](https://cheese-docs.deepmedchem.com/).

You will need a `.env` file with the `CHEESE_API_KEY`. You need to place this file in the `sars-cov-2-chemspace` directory. The file should have the following information:
```text
CHEESE_API_KEY="your_api_key"
```
To obtain a this key, please follow instructions from the [CHEESE documentation](https://cheese-docs.deepmedchem.com/getting_started-api/). CHEESE is a remarkable resource developed by a small team. Please acknowledge this resource if you find it useful and be mindful of resource consumption when querying their tool.

## Notebooks
The `notebooks` folder contains a relatively unstructure set of Jupyter notebooks showcasing data analysis recipes and publication-ready figure generation. Figures are stored in the `figures` folder.

## Citation
We are currently working on a draft related to this project. Feel free to read this draft in this [online document](https://docs.google.com/document/d/1YgLPCFoM3Zh1ZCTlY_I0X4r0wkm8OyQ1/edit?usp=sharing&ouid=114775674178390159004&rtpof=true&sd=true).

## License
The code in this repository is available under a GPLv3 license and the data and figures under a CC-BY-4 License.

## About us
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit research organisation aimed at developing data science capacity in the global south.
