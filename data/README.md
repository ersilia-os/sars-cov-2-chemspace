# Raw data

The `all` folder was downloaded from Ersilia's Google Drive as provided by Prof. Ntie-Kang on the 7th of July, 2024.
This list corresponds to manually-curated SARS-Cov-2 compounds. Manual curation was done by Prof. Ntie-Kang's team at University of Buea.

# Molecules

The `all_molecules.csv` file contains all molecules as processed by `00_ParseManuallyCuratedData.ipynb`.

# Ersilia models

We have then generated files for the following Ersilia models:
- [eos7d58](https://github.com/ersilia-os/eos7d58): `admet-ai-prediction`
- [eos24ur](https://github.com/ersilia-os/eos24ur): `whales-scaled`
- [eos4wt0](https://github.com/ersilia-os/eos4wt0): `morgan-fps`
- [eos2gw4](https://github.com/ersilia-os/eos2gw4): `eosce`
- [eos9yui](https://github.com/ersilia-os/eos9yui): `natural-product-likeness`
- [eos8ioa](https://github.com/ersilia-os/eos8ioa): `natural-product-score`
- [eos2r5a](https://github.com/ersilia-os/eos2r5a): `retrosynthetic-accessibility`
- [eos7pw8](https://github.com/ersilia-os/eos7pw8): `syba-synthetic-accessibility`
- [eos9ei3](https://github.com/ersilia-os/eos9ei3): `sa-score`

Results as saved as `all_molecules_$MODEL_ID.csv`.

# MOE data

The files `np2_moe.csv` and `sd2_moe.csv` have been provided by Prof. Ntie-Kang on the 9th of July 2024.

