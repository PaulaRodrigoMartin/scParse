# PARSE sequencing data analysis --> PARSE pipeline and STARsolo

## Overview
This project involves preprocessing, analyzing and processing single-cell sequencing data. This README provides guidelines on running the necessary pipelines and scripts for the project.

## Prerequisites
Ensure you have the following dependencies installed:
- Python 3.10
- Jupyter Notebook
- Slurm (for cluster computing)
- Required Python libraries (environment.yaml)

## Setup
To run a Jupyter notebook in the cluster, execute the following command in your local terminal:
```sh
slurm-jupyter -e jupyter -A account_of_your_project -u your_username
```
For more details, refer to [Slurm Jupyter documentation](https://slurm-jupyter.readthedocs.io/).

## Workflow

Steps 2,3 and 4 are the steps that are included in workflow.py, so the information given here is merely explanatory. 

### 1. Preprocess FastQ Files
#### Extract and Concatenate Files
1. Untar files from:
   ```sh
   /testis_singlecell/backup/PrimaryData/PARSE
   ```
2. Concatenate lanes per sublibrary:
   ```sh
   cat $(ls | grep "_1.fq.gz" | sort) > PARSE2_UDI_WT_2_R1.fq.gz
   cat $(ls | grep "_2.fq.gz" | sort) > PARSE2_UDI_WT_2_R2.fq.gz
   ```
3. Move files to `primary_data`:
   ```sh
   mv PARSE2_UDI_WT_2_R{1,2}.fq.gz /testis_singlecell/Workspaces/{your_username}/primary_data
   ```

### 2. Split FastQ Files by Sample
Execute the function `split_fq_by_well()` in `workflow.py`.

### 3. Run the PARSE Pipeline
For each sublibrary and sample:
```sh
reference_parse()
spipe_parse()
```
To combine sublibraries:
```sh
comb_parse()
```

### 4. Run the Starsolo Pipeline
For each sublibrary and sample:
```sh
reference_starsolo()
mapping_starsolo()
```
To combine sublibraries:
```sh
# Combine filtered matrices, barcodes, and features
combined_fil == GeneFull/filtered
combined == GeneFull/raw
Velocyto == Velocyto
```

### 5. Process Files for Further Analysis
Run the following scripts, ensuring paths are correctly adjusted:
```sh
python3.10 analysis/combine_files.py
python3.10 analysis/after_mapping_tailored.sh
python3.10 scripts/unite_tables.py
```

### 6. Processing Output Files
- **Collapse random hexamers and polyT:** Run notebooks in `/jupyter/collapse/collapse_primers_{sp}.ipynb`.
- **Quality Control (QC):** Run `/jupyter/QC/preproc_QC_{sp}_starsolo.ipynb`.
- **Clustering:**
  - For all species: `/jupyter/cluster/cluster_intergration.ipynb`
  - For individual species: `/jupyter/cluster/cluster_{sp}.ipynb`

## Contribution
Feel free to submit pull requests or open issues for improvements.

## License
This project is licensed under [Your License].
