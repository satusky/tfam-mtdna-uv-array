# TFAM-mtDNA protein binding microarray (PBM) analysis
This repository contains the code used to derive data and figures in [https://doi.org/10.7554/eLife.108862.3]
>King Dillon E, Beard Emily E, Satusky Matthew J, George Alex, Ryde Ian, Johnson Caitlin, Dolan Emma L, Zhang Yuning, Zhu Wei, Wilkins Hunter, Corden Evan, Murphy Susan K, Erie Dorothy, Gordân Raluca, Meyer Joel N (2025) UV irradiation alters TFAM binding specificity and compaction of DNA eLife2026; 14:RP108862

## Requirements
Running notebooks locally requires (jupyter)[https://jupyter.org/]. Original figures were generated using Python 3.11, but the code is compatible with Python >= 3.10, < 3.13.
Using a virtual environment is recommended. To set up the environment, use one of the following options:
venv:
```bash
python -m venv tfam-uv-pbm
source tfam-uv-pbm/bin/activate
pip install -r requirements.txt
```
conda:
```bash
conda env create -f environment.yml
conda activate tfam-uv-pbm
```

## PBM dataframe creation
- file: `tfam_uv_pbm_df_creation.ipynb`
- data: `data/df_creation/`
- figures in paper:
    - Figure 3 figure supplement 2
- description: Notebook used for creating a Pandas DataFrame from the PBM outputs


## PBM analysis
- file: `tfam_uv_pbm_analysis.ipynb`
- data: `data/`
- figures in paper:
    - Figure 3 main
    - Figure 3 figure supplement 3
    - Figure 3 figure supplement 4
    - Figure 3 figure supplement 5
    - Figure 3 figure supplement 8
    - Figure 3 figure supplement 9
- description: Notebook used for main PBM analysis and figure generation

## fiber-seq comparison
- file: `fiber_seq_data_comparison.ipynb`
- data: `data/`, `data/fiber_seq_comparison/`
- figures in paper:
    - Figure 3 main
    - Figure 3 figure supplement 3
    - Figure 3 figure supplement 4
    - Figure 3 figure supplement 5
    - Figure 3 figure supplement 8
    - Figure 3 figure supplement 9
- description: Notebook used for comparing high and low affinity regions identified by [https://doi.org/10.1038/s41594-024-01225-6] 
>Isaac RS, Tullius TW, Hansen KG, Dubocanin D, Couvillion M, Stergachis AB, Churchman LS. 2024. Single-nucleoid architecture reveals heterogeneous packaging of mitochondrial DNA. Nature Structural & Molecular Biology 31:568–577. 10.1038/s41594-024-01225-6, 38347148

## Footprinting simulation
- file: `footprinting_simulation.py`
- data: `data/`, `data/footprinting_comparison`
- figures in paper:
    - Figure 3 figure supplement 6
- description: Script comparing high affinity regions identified by [https://doi.org/10.1101/gr.230409.117] 
>Blumberg A, Danko CG, Kundaje A, Mishmar D. 2018. A common pattern of DNase I footprinting throughout the human mtDNA unveils clues for a chromatin-like organization. Genome Research 28:1158–1168. 10.1101/gr.230409.117, 30002158


## GEO creation
- file: `create_geo_dataset.py`
- data: `data/`
- description: Converts CSV to GEO dataset