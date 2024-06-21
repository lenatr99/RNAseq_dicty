# Analysis of *Dictyostelium discoideum* RNA-seq data

This repository contains the code and data for the analysis of RNA-seq data from *Dictyostelium discoideum*.

## Differential expression analysis
To run the code, first navigate to the directory and make a conda environment with the required packages:

```bash
cd differential_expression
conda create --name differential_expression python=3.11
conda activate differential_expression
pip install -r requirements.txt
```

The website with the acquired results is available [here](https://lenatr99.pythonanywhere.com/). But if you want to run the code locally, you can do it by running the following command:

```bash
python main.py
```
The website will be available at http://127.0.0.1:5000.


-------------------

To reproduce the results, you can run the scripts in directory `differential_expression/code`. The analysis is divided into pairwise and time series comparisons. The pairwise comparisons are done by comparing the expression of each gene between wild type and each mutant for each time point. The time series comparisons are done by comparing the time series expression of each gene between wild type and each mutant using the spline regression model.

### Pairwise analysis
The pairwise comparisons are done by comparing the expression of each gene between wild type and each mutant for each time point. The results are stored in the `differential_expression/data/pairwise` directory. The code for the pairwise comparisons is in the `code` directory and have `pw` prefix.

Run the following command to reproduce the results:

```bash
cd code
python pw_01_DE_genes.py # to get the differentially expressed genes
python pw_02_DE_genesets.py # to get the differentially expressed gene sets
python pw_03_DE_network.py # to get the infomap/louvain communities
python pw_04_DE_network_analysis.py # to get the network visualization
```

Before running the last script, make sure to use ChatGPT or another tool to generate the common function of the gene sets in the communities. The template for that is saved in the `communities_desc_template.json` file. Save the output of the ChatGPT model in the `communities_desc.json` file.

### Time series analysis
The time series comparisons are done by comparing the time series expression of each gene between wild type and each mutant using the spline regression model. The results are stored in the `differential_expression/data/time_series` directory. The code for the time series comparisons is in the `code` directory and have `ts` prefix.

To reproduce the results, tou must first run the R script `ts_01_DE_spline_r_script.R` to get the differentially expressed genes. Then run the following Python script to get the differentially expressed gene sets:

```bash
cd code
python ts_02_DE_genes.py # to get the differentially expressed gene sets
python ts_03_DE_network.py # to get the infomap/louvain communities
python ts_04_DE_network_analysis.py # to get the network visualization
```
Like in the pairwise analysis, make sure to use ChatGPT or another tool to generate the common function of the gene sets in the communities. The template for that is saved in the `communities_desc_template.json` file. Save the output of the ChatGPT model in the `communities_desc.json` file.

