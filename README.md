# Analysis of *Dictyostelium discoideum* RNA-seq Data

This repository contains the code and data for the analysis of RNA-seq data from *Dictyostelium discoideum* mutants.

First, create a conda environment with the required packages:

```bash
conda create --name dicty_analysis python=3.11
conda activate dicty_analysis
pip install -r requirements.txt
```

## MDS-like Analysis

To run the MDS-like analysis code and generate the plots, run the following commands:

```bash
cd mds_like_analysis/code
python all_genes.py # MDS-like analysis for all genes
python milestones.py # MDS-like analysis across different milestone gene sets
python milestones_updown.py # MDS-like analysis across different upregulated and downregulated milestone gene sets
python groups.py # MDS-like analysis across different special groups of genes
```

The results will be saved in the `mds_like_analysis/results` directory.

## Differential Expression Analysis

Reproduce the results by running the scripts in the `differential_expression/code` directory. The analysis is divided into pairwise and time series comparisons.

### Pairwise Analysis

Pairwise comparisons involve comparing the expression of each gene between wild type and each mutant at each time point. The results are stored in the `differential_expression/results/pairwise` directory. The code for the pairwise comparisons is in the `differential_expression/code` directory and has the `pw` prefix.

Run the following commands to reproduce the results:

```bash
cd differential_expression/code
python pw_01_DE_genes.py # Identify differentially expressed genes
python pw_02_DE_genesets.py # Identify differentially expressed gene sets
python pw_03_DE_network.py # Generate infomap/louvain communities
python pw_04_DE_network_analysis.py # Visualize the network
```
Before running the last script, use ChatGPT or another tool to generate the common functions of the gene sets in the communities. The template for this is saved in the `communities_desc_template.json` file. Save the output in the `communities_desc.json` file. The tool used in our case was [ChatGPT, version 4](https://openai.com/index/gpt-4/). 

### Time Series Analysis

Time series comparisons involve comparing the time series expression of each gene between wild type and each mutant using a spline regression model. The results are stored in the `differential_expression/results/time_series` directory. The code for the time series comparisons is in the `differential_expression/code` directory and has the `ts` prefix.

To reproduce the results, first run the R script `ts_01_DE_spline_r_script.R` in R Studio (or another R environment) to identify the differentially expressed genes. Then run the following Python scripts to get the differentially expressed gene sets and visualize the network:

```bash
cd differential_expression/code
python ts_02_DE_genes.py # Identify differentially expressed gene sets
python ts_03_DE_network.py # Generate infomap/louvain communities
python ts_04_DE_network_analysis.py # Visualize the network
```

Like in the pairwise analysis, use ChatGPT or another tool to generate the common functions of the gene sets in the communities. The template for this is saved in the `communities_desc_template.json` file. Save the output in the `communities_desc.json` file. The tool used in our case was [ChatGPT, version 4](https://openai.com/index/gpt-4/). 

----------------

The website with the acquired results is available [here](https://lenatr99.pythonanywhere.com/). If you want to run the code locally, you can do so by running the following command:

```bash
python main.py
```

The website will be accessible at [http://127.0.0.1:5000](http://127.0.0.1:5000).