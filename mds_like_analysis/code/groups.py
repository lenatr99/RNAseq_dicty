"""
This Python script analyzes and visualizes gene expression data of different groups
by loading the data, extracting unique strains, and scaling each strain's data relative
to a reference strain ("AX4") using Z-score normalization. It calculates custom mappings
for each strain and scaling method based on a cost function to measure similarity over
time. The script then generates and saves plots to visualize these patterns, facilitating
comparative analysis of gene expression across different groups.
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from _constants import SCALING, DATA_PATH_ALL, DATA_PATH, ANNOT_PATH
from _functions import (
    scale_data,
    calculate_mds_like_mapping,
    calculate_mds_mapping,
    calculate_pca_mapping,
    plot_multiple_graphs_groups,
)

parser = argparse.ArgumentParser(description="Get mappings for all genes")
parser.add_argument(
    "--mapping_mode",
    type=str,
    default="mds-like",
    choices=["mds-like", "mds", "pca"],
    help="Choose mapping mode: mds-like, mds, or pca. Default: mds-like",
)
args = parser.parse_args()

mapping_mode = args.mapping_mode
if mapping_mode != "mds-like":
    print(f"Using {mapping_mode} for mapping")

plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["font.size"] = 12

PATH_RESULTS = "../results/groups"
if not os.path.exists(PATH_RESULTS):
    os.makedirs(PATH_RESULTS)

# Load data
data = pd.read_csv(DATA_PATH_ALL)
data_avg = pd.read_csv(DATA_PATH)

# Load gene annotation data
geneannot_data = pd.read_csv(ANNOT_PATH, index_col=0)
geneannot_data.drop(index=["string", "meta"], inplace=True)
geneannot_dict = {}
gene_annotations = {}
# get all unique values for each column
COLUMNS = [
    "tgr_disagg_gene",
    "dediff",
    "Parikh_cell-type",
    "Maeda_cell-type",
    "Maeda_detailed cell-type",
    "cAMP-pulse induced",
    "chemotaxis",
    "sig and  sigN genes",
    "gtaG dependent short protein",
]

for column in COLUMNS:
    geneannot_dict[column] = geneannot_data[column].unique()
    # get rid of nan
    geneannot_dict[column] = geneannot_dict[column][~pd.isnull(geneannot_dict[column])]
    gene_annotations[column] = {}
    if len(geneannot_dict[column]) > 1:
        gene_annotations[column]["entire_group"] = set()
    for subgroup in geneannot_dict[column]:
        gene_annotations[column][subgroup] = geneannot_data[
            geneannot_data[column] == subgroup
        ].index.tolist()
    if len(geneannot_dict[column]) > 1:
        for value in gene_annotations[column].values():
            gene_annotations[column]["entire_group"].update(value)
        gene_annotations[column]["entire_group"] = list(
            gene_annotations[column]["entire_group"]
        )


# Separate data for different strains
strains = data["Strain"].unique()
replicates = [
    "avg",
] + list(range(3))
strain_data_dict = {}
for replicate in replicates:
    strain_data_dict[replicate] = {}
    for scaling in SCALING:
        strain_data_dict[replicate][scaling] = {}
        for group in geneannot_dict:
            strain_data_dict[replicate][scaling][group] = {}
            for subgroup in gene_annotations[group].keys():
                strain_data_dict[replicate][scaling][group][subgroup] = {}

for replicate in replicates:
    for group in geneannot_dict:
        for subgroup in gene_annotations[group].keys():
            strain_df_AX4 = data_avg[data_avg["Strain"] == "AX4"]
            strain_df_AX4 = strain_df_AX4.drop(["Strain", "Unnamed: 0"], axis=1)
            for strain in strains:
                if replicate != "avg":
                    strain_df = data[data["Strain"] == strain]
                    strain_df = strain_df[strain_df["Replicate"] == replicate]
                    strain_df = strain_df.drop(
                        ["Replicate", "Unnamed: 0", "Strain"], axis=1
                    )
                else:
                    strain_df = data_avg[data_avg["Strain"] == strain]
                    strain_df = strain_df.drop(["Strain", "Unnamed: 0"], axis=1)
                strain_df = strain_df.astype(float)
                strain_df = strain_df.filter(
                    items=["Time"] + gene_annotations[group][subgroup]
                )
                strain_df.set_index("Time", inplace=True)
                for scaling in SCALING:
                    strain_data_dict[replicate][scaling][group][subgroup][strain] = (
                        scale_data(strain_df, strain_df_AX4, scaling)
                    )

t_values = np.array(
    strain_data_dict["avg"]["None"][COLUMNS[0]][geneannot_dict[COLUMNS[0]][0]][
        "AX4"
    ].index
)

mappings = {}
for replicate in replicates:
    mappings[replicate] = {}
    for scaling in SCALING:
        mappings[replicate][scaling] = {}
        for group in geneannot_dict:
            mappings[replicate][scaling][group] = {}
            for subgroup in gene_annotations[group].keys():
                ax4_t_dict = {
                    t: np.array(
                        strain_data_dict["avg"][scaling][group][subgroup]["AX4"].iloc[i]
                    )
                    for i, t in enumerate(t_values)
                }
                if mapping_mode == "mds-like":
                    mappings[replicate][scaling][group][subgroup] = (
                        calculate_mds_like_mapping(
                            strain_data_dict[replicate][scaling][group][subgroup],
                            t_values,
                            strains,
                            PATH_RESULTS
                            + f"/individual/{scaling}/{scaling}_{group}_{subgroup}_",
                            ax4_t_dict=ax4_t_dict,
                        )
                    )
                elif mapping_mode == "mds":
                    mappings[replicate][scaling][group][subgroup] = (
                        calculate_mds_mapping(
                            strain_data_dict[replicate][scaling][group][subgroup],
                            t_values,
                            strains,
                            ax4_t_dict=ax4_t_dict,
                        )
                    )
                elif mapping_mode == "pca":
                    mappings[replicate][scaling][group][subgroup] = (
                        calculate_pca_mapping(
                            strain_data_dict[replicate][scaling][group][subgroup],
                            t_values,
                            strains,
                            ax4_t_dict=ax4_t_dict,
                        )
                    )


for group in geneannot_dict:
    plot_multiple_graphs_groups(
        "m0s1",
        group,
        gene_annotations,
        mappings,
        strains,
        PATH_RESULTS,
        t_values,
        mapping_mode,
    )
    plot_multiple_graphs_groups(
        "None",
        group,
        gene_annotations,
        mappings,
        strains,
        PATH_RESULTS,
        t_values,
        mapping_mode,
    )


print("Your plots are saved in mds_like_analysis/results/groups/")
