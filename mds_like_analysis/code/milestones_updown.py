"""
This Python script analyzes and visualizes milestone data (differentiating between upregulation
and downregulation) by loading the data, extracting unique strains, and scaling each strain's
data relative to a reference strain ("AX4") using Z-score normalization. It calculates custom
mappings for each strain and scaling method based on a cost function to measure similarity over
time. The script then generates and saves plots to visualize these patterns, facilitating
comparative analysis of milestone expression across different strains
"""

import os

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from _constants import DATA_PATH, DATA_PATH_ALL, ANNOT_PATH, MILESTONES, SCALING
from _functions import (
    scale_data,
    calculate_mds_like_mapping,
    plot_multiple_graphs,
    calculate_mds_mapping,
    calculate_pca_mapping,
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

PATH_RESULTS = "../results/milestones_updown"
if not os.path.exists(PATH_RESULTS):
    os.makedirs(PATH_RESULTS)

# Load data
data = pd.read_csv(DATA_PATH_ALL)
data_avg = pd.read_csv(DATA_PATH)

# Load milestone data
milestone_data = pd.read_csv(ANNOT_PATH, index_col=0)
milestone_data.drop(index=["string", "meta"], inplace=True)
milestone_dict = {}
for gene in MILESTONES:
    milestone_dict[gene + "_up"] = milestone_data[
        (milestone_data[gene] == f"{gene}_up")
    ].index.tolist()
    milestone_dict[gene + "_down"] = milestone_data[
        (milestone_data[gene] == f"{gene}_down")
    ].index.tolist()

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
        for milestone in milestone_dict:
            strain_data_dict[replicate][scaling][milestone] = {}

for replicate in replicates:
    for milestone, genes in milestone_dict.items():
        for strain in strains:
            strain_df_AX4 = data_avg[data_avg["Strain"] == "AX4"]
            strain_df_AX4 = strain_df_AX4.drop(["Strain", "Unnamed: 0"], axis=1)

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
            strain_df = strain_df.filter(items=["Time"] + genes)
            strain_df.set_index("Time", inplace=True)
            for scaling in SCALING:
                strain_data_dict[replicate][scaling][milestone][strain] = scale_data(
                    strain_df, strain_df_AX4, scaling
                )

combinations = list(strain_data_dict.keys())[1:]
t_values = np.array(strain_data_dict["avg"]["None"][MILESTONES[0] + "_up"]["AX4"].index)

mappings = {}

for replicate in replicates:
    mappings[replicate] = {}
    for scaling in SCALING:
        mappings[replicate][scaling] = {}
        for milestone, genes in milestone_dict.items():
            ax4_t_dict = {}
            mappings[replicate][scaling][milestone] = {}
            if len(genes) > 0:
                ax4_t_dict = {
                    t: np.array(
                        strain_data_dict["avg"][scaling][milestone]["AX4"].iloc[i]
                    )
                    for i, t in enumerate(t_values)
                }
                if mapping_mode == "mds-like":
                    mappings[replicate][scaling][milestone] = (
                        calculate_mds_like_mapping(
                            strain_data_dict[replicate][scaling][milestone],
                            t_values,
                            strains,
                            PATH_RESULTS
                            + f"/individual/{scaling}/{scaling}_{milestone}_",
                            ax4_t_dict=ax4_t_dict,
                        )
                    )
                elif mapping_mode == "mds":
                    mappings[replicate][scaling][milestone] = calculate_mds_mapping(
                        strain_data_dict[replicate][scaling][milestone],
                        t_values,
                        strains,
                        ax4_t_dict=ax4_t_dict,
                    )
                elif mapping_mode == "pca":
                    mappings[replicate][scaling][milestone] = calculate_pca_mapping(
                        strain_data_dict[replicate][scaling][milestone],
                        t_values,
                        strains,
                        ax4_t_dict=ax4_t_dict,
                    )
            else:
                for strain in strains:
                    mappings[replicate][scaling][milestone][strain] = {}
plot_multiple_graphs(
    "None", strains, mappings, milestone_dict, PATH_RESULTS, t_values, mapping_mode
)
plot_multiple_graphs(
    "m0s1", strains, mappings, milestone_dict, PATH_RESULTS, t_values, mapping_mode
)

print("Your plots are saved in mds_like_analysis/results/milestones_updown/")
