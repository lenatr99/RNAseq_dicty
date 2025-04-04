"""
This Python script analyzes and visualizes gene expression data by loading the data,
extracting unique strains, and scaling each strain's data relative to a reference strain ("AX4")
using Z-score normalization. It calculates custom mappings for each strain and scaling method
based on a cost function to measure similarity over time. The script then generates and saves plots
to visualize these patterns, facilitating comparative analysis of gene expression across different
strains.
"""

import os
import warnings
from itertools import combinations

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from _constants import DATA_PATH_ALL, DATA_PATH, SCALING, NAME_DICT, LINE_COLORS
from _functions import (
    scale_data,
    calculate_mds_like_mapping,
    get_stars,
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

PATH_RESULTS = "../results/all_genes"
if not os.path.exists(PATH_RESULTS):
    os.makedirs(PATH_RESULTS)

# Load data
data = pd.read_csv(DATA_PATH_ALL)
data_avg = pd.read_csv(DATA_PATH)

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

for replicate in replicates:
    for strain in strains:
        strain_df_AX4 = data_avg[data_avg["Strain"] == "AX4"]
        strain_df_AX4 = strain_df_AX4.drop(["Strain", "Unnamed: 0"], axis=1)
        if replicate != "avg":
            strain_df = data[data["Strain"] == strain]
            strain_df = strain_df[strain_df["Replicate"] == replicate]
            strain_df = strain_df.drop(["Replicate", "Unnamed: 0", "Strain"], axis=1)
        else:
            strain_df = data_avg[data_avg["Strain"] == strain]
            strain_df = strain_df.drop(["Strain", "Unnamed: 0"], axis=1)
        strain_df = strain_df.astype(float)
        strain_df.set_index("Time", inplace=True)
        for scaling in SCALING:
            strain_data_dict[replicate][scaling][strain] = scale_data(
                strain_df, strain_df_AX4, scaling
            )

t_values = np.array(strain_data_dict["avg"]["None"]["AX4"].index)

mappings = {}
for replicate in replicates:
    mappings[replicate] = {}
    for scaling in SCALING:
        ax4_t_dict = {
            t: np.array(strain_data_dict["avg"][scaling]["AX4"].iloc[i])
            for i, t in enumerate(t_values)
        }
        if mapping_mode == "mds-like":
            mappings[replicate][scaling] = calculate_mds_like_mapping(
                strain_data_dict[replicate][scaling],
                t_values,
                strains,
                PATH_RESULTS + f"/individual/{scaling}/{scaling}_",
                ax4_t_dict=ax4_t_dict,
            )
        elif mapping_mode == "mds":
            mappings[replicate][scaling] = calculate_mds_mapping(
                strain_data_dict[replicate][scaling],
                t_values,
                strains,
                ax4_t_dict=ax4_t_dict,
            )
        elif mapping_mode == "pca":
            mappings[replicate][scaling] = calculate_pca_mapping(
                strain_data_dict[replicate][scaling],
                t_values,
                strains,
                ax4_t_dict=ax4_t_dict,
            )

# Calculate ANOVA and Tukey post-hoc tests to get p-values

dfs = {}
for scaling in SCALING:
    rows = []
    for time in t_values:
        for strain in strains:
            for rep in [0, 1, 2]:
                val = mappings[rep][scaling][strain][time]
                rows.append([strain, time, f"{strain}_{rep}", val])

    dfs[scaling] = pd.DataFrame(rows, columns=["Strain", "Time", "Replicate", "Value"])

df_results = {}
for scaling in SCALING:
    df = dfs[scaling]
    df["Strain_Time"] = df["Strain"].astype(str) + "_t" + df["Time"].astype(str)
    grouped = df.groupby("Strain_Time")
    groups = list(grouped.groups.keys())
    pvals = []
    comparisons = []

    for g1, g2 in combinations(groups, 2):
        vals1 = grouped.get_group(g1)["Value"]
        vals2 = grouped.get_group(g2)["Value"]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            if np.var(vals1, ddof=1) < 1e-10 and np.var(vals2, ddof=1) < 1e-10:
                pval = 1.0
            else:
                _, pval = ttest_ind(vals1, vals2, equal_var=False)

        if np.isnan(pval):
            pval = 1.0
        pvals.append(pval)
        comparisons.append((g1, g2))

    # Correct p-values
    _, pvals_corrected, _, _ = multipletests(pvals, method="fdr_bh")

    # Save results as DataFrame
    df_results[scaling] = pd.DataFrame(
        {
            "group1": [c[0] for c in comparisons],
            "group2": [c[1] for c in comparisons],
            "p-raw": pvals,
            "p-adj": pvals_corrected,
        }
    )


p_values_dict = {}
for scaling in SCALING:
    filtered_dataframe = pd.DataFrame()
    p_values_dict[scaling] = {}
    for strain in strains:
        p_values_dict[scaling][strain] = {}
        for time in t_values:
            mask = (
                df_results[scaling]["group1"].str.startswith(f"AX4_t{time}", na=False)
                & df_results[scaling]["group2"].str.startswith(
                    f"{strain}_t{time}", na=False
                )
            ) | (
                df_results[scaling]["group1"].str.startswith(
                    f"{strain}_t{time}", na=False
                )
                & df_results[scaling]["group2"].str.startswith(f"AX4_t{time}", na=False)
            )
            match = df_results[scaling][mask]
            if not match.empty and match.iloc[0]["group1"].startswith(
                f"{strain}_t{time}"
            ):
                match = match.rename(columns={"group1": "group2", "group2": "group1"})
            filtered_dataframe = pd.concat([filtered_dataframe, match])
            if not match.empty:
                p_values_dict[scaling][strain][time] = match.iloc[0]["p-adj"]
            else:
                p_values_dict[scaling][strain][time] = 1.0
    if mapping_mode == "mds-like":
        filtered_dataframe.to_csv(f"{PATH_RESULTS}/{scaling}_p_values.csv", index=False)

for scaling in SCALING:
    fig, axs = plt.subplots(1, 2, figsize=(9, 4))
    split_strains = [list(strains[:4]), list(strains[4:])]
    split_strains[1].insert(0, "AX4")
    for i, strains_sub in enumerate(split_strains):
        ax = axs[i]
        y_offset = (len(strains_sub) * 1.5) / 2
        ax.set_xticks(range(0, 21, 4))
        ax.set_yticks(range(0, 21, 4))
        ax.set_xticklabels(range(0, 21, 4))
        ax.set_yticklabels(range(0, 21, 4))

        plot_bounds = ax.get_xlim() + ax.get_ylim()
        existing_labels = []
        p_values_storage = []
        time_indices = []

        for j, strain in enumerate(strains_sub):
            sorted_x_min = []
            sorted_x_max = []
            sorted_x_avg = []
            replicate_values_at_times = []
            for time in t_values:
                replicate_values = []
                for rep in range(3):
                    replicate_values.append(mappings[rep][scaling][strain][time])

                replicate_values = np.array(replicate_values, dtype=float)
                std_val = np.std(replicate_values, ddof=1) / np.sqrt(
                    len(replicate_values)
                )
                sorted_x_min.append(mappings["avg"][scaling][strain][time] - std_val)
                sorted_x_max.append(mappings["avg"][scaling][strain][time] + std_val)
            sorted_x, sorted_y = list(mappings["avg"][scaling][strain].values()), list(
                mappings["avg"][scaling][strain].keys()
            )

            label_y_pos = sorted_x[-1]
            tr = mtrans.offset_copy(
                ax.transData, fig=fig, x=0.0, y=y_offset, units="points"
            )
            ax.plot(
                sorted_y,
                sorted_x,
                label=f"{NAME_DICT[strain]}",
                color=LINE_COLORS[strain],
                transform=tr,
            )
            min_values = [
                sorted_x[tp] - sorted_x_min[tp] for tp in range(len(sorted_x))
            ]
            max_values = [
                sorted_x_max[tp] - sorted_x[tp] for tp in range(len(sorted_x))
            ]
            for time in t_values:
                p_adj = p_values_dict[scaling][strain][time]
                ax.text(
                    time,
                    y_offset * 0.2 - 2.7,
                    get_stars(p_adj),
                    ha="center",
                    va="bottom",
                    color=LINE_COLORS[strain],
                    fontsize=8,
                    fontweight="bold",
                    transform=tr,
                    clip_on=False,
                )

            ax.errorbar(
                sorted_y,
                sorted_x,
                yerr=[
                    min_values,
                    max_values,
                ],
                fmt="none",
                color=LINE_COLORS[strain],
                transform=tr,
                capsize=4,
                capthick=1,
                alpha=1,
            )

            # # Compute lower and upper bounds for shading
            # lower_bound = [
            #     x - min_val for x, min_val in zip(sorted_x, min_values)
            # ]
            # upper_bound = [x + max_val for x, max_val in zip(sorted_x, max_values)]

            # # Plot shaded area for confidence interval
            # ax.fill_between(
            #     sorted_y,
            #     lower_bound,
            #     upper_bound,
            #     color=LINE_COLORS[strain],
            #     alpha=0.25,
            #     transform=tr,
            #     linewidth=0,
            # )

            y_offset -= 1.5

        ax.legend(
            fontsize=6, title_fontsize=7, labelspacing=0.2, loc="best", frameon=False
        )
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.set_xlim(xlim[0] - 3, xlim[1] + 3)
        ax.set_ylim(ylim[0] - 3, ylim[1] + 3)

    # Set a single x-label and y-label for the entire figure
    fig.text(0.5, 0.04, "hours of mutant development", ha="center", fontsize=12)
    fig.text(
        0.04,
        0.5,
        "hours of AX4 development",
        va="center",
        rotation="vertical",
        fontsize=12,
    )

    plt.subplots_adjust(left=0.1, bottom=0.15, right=0.88, wspace=0.34)
    if not os.path.exists(f"{PATH_RESULTS}/{scaling}"):
        os.makedirs(f"{PATH_RESULTS}/{scaling}")
    if mapping_mode != "mds-like":
        plt.savefig(f"{PATH_RESULTS}/{scaling}/{mapping_mode}_{scaling}_all_genes.pdf")
    else:
        plt.savefig(f"{PATH_RESULTS}/{scaling}/{scaling}_all_genes.pdf")
    plt.close()


print("Your plots are saved in mds_like_analysis/results/all_genes/")
