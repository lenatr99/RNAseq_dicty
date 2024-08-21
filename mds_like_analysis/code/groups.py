"""
This Python script analyzes and visualizes gene expression data of different groups by loading the data, extracting unique strains, and scaling each strain's data relative to a reference strain ("AX4") using Z-score normalization. It calculates custom mappings for each strain and scaling method based on a cost function to measure similarity over time. The script then generates and saves plots to visualize these patterns, facilitating comparative analysis of gene expression across different groups.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
import os
from _constants import *
from _functions import *
from svgutils.compose import Figure, SVG
import cairosvg

plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["font.size"] = 12

PATH_RESULTS = "../results/groups"
if not os.path.exists(PATH_RESULTS):
    os.makedirs(PATH_RESULTS)


def plot_multiple_graphs(scaling, group):
    plots = []
    for i, subgroup in enumerate(gene_annotations[group].keys()):
        fig, axs = plt.subplots(1, 2, figsize=(9, 4))
        split_strains = [list(strains[:4]), list(strains[4:])]
        split_strains[1].insert(0, "AX4")
        fig.suptitle(f"{subgroup} ({len(gene_annotations[group][subgroup])})")
        for k, strains_sub in enumerate(split_strains):
            ax = axs[k]
            y = (len(strains_sub) * 1.5) / 2
            ax.set_xticks(range(0, 21, 4))
            ax.set_yticks(range(0, 21, 4))

            for strain in strains_sub:
                X_times, Y_times = list(
                    mappings[scaling][group][subgroup][strain].values()
                ), list(mappings[scaling][group][subgroup][strain].keys())
                label_y_pos = X_times[-1]
                tr = mtrans.offset_copy(
                    ax.transData, fig=fig, x=0.0, y=y, units="points"
                )
                ax.plot(
                    Y_times,
                    X_times,
                    label=f"{NAME_DICT[strain]}",
                    color=LINE_COLORS[strain],
                    transform=tr,
                )
                y -= 1.5

            ax.legend(fontsize=6, title_fontsize=7, labelspacing=0.2, frameon=False)

        fig.text(0.5, 0.04, f"hours of mutant development", ha="center", fontsize=12)
        fig.text(
            0.04,
            0.5,
            f"hours of AX4 development",
            va="center",
            rotation="vertical",
            fontsize=12,
        )

        plt.subplots_adjust(left=0.1, bottom=0.15, right=0.88, wspace=0.34)
        if not os.path.exists(f"{PATH_RESULTS}/{scaling}/subplots"):
            os.makedirs(f"{PATH_RESULTS}/{scaling}/subplots")
        plt.savefig(
            f"{PATH_RESULTS}/{scaling}/subplots/{subgroup}_{scaling}.svg", format="svg"
        )
        plots.append(f"{PATH_RESULTS}/{scaling}/subplots/{subgroup}_{scaling}.svg")
        plt.close()

    figures = [SVG(plot) for plot in plots]
    for i, figure in enumerate(figures):
        figure.move(0, i * 288)
    concatenated_figure = Figure(f"648pt", f"{288*len(figures)}pt", *figures)
    concatenated_svg_path = f"{PATH_RESULTS}/{scaling}/{scaling}_groups_{group}.svg"
    concatenated_figure.save(concatenated_svg_path)
    pdf_path = f"{PATH_RESULTS}/{scaling}/{scaling}_groups_{group}.pdf"
    cairosvg.svg2pdf(url=concatenated_svg_path, write_to=pdf_path)


# Load data
data = pd.read_csv(DATA_PATH)

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
strains = data["Strain"].unique()[:-1]
strain_data_dict = {}
for scaling in SCALING:
    strain_data_dict[scaling] = {}
    for group in geneannot_dict.keys():
        strain_data_dict[scaling][group] = {}
        for subgroup in gene_annotations[group].keys():
            strain_data_dict[scaling][group][subgroup] = {}

for group in geneannot_dict.keys():
    for subgroup in gene_annotations[group].keys():
        strain_df_AX4 = data[data["Strain"] == "AX4"]
        for strain in strains:
            strain_df = data[data["Strain"] == strain]
            strain_df = strain_df.drop(["Strain", "Unnamed: 0"], axis=1)
            strain_df = strain_df.astype(float)
            strain_df = strain_df.filter(
                items=["Time"] + gene_annotations[group][subgroup]
            )
            strain_df.set_index("Time", inplace=True)
            for scaling in SCALING:
                strain_data_dict[scaling][group][subgroup][strain] = scale_data(
                    strain_df, strain_df_AX4, scaling
                )

t_values = np.array(
    strain_data_dict["None"][COLUMNS[0]][geneannot_dict[COLUMNS[0]][0]]["AX4"].index
)

mappings = {}
euclidean_distances = {}
for scaling in SCALING:
    euclidean_distances[scaling] = {}
    mappings[scaling] = {}
    ax4_t_dict = {}
    for group in geneannot_dict.keys():
        ax4_t_dict[group] = {}
        euclidean_distances[scaling][group] = {}
        mappings[scaling][group] = {}
        for subgroup in gene_annotations[group].keys():
            mappings[scaling][group][subgroup] = calculate_MDS_mapping(
                strain_data_dict[scaling][group][subgroup],
                scaling,
                t_values,
                strains,
                PATH_RESULTS
                + f"/individual/{scaling}/{scaling}_{group}_{subgroup}_",
            )


for group in geneannot_dict.keys():
    plot_multiple_graphs("m0s1", group)
    plot_multiple_graphs("None", group)


print("Your plots are saved in mds_like_analysis/results/groups/")