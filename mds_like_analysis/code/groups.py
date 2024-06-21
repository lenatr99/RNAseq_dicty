import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
from mds_like_analysis.code._constants import *
from mds_like_analysis.code._functions import *
from PIL import Image

plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["font.size"] = 12

PATH_DATA = "Data/"
data_file = "B1C1rgB_AX4L846F_B1L846F_ave.csv"

PATH_RESULTS = "all_results/Results_mds/"


def plot_multiple_graphs(scaling, group):
    plots = []
    for i, subgroup in enumerate(gene_annotations[group].keys()):
        fig, axs = plt.subplots(1, 2, figsize=(9, 4))
        split_strains = [list(strains[:4]), list(strains[4:])]
        split_strains[1].insert(0, "AX4")
        text_buffer = 0.8
        fig.suptitle(f"{subgroup} ({len(gene_annotations[group][subgroup])})")
        for k, strains_sub in enumerate(split_strains):
            ax = axs[k]
            y = (len(strains_sub) * 1.5) / 2
            text_positions = []
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

            ax.legend(title="strain", fontsize=6, title_fontsize=7, labelspacing=0.2)

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
        plt.savefig(
            f"{PATH_RESULTS}/groups/{scaling}/subplots/{subgroup}_{scaling}.png",
            dpi=300,
        )
        plots.append(
            f"{PATH_RESULTS}/groups/{scaling}/subplots/{subgroup}_{scaling}.png"
        )
        plt.close()

    images = [Image.open(image_path) for image_path in plots]
    total_height = sum(image.size[1] for image in images)
    max_width = max(image.size[0] for image in images)
    concatenated_image = Image.new("RGB", (max_width, total_height))
    current_height = 0
    for image in images:
        concatenated_image.paste(image, (0, current_height))
        current_height += image.size[1]
    concatenated_image.save(
        f"{PATH_RESULTS}/groups/{scaling}/MDS_{group}_{scaling}.pdf",
        "PDF",
        resolution=300.0,
    )


# Load data
data = pd.read_csv(PATH_DATA + data_file)

# Load gene annotation data
geneannot_data = pd.read_csv(PATH_DATA + "DictyGeneAnnotations_3504.csv", index_col=0)
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
                + f"groups/individual/{scaling}/{scaling}_{group}_{subgroup}_",
                modified=False,
                plot=True,
            )


for group in geneannot_dict.keys():
    plot_multiple_graphs("m0s1", group)
    plot_multiple_graphs("None", group)
