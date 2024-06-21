import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
from _constants import *
from _functions import *


plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["font.size"] = 12

PATH_RESULTS = "../results/milestones_updown"


def plot_multiple_graphs(scaling):
    """
    Plot the similarity mappings for multiple strains.

    Parameters
    ----------
    scaling : str
        The scaling method to use.
    """

    split_strains = [list(strains[:4]), list(strains[4:])]
    split_strains[1].insert(0, "AX4")
    text_buffer = 0.8

    for i, milestone in enumerate(milestone_dict.keys()):
        fig, axs = plt.subplots(1, 2, figsize=(9, 4))
        fig.suptitle(f"{milestone} ({len(milestone_dict[milestone])})", fontsize=12)

        if len(milestone_dict[milestone]) > 0:

            for k, strains_sub in enumerate(split_strains):
                ax = axs[k]
                y_offset = (len(strains_sub) * 1.5) / 2
                ax.set_xticks(range(0, 21, 4))
                ax.set_yticks(range(0, 21, 4))

                plot_bounds = ax.get_xlim() + ax.get_ylim()
                times_fine = np.linspace(0, 1, SPACE)

                all_interpolated_curves = [
                    interpolate_curve(
                        (
                            list(mappings[scaling][milestone][strain].values()),
                            list(mappings[scaling][milestone][strain].keys()),
                        ),
                        times_fine,
                    )
                    for strain in strains_sub
                ]
                existing_labels = []

                for j, strain in enumerate(strains_sub):
                    sorted_X, sorted_Y = list(
                        mappings[scaling][milestone][strain].values()
                    ), list(mappings[scaling][milestone][strain].keys())

                    # original_curve = interpolate_curve((sorted_X, sorted_Y), times_fine)
                    # above_curve = interpolate_curve(
                    #     adjust_curve_points((sorted_X, sorted_Y), 1, -1), times_fine
                    # )
                    # below_curve = interpolate_curve(
                    #     adjust_curve_points((sorted_X, sorted_Y), -1, 1), times_fine
                    # )
                    # x_interpolated, y_interpolated = zip(*original_curve)

                    # name_length = LENGTH_DICT[strain]
                    # combined_label_positions = above_curve + below_curve
                    # losses = [
                    #     loss_function(
                    #         index,
                    #         existing_labels,
                    #         combined_label_positions,
                    #         plot_bounds,
                    #         all_interpolated_curves,
                    #         name_length,
                    #     )[0]
                    #     for index, _ in enumerate(combined_label_positions)
                    # ]
                    # min_loss_index = losses.index(min(losses))

                    # optimal_label_pos = combined_label_positions[min_loss_index]
                    # above = min_loss_index <= SPACE
                    # existing_labels.append(
                    #     (optimal_label_pos, name_length, above),
                    # )

                    # ha, va = (
                    #     ("right", "bottom")
                    #     if min_loss_index < SPACE
                    #     else ("left", "top")
                    # )
                    # min_loss_index_curv = (
                    #     min_loss_index
                    #     if min_loss_index < SPACE
                    #     else min_loss_index - SPACE
                    # )

                    # # Drawing and labeling
                    # if above:
                    #     ax.plot(
                    #         optimal_label_pos[1] - name_length,
                    #         optimal_label_pos[0] + SPACE_VAR * 2,
                    #         "o",
                    #         color="white",
                    #         transform=mtrans.offset_copy(
                    #             ax.transData, fig=fig, x=0.0, y=y_offset, units="points"
                    #         ),
                    #     )
                    # else:
                    #     ax.plot(
                    #         optimal_label_pos[1] + name_length,
                    #         optimal_label_pos[0] - SPACE_VAR * 2,
                    #         "o",
                    #         color="white",
                    #         transform=mtrans.offset_copy(
                    #             ax.transData, fig=fig, x=0.0, y=y_offset, units="points"
                    #         ),
                    #     )
                    # ax.plot(
                    #     y_interpolated,
                    #     x_interpolated,
                    #     color=LINE_COLORS[strain],
                    #     linestyle="-",
                    #     linewidth=1.5,
                    #     label=NAME_DICT[strain],
                    #     transform=mtrans.offset_copy(
                    #         ax.transData, fig=fig, x=0.0, y=y_offset, units="points"
                    #     ),
                    # )
                    # ax.plot(
                    #     [y_interpolated[min_loss_index_curv], optimal_label_pos[1]],
                    #     [x_interpolated[min_loss_index_curv], optimal_label_pos[0]],
                    #     color=LINE_COLORS[strain],
                    #     linestyle="-",
                    #     linewidth=1,
                    #     transform=mtrans.offset_copy(
                    #         ax.transData, fig=fig, x=0.0, y=y_offset, units="points"
                    #     ),
                    # )
                    # ax.text(
                    #     optimal_label_pos[1],
                    #     optimal_label_pos[0],
                    #     NAME_DICT[strain],
                    #     ha=ha,
                    #     va=va,
                    #     fontsize=6,
                    #     transform=mtrans.offset_copy(
                    #         ax.transData, fig=fig, x=0.0, y=y_offset, units="points"
                    #     ),
                    #     color=LINE_COLORS[strain],
                    # )
                    # ax.autoscale(enable=False)
                    label_y_pos = sorted_X[-1]
                    tr = mtrans.offset_copy(
                        ax.transData, fig=fig, x=0.0, y=y_offset, units="points"
                    )
                    ax.plot(
                        sorted_Y,
                        sorted_X,
                        label=f"{NAME_DICT[strain]}",
                        color=LINE_COLORS[strain],
                        transform=tr,
                    )

                    y_offset -= 1.5

                ax.legend(
                    title="strain",
                    fontsize=6,
                    title_fontsize=7,
                    labelspacing=0.2,
                    loc="best",
                )

        # Set a single x-label and y-label for the entire figure
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
            f"{PATH_RESULTS}/{scaling}/{scaling}_milestones_updown_{milestone}.pdf",
            dpi=300,
        )
        plt.close()


# Load data
data = pd.read_csv(DATA_PATH)

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
strain_data_dict = {}
for scaling in SCALING:
    strain_data_dict[scaling] = {}
    for milestone in milestone_dict.keys():
        strain_data_dict[scaling][milestone] = {}

for milestone in milestone_dict.keys():
    for strain in strains:
        strain_df_ax4 = data[data["Strain"] == "AX4"]
        strain_df = data[data["Strain"] == strain]
        strain_df = strain_df.drop(["Strain", "Unnamed: 0"], axis=1)
        strain_df = strain_df.astype(float)
        strain_df = strain_df.filter(items=["Time"] + milestone_dict[milestone])
        strain_df.set_index("Time", inplace=True)
        for scaling in SCALING:
            strain_data_dict[scaling][milestone][strain] = scale_data(
                strain_df, strain_df_ax4, scaling
            )


t_values = np.array(strain_data_dict["None"][MILESTONES[0] + "_up"]["AX4"].index)

mappings = {}
for scaling in SCALING:
    mappings[scaling] = {}
    for milestone in milestone_dict.keys():
        ax4_t_dict = {}
        mappings[scaling][milestone] = {}
        if len(milestone_dict[milestone]) > 0:
            mappings[scaling][milestone] = calculate_MDS_mapping(
                strain_data_dict[scaling][milestone],
                scaling,
                t_values,
                strains,
                PATH_RESULTS
                + f"/individual/{scaling}/{scaling}_{milestone}_",
            )
        else:
            for strain in strains:
                mappings[scaling][milestone][strain] = {}


plot_multiple_graphs("None")
plot_multiple_graphs("m0s1")

print("Your plots are saved in mds_like_analysis/results/milestones_updown/")
