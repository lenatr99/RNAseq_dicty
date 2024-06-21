import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
from _constants import *
from _functions import *


plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["font.size"] = 12

PATH_RESULTS = "all_results/Results_mds/"

# Load data
data = pd.read_csv(DATA_PATH)
strains = data["Strain"].unique()
strain_data_dict = {}
for scaling in SCALING:
    strain_data_dict[scaling] = {}

for strain in strains:
    strain_df_AX4 = data[data["Strain"] == "AX4"]
    strain_df = data[data["Strain"] == strain]
    strain_df = strain_df.drop(["Strain", "Unnamed: 0"], axis=1)
    strain_df = strain_df.astype(float)
    strain_df.set_index("Time", inplace=True)
    for scaling in SCALING:
        strain_data_dict[scaling][strain] = scale_data(
            strain_df, strain_df_AX4, scaling
        )


t_values = np.array(strain_data_dict["None"]["AX4"].index)

mappings = {}
for scaling in SCALING:
    mappings[scaling] = calculate_MDS_mapping(
        strain_data_dict[scaling],
        scaling,
        t_values,
        strains,
        PATH_RESULTS + f"all_strains/individual/{scaling}/{scaling}_",
    )

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
    times_fine = np.linspace(0, 1, SPACE)

    all_interpolated_curves = [
        interpolate_curve(
            (
                list(mappings["None"][strain].values()),
                list(mappings["None"][strain].keys()),
            ),
            times_fine,
        )
        for strain in strains_sub
    ]
    existing_labels = []

    for j, strain in enumerate(strains_sub):
        sorted_X, sorted_Y = list(mappings["None"][strain].values()), list(
            mappings["None"][strain].keys()
        )

        # # Interpolate original, above, and below curves
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
        #     )
        #     for index, _ in enumerate(combined_label_positions)
        # ]
        # min_loss_index = losses.index(min(losses))

        # optimal_label_pos = combined_label_positions[min_loss_index]
        # above = min_loss_index <= SPACE
        # existing_labels.append(
        #     (optimal_label_pos, name_length, above),
        # )

        # ha, va = ("right", "bottom") if min_loss_index < SPACE else ("left", "top")
        # min_loss_index_curv = (
        #     min_loss_index if min_loss_index < SPACE else min_loss_index - SPACE
        # )

        # # Drawing and labeling
        # if above:
        #     ax.plot(
        #         optimal_label_pos[1] - name_length,
        #         optimal_label_pos[0] + SPACE_VAR,
        #         "o",
        #         alpha=0.0,
        #         transform=mtrans.offset_copy(
        #             ax.transData, fig=fig, x=0.0, y=y_offset, units="points"
        #         ),
        #     )
        # else:
        #     ax.plot(
        #         optimal_label_pos[1] + name_length,
        #         optimal_label_pos[0] - SPACE_VAR,
        #         "o",
        #         alpha=0.0,
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
        title="strain", fontsize=6, title_fontsize=7, labelspacing=0.2, loc="best"
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
plt.savefig(f"{PATH_RESULTS}all_strains/MDS_None.pdf")
plt.close()
