import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.spatial import distance
from _constants import *


def scale_data(data, wildtype_data, scaling):
    """
    Scale the data of variant strains against the wild type reference using Z-score normalization.

    Parameters
    ----------
    data : pandas.DataFrame
        Data of the variant strain to be scaled.
    wildtype_data : pandas.DataFrame
        Reference wild type data used for calculating the mean and standard deviation.
    scaling : str
        Indicates the type of scaling. For this function, it should always be 'm0s1'.

    Returns
    -------
    scaled_data : pandas.DataFrame
        Scaled data.
    """
    scaled_data = data.copy()
    if scaling == "m0s1":
        for column in data.columns:
            if wildtype_data[column].std() != 0:
                scaled_data[column] = (
                    data[column] - wildtype_data[column].mean()
                ) / wildtype_data[column].std()
            else:
                scaled_data[column] = 0
    return scaled_data


def plot_graph(
    x_values, costs, strain, time, i, euclidean_distances, scaling, PATH_RESULTS
):
    """
    Plot the cost function for a given strain and time point to visualize if there's a minimum.

    Parameters
    ----------
    x_values : numpy.ndarray
        The x values to plot.
    costs : numpy.ndarray
        The cost values to plot.
    strain : str
        The strain to plot.
    time : float
        The time point to plot.
    i : int
        The index of the time point.
    """

    plt.figure(figsize=(6, 4))
    plt.plot(x_values, costs, label="Cost Function")
    plt.xticks(range(0, 21, 4))
    plt.xlabel("$p_x$")
    plt.ylabel("$\sum_t (||a_t - x|| - |t - p_x|)^2$")
    list_str = ", ".join(map(str, np.round(euclidean_distances[strain][i].tolist(), 2)))
    plt.text(
        1,
        costs[0],
        f"minmax D: {list_str}",
        horizontalalignment="left",
        verticalalignment="bottom",
        fontsize=8,
    )
    plt.savefig(f"{PATH_RESULTS}{strain}_{time}.pdf")
    plt.close()


def line_plot(strains, dict):
    """
    Plot the line plots for multiple strains.

    Parameters
    ----------
    strains : list
        List of strains to compare and plot.
    dict : dict
        Dictionary containing the gene expressions.
    method : str
        Method to use for calculating similarity score. Options are 'knn' (k-nearest neighbors), 'pearson' (Pearson correlation coefficient) and
        'mds' (multidimensional scaling).

    """
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 14

    MILESTONE_TO_CHECK = ["tag_tip_up", "tag_tip_down"]
    for scaling in SCALING:
        for strain in strains:
            fig, ax = plt.subplots(figsize=(9, 4))
            ax.set_xticks(range(0, 20 + 1, 4))
            ax.set_xticklabels(range(0, 20 + 1, 4))

            added_labels = set()
            custom_handles = []
            for j, milestone in enumerate(MILESTONE_TO_CHECK):
                for k in range(len(dict[scaling][milestone][strain].values[0])):
                    X_times = dict[scaling][milestone][strain].index
                    Y_times = dict[scaling][milestone][strain].values[:, k]
                    if milestone not in added_labels:
                        ax.plot(
                            X_times,
                            Y_times,
                            label=f"{milestone}",
                            color=MILESTONE_COLORS[milestone],
                            alpha=0.2,
                            linewidth=0.5,
                        )
                        custom_handles.append(
                            Line2D(
                                [0],
                                [0],
                                color=MILESTONE_COLORS[milestone],
                                lw=2,
                                label=milestone,
                            )
                        )
                        added_labels.add(milestone)
                    else:
                        ax.plot(
                            X_times,
                            Y_times,
                            color=MILESTONE_COLORS[milestone],
                            alpha=0.2,
                            linewidth=0.5,
                        )

            ax.set_title(f"{scaling}", fontsize=10)
            ax.legend(
                handles=custom_handles,
                title="milestone",
                fontsize=9,
                title_fontsize=10,
                loc="upper right",
            )

            # Set a single x-label and y-label for the entire figure
            plt.xlabel("hours of mutant development", fontsize=12)
            plt.tight_layout()

            plt.savefig(f"{PATH_RESULTS}/line_plots/line_plot_{scaling}_{strain}.pdf")
            plt.close()


def calculate_MDS_mapping(
    strain_data_dict,
    scaling,
    time_values,
    strains,
    results_path,
    plot=False,
    reference_strain="AX4",
    modified=False,
    radius=2,
):
    mappings = {}
    euclidean_distances = {}
    ax4_t_dict = {
        t: np.array(strain_data_dict[reference_strain].iloc[i])
        for i, t in enumerate(time_values)
    }

    for strain in strains:
        euclidean_distances[strain] = {}
        mappings[strain] = {}

        for i, time in enumerate(time_values):
            x = np.array(strain_data_dict[strain].iloc[i])

            # Calculate euclidean distances between
            distances = [np.linalg.norm(ax4_t_dict[t] - x) for t in time_values]
            distances = np.array(distances)
            min_value, max_value = np.min(distances), np.max(distances)
            normalized_distances = (
                (distances - min_value) / (max_value - min_value)
                if max_value != min_value
                else np.zeros(len(distances))
            )
            euclidean_distances[strain][i] = normalized_distances

            # Calculate px values and corresponding distance sums
            if modified:
                if i > 0:
                    px_values = np.linspace(
                        max(0, mappings[strain][(i - 1) * 4] / 4 - radius),
                        min(
                            len(time_values) - 1,
                            mappings[strain][(i - 1) * 4] / 4 + radius,
                        ),
                        (len(time_values) - 1) * 100 + 1,
                    )
                else:
                    px_values = np.linspace(
                        max(0, i - radius),
                        min(len(time_values) - 1, i + radius),
                        (len(time_values) - 1) * 100 + 1,
                    )
            else:
                px_values = np.linspace(
                    0, len(time_values) - 1, (len(time_values) - 1) * 100 + 1
                )

            dists = {}
            for px in px_values:
                dists[px] = [abs(j - px) for j in range(len(time_values))]
                min_dist, max_dist = np.min(dists[px]), np.max(dists[px])
                dists[px] = (
                    (dists[px] - min_dist) / (max_dist - min_dist)
                    if max_dist != min_dist
                    else np.zeros(len(dists[px]))
                )

            sums = [
                sum(
                    [
                        (euclidean_distance - dists[px][t]) ** 2
                        for t, euclidean_distance in enumerate(normalized_distances)
                    ]
                )
                for px in px_values
            ]

            # Find the minimum sum and update the mapping
            mappings[strain][i * 4] = px_values[np.argmin(sums)] * 4
            if plot:
                plot_graph(
                    px_values * 4,
                    sums,
                    strain,
                    time,
                    i,
                    euclidean_distances,
                    scaling,
                    results_path,
                )

    return mappings


def interpolate_curve(points, times_fine):
    """Interpolate x and y points."""
    t = np.linspace(0, 1, len(points[0]))  # Assume points as (x_points, y_points)
    interpolated_x = interp1d(t, points[0])(times_fine)
    interpolated_y = interp1d(t, points[1])(times_fine)
    return list(zip(interpolated_x, interpolated_y))


def adjust_curve_points(points, adjust_x=0, adjust_y=0):
    """Adjust curve points."""
    return ([x + adjust_x for x in points[0]], [y + adjust_y for y in points[1]])


def uniqueness_penalty(x, y, all_interpolated_curves, above=True):
    """
    Calculate the uniqueness of the curve at the given x by considering the difference
    in y from other curves at the same x. The more unique (i.e., different in y),
    the lower the penalty should be.
    """
    if above:
        x = x - 1
        y = y + 1
    else:
        x = x + 1
        y = y - 1
    penalties = []
    for curve in all_interpolated_curves:

        # Find the y value on the current curve at the same x position
        y_values_at_same_x = [point[0] for point in curve if point[1] == y]
        # print(curve)
        # If there is no exact match in x, we can interpolate or skip
        if not y_values_at_same_x:
            interpolated_y = np.interp(y, [p[1] for p in curve], [p[0] for p in curve])
            y_values_at_same_x = [interpolated_y]

        # Calculate penalties based on the y difference
        for other_y in y_values_at_same_x:
            # The larger the difference, the smaller the penalty
            penalties.append(abs(x - other_y))

    penalties.remove(min(penalties))
    # calculate standard deviation of penalties
    penalties_std = np.std(penalties)
    penalties_mean = np.mean(penalties)

    if not penalties:
        return 0
    else:
        return penalties_mean  # Smaller differences in y result in larger penalties


def loss_function(
    label_index,
    other_labels,
    combined_label_positions,
    plot_bounds,
    all_interpolated_curves,
    name_length=10,
    space=SPACE,
    label_buffer=3,
    bounds=2,
    curve_buffer=1.41,
):
    int_label_index = int(np.round(label_index))  # Round and convert to int
    above = True

    # Ensure the index is within the bounds of the label positions list
    if int_label_index >= len(combined_label_positions):
        int_label_index = len(combined_label_positions) - 1
    elif int_label_index < 0:
        int_label_index = 0

    if int_label_index >= space:
        above = False

    x, y = combined_label_positions[int_label_index]

    loss = 0

    # Penalty for being outside the plot
    if above:
        if (
            x < plot_bounds[0] + 1.5
            or x > plot_bounds[1] - 1.5
            or y < plot_bounds[2] + name_length
            or y > plot_bounds[3] - 1.5
        ):
            loss += 100000
    else:
        if (
            x < plot_bounds[0] + 1.5
            or x > plot_bounds[1] - 1.5
            or y < plot_bounds[2] + 1.5
            or y > plot_bounds[3] - name_length
        ):
            loss += 100000

    # Penalty for overlapping with other labels
    for other_label in other_labels:
        space_adjustment = 1 if other_label[2] else -1
        for other_letter in np.arange(0, other_label[1], 0.5):
            y_adjust = (
                (other_label[0][1] + other_letter * space_adjustment)
                if other_label[2]
                else (other_label[0][1] - other_letter)
            )
            dist = distance.euclidean((x, y), (other_label[0][0], y_adjust))
            if dist < label_buffer:
                loss += (label_buffer - dist) * 100

        direction = -1 if above else 1
        for letter in np.arange(0, name_length, 0.5):
            y_letter_adjust = y + letter * direction
            dist = distance.euclidean((x, y_letter_adjust), other_label[0])
            if dist < label_buffer:
                loss += (label_buffer - dist) * 100

    # Penalty for overlapping with other curves
    for curve in all_interpolated_curves:
        for point in curve:
            dist = distance.euclidean((x, y), point)
            if dist < curve_buffer:
                loss += curve_buffer - dist
            if above:
                for letter in np.arange(0, name_length, 0.25):
                    dist = distance.euclidean((x, y - letter), point)
                    if dist < curve_buffer:
                        loss += (curve_buffer - dist) * 40
            else:
                for letter in np.arange(0, name_length, 0.5):
                    dist = distance.euclidean((x, y + letter), point)
                    if dist < curve_buffer:
                        loss += (curve_buffer - dist) * 40

    # Incorporate uniqueness based on y-value differences at this x
    uniqueness = uniqueness_penalty(x, y, all_interpolated_curves, above)
    loss -= uniqueness * 50  # Add the uniqueness penalty to the loss

    return loss, uniqueness
