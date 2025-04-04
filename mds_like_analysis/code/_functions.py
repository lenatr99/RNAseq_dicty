"""
This Python script contains all the functions used in the main script to analyze and visualize
gene expression data. The functions include scaling data, plotting graphs, calculating bootstrap
confidence intervals, and calculating mappings for each strain based on a cost function. The script
also imports necessary libraries and constants from other files to facilitate the analysis and
visualization of gene expression data.
"""

import os
import warnings
from collections import Counter
from itertools import combinations

import cairosvg
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
import numpy as np
import pandas as pd
from _constants import LINE_COLORS, NAME_DICT, SPACE
from custom_mds import CustomMDS
from scipy.interpolate import interp1d
from scipy.spatial import distance
from scipy.stats import ttest_ind
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
from svgutils.compose import SVG, Figure


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
            if wildtype_data[column].std() != 0:  # Avoid division by zero
                scaled_data[column] = (
                    data[column] - wildtype_data[column].mean()
                ) / wildtype_data[column].std()
            else:
                scaled_data[column] = 0
    return scaled_data


def plot_graph(x_values, costs, strain, time, i, euclidean_distances, path_results):
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
    plt.ylabel(r"$\sum_t (||a_t - x|| - |t - p_x|)^2$")
    list_str = ", ".join(map(str, np.round(euclidean_distances[strain][i].tolist(), 2)))
    plt.text(
        1,
        costs[0],
        f"minmax D: {list_str}",
        horizontalalignment="left",
        verticalalignment="bottom",
        fontsize=8,
    )
    plt.savefig(f"{path_results}{strain}_{time}.pdf")
    plt.close()


def bootstrap_confidence_interval(
    data, num_bootstrap=10000, ci=95, return_samples=False
):
    """
    Calculate bootstrap confidence intervals for the mean of the data.

    Args:
        data: Array-like, the data to bootstrap.
        num_bootstrap: Number of bootstrap samples.
        ci: Confidence interval percentage (default is 95%).

    Returns:
        mean: Mean of the original data.
        lower: Lower bound of the confidence interval.
        upper: Upper bound of the confidence interval.
    """
    bootstrap_means = []
    bootstrap_samples = np.random.choice(data, (num_bootstrap, len(data)), replace=True)

    if return_samples:
        ordered_samples = np.sort(bootstrap_samples, axis=1)
        rows_as_tuples = [tuple(row) for row in ordered_samples]
        counts = Counter(rows_as_tuples)
        return counts
    bootstrap_means = np.mean(bootstrap_samples, axis=1)

    lower = np.percentile(bootstrap_means, (100 - ci) / 2)
    upper = np.percentile(bootstrap_means, 100 - (100 - ci) / 2)
    return [np.mean(data), lower, upper, bootstrap_means]


def calculate_mds_like_mapping(
    strain_data_dict,
    time_values,
    strains,
    results_path,
    plot=False,
    reference_strain="AX4",
    modified=True,
    radius=2,
    ax4_t_dict=None,
):
    """
    Calculate the MDS-like mapping for each strain based on the cost function. This is different
    from the MDS algorithm as distances from specific time points are compared and weighed to find
    best pseudo-time mapping.
    """
    mappings = {}
    euclidean_distances = {}
    if ax4_t_dict is None:
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

            # Calculate px values and corresponding distance sums. if modified,
            # we only consider a small window around the previous px value
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
                    results_path,
                )

    return mappings


def calculate_mds_mapping(
    strain_data_dict,
    time_values,
    strains,
    reference_strain="AX4",
    ax4_t_dict=None,
):
    """
    Calculate the MDS mapping for each strain. This uses Classical MDS to fit on
    reference strain data and then project other strains onto the same space
    using custom function that aproapproximates using distances + Procrustes.

    """
    strictly_ascending = False
    strictly_descending = False
    mappings = {}
    if ax4_t_dict is None:
        ax4_t_dict = {
            t: np.array(strain_data_dict[reference_strain].iloc[i])
            for i, t in enumerate(time_values)
        }

    mds = CustomMDS(n_components=1)
    mds_fit = mds.fit(np.array(list(ax4_t_dict.values())))
    reference_mapping = mds_fit.transform(np.array(list(ax4_t_dict.values())))

    diffs = np.diff(reference_mapping.flatten())
    if np.all(diffs > 0):
        strictly_ascending = True
    elif np.all(diffs < 0):
        strictly_descending = True

    for strain in strains:
        strain_mds = mds_fit.transform(np.array(strain_data_dict[strain]))
        if strain not in mappings:
            mappings[strain] = {}
            for i, time in enumerate(time_values):
                mappings[strain][time] = None
                value = strain_mds[i]
                if strictly_ascending:
                    # Check if within AX4 range
                    if value < reference_mapping[0]:
                        # Extrapolate below first time point
                        x0, x1 = reference_mapping[0], reference_mapping[1]
                        t0, t1 = time_values[0], time_values[1]
                        projected_time = ((value - x0) / (x1 - x0)) * (t1 - t0) + t0
                        mappings[strain][time] = projected_time.item()

                    elif value > reference_mapping[-1]:
                        # Extrapolate above last time point
                        x0, x1 = reference_mapping[-2], reference_mapping[-1]
                        t0, t1 = time_values[-2], time_values[-1]
                        projected_time = ((value - x1) / (x1 - x0)) * (t1 - t0) + t1
                        mappings[strain][time] = projected_time.item()

                    else:
                        # Interpolate within AX4 range
                        for j in range(1, len(reference_mapping)):
                            if (
                                reference_mapping[j - 1]
                                <= value
                                <= reference_mapping[j]
                            ):
                                t0, t1 = time_values[j - 1], time_values[j]
                                x0, x1 = reference_mapping[j - 1], reference_mapping[j]
                                projected_time = ((value - x0) / (x1 - x0)) * (
                                    t1 - t0
                                ) + t0
                                mappings[strain][time] = projected_time.item()
                                break
                elif strictly_descending:
                    # Check if within AX4 range
                    if value > reference_mapping[0]:
                        # Extrapolate above first time point
                        x0, x1 = reference_mapping[0], reference_mapping[1]
                        t0, t1 = time_values[0], time_values[1]
                        projected_time = ((value - x0) / (x1 - x0)) * (t1 - t0) + t0
                        mappings[strain][time] = projected_time.item()

                    elif value < reference_mapping[-1]:
                        # Extrapolate below last time point
                        x0, x1 = reference_mapping[-2], reference_mapping[-1]
                        t0, t1 = time_values[-2], time_values[-1]
                        projected_time = ((value - x1) / (x1 - x0)) * (t1 - t0) + t1
                        mappings[strain][time] = projected_time.item()

                    else:
                        # Interpolate within AX4 range
                        for j in range(1, len(reference_mapping)):
                            if (
                                reference_mapping[j - 1]
                                >= value
                                >= reference_mapping[j]
                            ):
                                t0, t1 = time_values[j - 1], time_values[j]
                                x0, x1 = reference_mapping[j - 1], reference_mapping[j]
                                projected_time = ((value - x0) / (x1 - x0)) * (
                                    t1 - t0
                                ) + t0
                                mappings[strain][time] = projected_time.item()
                                break
                else:
                    # the MDS1 value of the reference strain is not strictly ascending
                    # or descending therefore we cannot interpolate between time points
                    # so we find the closest time point to the value
                    closest = np.argmin(np.abs(reference_mapping - value))
                    mappings[strain][time] = time_values[closest]

    return mappings


def calculate_pca_mapping(
    strain_data_dict,
    time_values,
    strains,
    reference_strain="AX4",
    ax4_t_dict=None,
):
    """
    Calculate the PCA mapping for each strain. This uses PCA to fit on
    reference strain data and then project other strains onto the same space.

    """
    strictly_ascending = False
    strictly_descending = False
    mappings = {}
    if ax4_t_dict is None:
        ax4_t_dict = {
            t: np.array(strain_data_dict[reference_strain].iloc[i])
            for i, t in enumerate(time_values)
        }

    pca = PCA(n_components=1, random_state=42)
    pca_fit = pca.fit(np.array(list(ax4_t_dict.values())))
    reference_mapping = pca_fit.transform(np.array(list(ax4_t_dict.values())))

    diffs = np.diff(reference_mapping.flatten())
    if np.all(diffs > 0):
        strictly_ascending = True
    elif np.all(diffs < 0):
        strictly_descending = True

    for strain in strains:
        strain_mds = pca_fit.transform(np.array(strain_data_dict[strain]))
        if strain not in mappings:
            mappings[strain] = {}
            for i, time in enumerate(time_values):
                mappings[strain][time] = None
                value = strain_mds[i]
                if strictly_ascending:
                    # Check if within AX4 range
                    if value < reference_mapping[0]:
                        # Extrapolate below first time point
                        x0, x1 = reference_mapping[0], reference_mapping[1]
                        t0, t1 = time_values[0], time_values[1]
                        projected_time = ((value - x0) / (x1 - x0)) * (t1 - t0) + t0
                        mappings[strain][time] = projected_time.item()

                    elif value > reference_mapping[-1]:
                        # Extrapolate above last time point
                        x0, x1 = reference_mapping[-2], reference_mapping[-1]
                        t0, t1 = time_values[-2], time_values[-1]
                        projected_time = ((value - x1) / (x1 - x0)) * (t1 - t0) + t1
                        mappings[strain][time] = projected_time.item()

                    else:
                        # Interpolate within AX4 range
                        for j in range(1, len(reference_mapping)):
                            if (
                                reference_mapping[j - 1]
                                <= value
                                <= reference_mapping[j]
                            ):
                                t0, t1 = time_values[j - 1], time_values[j]
                                x0, x1 = reference_mapping[j - 1], reference_mapping[j]
                                projected_time = ((value - x0) / (x1 - x0)) * (
                                    t1 - t0
                                ) + t0
                                mappings[strain][time] = projected_time.item()
                                break
                elif strictly_descending:
                    # Check if within AX4 range
                    if value > reference_mapping[0]:
                        # Extrapolate above first time point
                        x0, x1 = reference_mapping[0], reference_mapping[1]
                        t0, t1 = time_values[0], time_values[1]
                        projected_time = ((value - x0) / (x1 - x0)) * (t1 - t0) + t0
                        mappings[strain][time] = projected_time.item()

                    elif value < reference_mapping[-1]:
                        # Extrapolate below last time point
                        x0, x1 = reference_mapping[-2], reference_mapping[-1]
                        t0, t1 = time_values[-2], time_values[-1]
                        projected_time = ((value - x1) / (x1 - x0)) * (t1 - t0) + t1
                        mappings[strain][time] = projected_time.item()

                    else:
                        # Interpolate within AX4 range
                        for j in range(1, len(reference_mapping)):
                            if (
                                reference_mapping[j - 1]
                                >= value
                                >= reference_mapping[j]
                            ):
                                t0, t1 = time_values[j - 1], time_values[j]
                                x0, x1 = reference_mapping[j - 1], reference_mapping[j]
                                projected_time = ((value - x0) / (x1 - x0)) * (
                                    t1 - t0
                                ) + t0
                                mappings[strain][time] = projected_time.item()
                                break
                else:
                    # the MDS1 value of the reference strain is not strictly ascending
                    # or descending therefore we cannot interpolate between time points
                    # so we find the closest time point to the value
                    closest = np.argmin(np.abs(reference_mapping - value))
                    mappings[strain][time] = time_values[closest]
    return mappings


def interpolate_curve(points, times_fine):
    """Interpolate x and y points."""
    t = np.linspace(0, 1, len(points[0]))  # Assume points as (x_points, y_points)
    interpolated_x = interp1d(t, points[0])(times_fine)
    interpolated_y = interp1d(t, points[1])(times_fine)
    return list(zip(interpolated_x, interpolated_y))


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
    curve_buffer=1.41,
):
    """
    Calculate the loss function for the given label index and other labels. The loss function
    considers the position of the label, the distance from other labels, and the distance from
    other curves.
    """
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


def plot_multiple_graphs(
    scaling,
    strains,
    mappings,
    milestone_dict,
    path_results,
    t_values,
    mapping_mode,
    show_error_bars=True,
    show_significance=True,
):
    """
    Plot the similarity mappings for multiple strains.

    Parameters
    ----------
    scaling : str
        The scaling method to use.
    """

    dfs = {}
    for milestone in milestone_dict.keys():
        if len(milestone_dict[milestone]) > 0:
            rows = []
            for time in t_values:
                for strain in strains:
                    for rep in mappings.keys():
                        if rep != "avg":
                            val = mappings[rep][scaling][milestone][strain][time]
                            rows.append([strain, time, rep, val])

            dfs[milestone] = pd.DataFrame(
                rows, columns=["Strain", "Time", "Replicate", "Value"]
            )

    df_results = {}
    for milestone in milestone_dict.keys():
        if len(milestone_dict[milestone]) > 0:
            df = dfs[milestone]
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
            correction_method = "fdr_bh"  # this is the Benjamini/Hochberg correction
            _, pvals_corrected, _, _ = multipletests(pvals, method=correction_method)

            # Save results as DataFrame
            df_results[milestone] = pd.DataFrame(
                {
                    "group1": [c[0] for c in comparisons],
                    "group2": [c[1] for c in comparisons],
                    "p-raw": pvals,
                    "p-adj": pvals_corrected,
                }
            )

    p_values_dict = {}
    filtered_dataframe = pd.DataFrame()
    for milestone in milestone_dict.keys():
        if len(milestone_dict[milestone]) > 0:
            p_values_dict[milestone] = {}
            for strain in strains:
                p_values_dict[milestone][strain] = {}
                for time in t_values:
                    mask = (
                        df_results[milestone]["group1"].str.startswith(
                            f"AX4_t{time}", na=False
                        )
                        & df_results[milestone]["group2"].str.startswith(
                            f"{strain}_t{time}", na=False
                        )
                    ) | (
                        df_results[milestone]["group1"].str.startswith(
                            f"{strain}_t{time}", na=False
                        )
                        & df_results[milestone]["group2"].str.startswith(
                            f"AX4_t{time}", na=False
                        )
                    )
                    match = df_results[milestone][mask]
                    if not match.empty and match.iloc[0]["group1"].startswith(
                        f"{strain}_t{time}"
                    ):
                        match = match.rename(
                            columns={"group1": "group2", "group2": "group1"}
                        )
                    match.insert(0, "milestone", milestone)
                    filtered_dataframe = pd.concat([filtered_dataframe, match])
                    if not match.empty:
                        p_values_dict[milestone][strain][time] = match.iloc[0]["p-adj"]
                    else:
                        p_values_dict[milestone][strain][time] = 1.0

    if mapping_mode == "mds-like":
        filtered_dataframe.to_csv(f"{path_results}/{scaling}_p_values.csv", index=False)

    split_strains = [list(strains[:4]), list(strains[4:])]
    split_strains[1].insert(0, "AX4")

    for milestone in milestone_dict.keys():
        fig, axs = plt.subplots(1, 2, figsize=(9, 4))
        fig.suptitle(f"{milestone} ({len(milestone_dict[milestone])})", fontsize=12)

        if len(milestone_dict[milestone]) > 0:

            for k, strains_sub in enumerate(split_strains):
                ax = axs[k]
                y_offset = (len(strains_sub) * 1.5) / 2
                ax.set_xticks(range(0, 21, 4))
                ax.set_yticks(range(0, 21, 4))

                for strain in strains_sub:
                    sorted_x_min = []
                    sorted_x_max = []
                    sorted_x_avg = []
                    for time in mappings["avg"][scaling][milestone][strain].keys():
                        replicate_values = [
                            mappings[combination][scaling][milestone][strain][time]
                            for combination in [
                                comb for comb in mappings.keys() if comb != "avg"
                            ]
                        ]
                        sorted_x_avg.append(np.mean(replicate_values))

                        std_val = np.std(replicate_values, ddof=1) / np.sqrt(
                            len(replicate_values)
                        )
                        sorted_x_min.append(
                            mappings["avg"][scaling][milestone][strain][time] - std_val
                        )
                        sorted_x_max.append(
                            mappings["avg"][scaling][milestone][strain][time] + std_val
                        )
                    sorted_x, sorted_y = (
                        list(mappings["avg"][scaling][milestone][strain].values()),
                        list(mappings["avg"][scaling][milestone][strain].keys()),
                    )

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
                    if show_significance:
                        for time in t_values:
                            p_adj = p_values_dict[milestone][strain][time]
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
                    if show_error_bars:
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

                    y_offset -= 1.5

                ax.legend(
                    fontsize=6,
                    title_fontsize=7,
                    labelspacing=0.2,
                    loc="best",
                    frameon=False,
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
        if not os.path.exists(f"{path_results}/{scaling}"):
            os.makedirs(f"{path_results}/{scaling}")
        if mapping_mode != "mds-like":
            plt.savefig(
                f"{path_results}/{scaling}/{mapping_mode}_{scaling}_milestones_updown_{milestone}.pdf",
                dpi=300,
            )
        else:
            plt.savefig(
                f"{path_results}/{scaling}/{scaling}_milestones_updown_{milestone}.pdf",
                dpi=300,
            )
        plt.close()


def plot_multiple_graphs_groups(
    scaling,
    group,
    gene_annotations,
    mappings,
    strains,
    path_results,
    t_values,
    mapping_mode,
    show_error_bars=True,
    show_significance=True,
):
    """
    Plot the similarity mappings for multiple strains in groups.
    """
    dfs = {}
    for subgroup in gene_annotations[group].keys():
        rows = []
        for time in t_values:
            for strain in strains:
                for rep in mappings.keys():
                    if rep != "avg":
                        val = mappings[rep][scaling][group][subgroup][strain][time]
                        rows.append([strain, time, f"{time}_{rep}", val])

        dfs[subgroup] = pd.DataFrame(
            rows, columns=["Strain", "Time", "Replicate", "Value"]
        )

    df_results = {}
    for subgroup in gene_annotations[group].keys():
        df = dfs[subgroup]
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
        correction_method = "fdr_bh"  # this is the Benjamini/Hochberg correction
        _, pvals_corrected, _, _ = multipletests(pvals, method=correction_method)

        # Save results as DataFrame
        df_results[subgroup] = pd.DataFrame(
            {
                "group1": [c[0] for c in comparisons],
                "group2": [c[1] for c in comparisons],
                "p-raw": pvals,
                "p-adj": pvals_corrected,
            }
        )

    p_values_dict = {}
    filtered_dataframe = pd.DataFrame()
    for subgroup in gene_annotations[group].keys():
        p_values_dict[subgroup] = {}
        for strain in strains:
            p_values_dict[subgroup][strain] = {}
            for time in t_values:
                mask = (
                    df_results[subgroup]["group1"].str.startswith(
                        f"AX4_t{time}", na=False
                    )
                    & df_results[subgroup]["group2"].str.startswith(
                        f"{strain}_t{time}", na=False
                    )
                ) | (
                    df_results[subgroup]["group1"].str.startswith(
                        f"{strain}_t{time}", na=False
                    )
                    & df_results[subgroup]["group2"].str.startswith(
                        f"AX4_t{time}", na=False
                    )
                )
                match = df_results[subgroup][mask]
                if not match.empty and match.iloc[0]["group1"].startswith(
                    f"{strain}_t{time}"
                ):
                    match = match.rename(
                        columns={"group1": "group2", "group2": "group1"}
                    )
                match.insert(0, "subgroup", subgroup)
                match.insert(0, "group", group)
                filtered_dataframe = pd.concat([filtered_dataframe, match])
                if not match.empty:
                    p_values_dict[subgroup][strain][time] = match.iloc[0]["p-adj"]
                    pvals.append(match.iloc[0]["p-raw"])
                else:
                    p_values_dict[subgroup][strain][time] = 1.0
        

    if mapping_mode == "mds-like":
        filtered_dataframe.to_csv(
            f"{path_results}/{scaling}_{group}_p_values.csv", index=False
        )

    plots = []
    for i, subgroup in enumerate(gene_annotations[group].keys()):
        fig, axs = plt.subplots(1, 2, figsize=(9, 4))
        split_strains = [list(strains[:4]), list(strains[4:])]
        split_strains[1].insert(0, "AX4")
        fig.suptitle(f"{subgroup} ({len(gene_annotations[group][subgroup])})")
        for k, strains_sub in enumerate(split_strains):
            ax = axs[k]
            y_offset = (len(strains_sub) * 1.5) / 2
            ax.set_xticks(range(0, 21, 4))
            ax.set_yticks(range(0, 21, 4))

            for strain in strains_sub:
                sorted_x_min = []
                sorted_x_max = []
                sorted_x_avg = []
                for time in mappings["avg"][scaling][group][subgroup][strain].keys():
                    replicate_values = [
                        mappings[combination][scaling][group][subgroup][strain][time]
                        for combination in range(3)
                    ]

                    std_val = np.std(replicate_values, ddof=1) / np.sqrt(
                        len(replicate_values)
                    )
                    sorted_x_avg.append(np.mean(replicate_values))
                    sorted_x_min.append(
                        mappings["avg"][scaling][group][subgroup][strain][time]
                        - std_val
                    )
                    sorted_x_max.append(
                        mappings["avg"][scaling][group][subgroup][strain][time]
                        + std_val
                    )
                sorted_x, sorted_y = list(
                    mappings["avg"][scaling][group][subgroup][strain].values()
                ), list(mappings["avg"][scaling][group][subgroup][strain].keys())

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
                if show_significance:
                    for time in t_values:
                        p_adj = p_values_dict[subgroup][strain][time]
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
                if show_error_bars:
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

                y_offset -= 1.5

            ax.legend(fontsize=6, title_fontsize=7, labelspacing=0.2, frameon=False)
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            ax.set_xlim(xlim[0] - 3, xlim[1] + 3)
            ax.set_ylim(ylim[0] - 3, ylim[1] + 3)

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
        if not os.path.exists(f"{path_results}/{scaling}/subplots"):
            os.makedirs(f"{path_results}/{scaling}/subplots")
        if mapping_mode != "mds-like":
            plt.savefig(
                f"{path_results}/{scaling}/subplots/{mapping_mode}_{subgroup}_{scaling}.svg",
                format="svg",
            )
            plots.append(
                f"{path_results}/{scaling}/subplots/{mapping_mode}_{subgroup}_{scaling}.svg"
            )
        else:
            plt.savefig(
                f"{path_results}/{scaling}/subplots/{subgroup}_{scaling}.svg",
                format="svg",
            )
            plots.append(f"{path_results}/{scaling}/subplots/{subgroup}_{scaling}.svg")
        plt.close()

    figures = [SVG(plot) for plot in plots]
    for i, figure in enumerate(figures):
        figure.move(0, i * 288)
    concatenated_figure = Figure("648pt", f"{288*len(figures)}pt", *figures)
    if mapping_mode != "mds-like":
        concatenated_svg_path = (
            f"{path_results}/{scaling}/{mapping_mode}_{scaling}_groups_{group}.svg"
        )
    else:
        concatenated_svg_path = f"{path_results}/{scaling}/{scaling}_groups_{group}.svg"
    concatenated_figure.save(concatenated_svg_path)
    if mapping_mode != "mds-like":
        pdf_path = (
            f"{path_results}/{scaling}/{mapping_mode}_{scaling}_groups_{group}.pdf"
        )
    else:
        pdf_path = f"{path_results}/{scaling}/{scaling}_groups_{group}.pdf"
    cairosvg.svg2pdf(url=concatenated_svg_path, write_to=pdf_path)


def get_stars(p):
    """
    Return a significance marker based on the p-value.

    Parameters:
    p (float): The adjusted p-value.

    Returns:
    str: A string of asterisks indicating significance.
    """
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
