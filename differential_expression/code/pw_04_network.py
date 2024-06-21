"""
Before running this script, make sure to use ChatGPT or another tool to generate the common function of the gene sets in the communities. The template
for that is saved in the communities_desc_template.json file. Save the output of the ChatGPT model in the communities_desc.json file.

This script visualizes network data in the context of different community modes and gene sets. It allows the user to choose a community mode 
(infomap or louvain), then loads corresponding network data, communities, and layout information. The script generates plots for each hour 
and gene set, displaying the network with color-coded nodes representing different communities. It includes interactive elements, such as 
draggable community descriptions and a button to save the plot. The visualization is tailored for exploring biological process (BP), molecular 
function (MF), and Dictybase Phenotypes (DB) within a specified dataset.
"""

import matplotlib.pyplot as plt
import pickle
import json
from matplotlib.widgets import Button
import matplotlib.patheffects as patheffects
import networkx as nx
import numpy as np
import os
from _constants import *

COMMUNITY_MODES = ["infomap", "louvain"]
community_mode = COMMUNITY_MODES[
    int(input(f"Enter community mode infomap - 0, louvain - 1: "))
]
DATA_PATH = f"../results/pairwise/communities/{community_mode}/"

# Check if json/communities_desc.json exists
if not os.path.exists(DATA_PATH + "json/communities_desc.json"):
    raise FileNotFoundError(
        "json/communities_desc.json not found. Before continuing, make sure to use ChatGPT or another tool to generate the common function of the gene sets in the communities. The template for that is saved in the communities_desc_template.json file. Save the output of the ChatGPT model in the communities_desc.json file in the same folder."
    )

strain = strains[
    int(
        input(
            f"Enter strain B1 - 0, C1 - 1, rgB - 2, B1_rgB - 3, C1_rgB - 4, AX4L846F - 5, B1_L846F - 6: "
        )
    )
]

color_list = ["royalblue", "red", "forestgreen", "orange", "k"]
N_communities = len(color_list)
NODE_SIZE = 15


def insert_line_breaks(text, max_line_length):
    words = text.split()
    lines = []
    current_line = []
    for word in words:
        if (
            sum(len(w) for w in current_line) + len(current_line) - 1 + len(word)
            > max_line_length
        ):
            lines.append(" ".join(current_line))
            current_line = [word]
        else:
            current_line.append(word)
    if current_line:
        lines.append(" ".join(current_line))
    return "\n".join(lines)


def xy_for_community_text(community_nodes, graph_layout, occupied_positions, delta=0.1):
    x, y = zip(*[graph_layout[node] for node in community_nodes])
    centroid = (np.mean(x), np.mean(y))

    best_pos = centroid
    min_distance = float("inf")
    for dx in np.arange(-1, 1, delta):
        for dy in np.arange(-1, 1, delta):
            test_pos = (centroid[0] + dx, centroid[1] + dy)
            distance = min(
                np.linalg.norm(np.array(test_pos) - np.array(pos))
                for pos in occupied_positions
            )
            if distance > min_distance:
                min_distance = distance
                best_pos = test_pos
    h_align = "center"
    if best_pos[0] < centroid[0]:
        h_align = "right"
    elif best_pos[0] > centroid[0]:
        h_align = "left"

    return best_pos, h_align


def save_plot(event):
    button_ax.set_visible(False)
    plt.draw()
    name = f"{community_mode}_{TERM}_{strain}_{hour}.png"
    plt.savefig(
        f"../templates/static/{community_mode}/{name}", bbox_inches="tight", dpi=300
    )
    print(f"Plot saved as /templates/static/{community_mode}/{name}!")
    button_ax.set_visible(True)
    plt.draw()


def calculate_button_position(fig, button_width=0.1, button_height=0.05, offset=0.02):
    bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig_width, fig_height = bbox.width, bbox.height
    left = offset
    bottom = offset
    return [left, bottom, button_width, button_height]


def load_data(file_path):
    with open(DATA_PATH + file_path, "rb" if file_path.endswith(".pkl") else "r") as f:
        return pickle.load(f) if file_path.endswith(".pkl") else json.load(f)


# CLASSES


class DraggableText:
    def __init__(self, text):
        self.text = text
        self.text.figure.canvas.mpl_connect("button_press_event", self.on_click)
        self.text.figure.canvas.mpl_connect("button_release_event", self.on_release)
        self.text.figure.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.is_dragging = False

    def on_click(self, event):
        if event.inaxes != self.text.axes:
            return
        contains, _ = self.text.contains(event)
        if contains:
            self.is_dragging = True
            self.mouse_x = event.xdata
            self.mouse_y = event.ydata

    def on_motion(self, event):
        if not self.is_dragging or not event.inaxes:
            return
        dx, dy = event.xdata - self.mouse_x, event.ydata - self.mouse_y

        # Adjust these scaling factors if the movement is still too fast or slow
        scale_factor_x = 0.5
        scale_factor_y = 0.5

        new_x = self.text.get_position()[0] + dx * scale_factor_x
        new_y = self.text.get_position()[1] + dy * scale_factor_y

        self.text.set_x(new_x)
        self.text.set_y(new_y)

        self.mouse_x, self.mouse_y = event.xdata, event.ydata
        self.text.figure.canvas.draw()

    def on_release(self, event):
        self.is_dragging = False


all_graph_data = load_data("pkl/graph_data.pkl")
all_communities = load_data("pkl/communities.pkl")
all_graph_data_cc = load_data("pkl/graph_data_cc.pkl")
all_communities_cc = load_data("pkl/communities_cc.pkl")
all_graph_layout = load_data("pkl/graph_layout.pkl")
communities_desc = load_data("json/communities_desc.json")

for hour in hours:
    for i in range(len(SETS)):
        hour, GENE_SET, TERM = hour, SETS[i], TERMS[i]
        graph = all_graph_data[(strain, hour, GENE_SET, TERM)]
        communities = all_communities[(strain, hour, GENE_SET, TERM)]
        graph_cc = all_graph_data_cc[(strain, hour, GENE_SET, TERM)]
        communities_cc = all_communities_cc[(strain, hour, GENE_SET, TERM)]
        graph_layout = all_graph_layout[(strain, hour, GENE_SET, TERM)]

        com_names_lgg = communities_desc[strain][TERM][str(hour)]
        f, ax = plt.subplots(1, 1, figsize=(7, 5))
        nx.draw_networkx(
            graph_cc,
            graph_layout,
            with_labels=False,
            node_color="lightgray",
            edge_color="lightgray",
            node_size=0,
            alpha=0.3,
        )
        for loc in ["right", "left", "top", "bottom"]:
            ax.spines[loc].set_visible(False)

        k = 0
        for com in communities_cc:
            if k >= 4:
                break
            nodelist = [x for x in com if x in graph_cc.nodes()]
            G_ = nx.k_core(graph_cc.subgraph(nodelist), 5)
            if len(G_.nodes()) != 0:
                nx.draw_networkx_nodes(
                    graph_cc,
                    graph_layout,
                    node_color=color_list[k],
                    nodelist=nodelist,
                    node_size=NODE_SIZE,
                    ax=ax,
                    alpha=0.0,
                )
                for name in G_.nodes():
                    nx.draw_networkx_nodes(
                        G_,
                        graph_layout,
                        node_color=color_list[k],
                        node_size=NODE_SIZE,
                        ax=ax,
                        alpha=1,
                        nodelist=[name],
                    )
                k += 1

        interesting_coms = np.arange(min(4, len(communities_cc)))

        max_line_length = 30
        for key in com_names_lgg:
            com_names_lgg[key] = insert_line_breaks(com_names_lgg[key], max_line_length)

        occupied_positions = []
        for node in graph_layout:
            occupied_positions.append(graph_layout[node])

        draggable_texts = []
        for i in range(len(com_names_lgg)):
            community_nodes = communities_cc[i]
            text_pos, h_align = xy_for_community_text(
                community_nodes, graph_layout, occupied_positions
            )
            text = ax.text(
                text_pos[0],
                text_pos[1],
                com_names_lgg[str(i)],
                color=color_list[i],
                fontsize=6,
                horizontalalignment=h_align,
            )
            text.set_path_effects(
                [patheffects.withStroke(linewidth=2, foreground="white", alpha=0.7)]
            )
            draggable_texts.append(DraggableText(text))
            occupied_positions.append(text_pos)

        button_pos = calculate_button_position(f)
        button_ax = plt.axes(button_pos)
        print(
            f"Now looking at hour {hour} and gene set {GENE_SET} for strain {strain} and term {TERM}"
        )
        button = Button(button_ax, "Save Plot")
        button.on_clicked(save_plot)

        plt.show()
