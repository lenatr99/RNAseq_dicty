"""
This script generates network data for different gene sets and community modes. It calculates the Jaccard similarity between gene sets and creates a graph based on the similarity matrix. The script then determines a threshold for adding edges to the graph and generates communities using the specified community mode (infomap or louvain). It also extracts the largest connected component of the graph and generates communities for this subgraph. Finally, it identifies the k-core communities and saves the resulting data for visualization.
"""

import pandas as pd
import numpy as np
import networkx as nx
from cdlib.algorithms import infomap, louvain
import json
import pickle
import os
from _constants import *


COMMUNITY_MODES = ["infomap", "louvain"]
community_mode = COMMUNITY_MODES[
    int(input(f"Enter community mode infomap - 0, louvain - 1: "))
]
if community_mode not in COMMUNITY_MODES:
    raise ValueError("Invalid community mode")
DATA_PATH = f"../results/time_series/communities/{community_mode}/"
os.makedirs(DATA_PATH, exist_ok=True)


layout_seed = 3


def calculate_jaccard(set1, set2):
    intersection = len(set(set1) & set(set2))
    union = len(set(set1) | set(set2))
    return intersection / union


def generate_graph_data(strain, GENE_SET, TERM):
    # Load GO terms data
    PATH = f"../results/time_series/gsea/json/{TERM}_{strain}.json"
    with open(PATH, "r") as f:
        go_terms_data = json.load(f)

    # Extract gene sets
    go_terms_list = [entry["Description"] for entry in go_terms_data]
    gene_sets = {}
    with open(f"../data/genesets/{GENE_SET}_mod.gmt", "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            go_term = parts[1][1:-1]
            genes = parts[2:]
            if go_term in go_terms_list:
                gene_sets[go_term] = genes

    # Create similarity matrix
    similarity_matrix = pd.DataFrame(
        index=gene_sets.keys(), columns=gene_sets.keys(), dtype=float
    )
    for set_name1, genes1 in gene_sets.items():
        for set_name2, genes2 in gene_sets.items():
            similarity = calculate_jaccard(genes1, genes2)
            similarity_matrix.loc[set_name1, set_name2] = similarity

    # Create graph
    G = nx.Graph()
    for set_name in gene_sets.keys():
        G.add_node(set_name)

    # Determine threshold for adding edges
    average_degree = 8
    N_nodes = len(G.nodes())
    N_links = N_nodes * average_degree / 2
    threshold = np.quantile(
        similarity_matrix.values.flatten(), 1 - N_links / (N_nodes * (N_nodes - 1) / 2)
    )

    # Add edges based on threshold
    for set_name1 in gene_sets.keys():
        for set_name2 in gene_sets.keys():
            if (
                set_name1 != set_name2
                and similarity_matrix.loc[set_name1, set_name2] > threshold
            ):
                G.add_edge(
                    set_name1,
                    set_name2,
                    weight=similarity_matrix.loc[set_name1, set_name2],
                )

    if community_mode == "infomap":
        communities = infomap(G).communities
    elif community_mode == "louvain":
        communities = louvain(G, weight="weight").communities

    graph_layout = nx.spring_layout(G, seed=layout_seed)

    G_cc = G.subgraph(max(nx.connected_components(G), key=lambda x: len(x)))

    if community_mode == "infomap":
        communities_cc = infomap(G_cc).communities
    elif community_mode == "louvain":
        communities_cc = louvain(G_cc, weight="weight").communities

    communities_kcore = {}
    communities_kcore_desc = {}

    k = 0
    for com in communities_cc:
        if k >= 4:
            break
        nodelist = [x for x in com if x in G_cc.nodes()]
        G_ = nx.k_core(G_cc.subgraph(nodelist), 5)
        if len(G_.nodes()) != 0:
            communities_kcore[k] = list(G_.nodes())
            communities_kcore_desc[k] = "summary"
            k += 1

    return (
        G,
        communities,
        G_cc,
        communities_cc,
        communities_kcore,
        communities_kcore_desc,
        graph_layout,
    )


all_graph_data = {}
all_communities = {}
all_graph_data_cc = {}
all_communities_cc = {}
all_communities_kcore = {}
all_communities_desc = {}
all_graph_layout = {}

for strain in strains:
    all_communities_kcore[strain] = {}
    all_communities_desc[strain] = {}
    for INDEX in range(len(SETS)):
        GENE_SET = SETS[INDEX]
        TERM = TERMS[INDEX]
        (
            graph,
            communities,
            graph_cc,
            communities_cc,
            communities_kcore,
            communities_kcore_desc,
            graph_layout,
        ) = generate_graph_data(strain, GENE_SET, TERM)
        all_graph_data[(strain, GENE_SET, TERM)] = graph
        all_communities[(strain, GENE_SET, TERM)] = communities
        all_graph_data_cc[(strain, GENE_SET, TERM)] = graph_cc
        all_communities_cc[(strain, GENE_SET, TERM)] = communities_cc
        all_communities_kcore[strain][TERM] = communities_kcore
        all_communities_desc[strain][TERM] = communities_kcore_desc
        all_graph_layout[(strain, GENE_SET, TERM)] = graph_layout


# Make directories if they don't exist
os.makedirs(DATA_PATH + "pkl", exist_ok=True)
os.makedirs(DATA_PATH + "json", exist_ok=True)

# Saving the graph data
with open(DATA_PATH + f"pkl/graph_data.pkl", "wb") as f:
    pickle.dump(all_graph_data, f)

with open(DATA_PATH + f"pkl/graph_data_cc.pkl", "wb") as f:
    pickle.dump(all_graph_data_cc, f)

with open(DATA_PATH + f"pkl/communities.pkl", "wb") as f:
    pickle.dump(all_communities, f)

with open(DATA_PATH + f"pkl/communities_cc.pkl", "wb") as f:
    pickle.dump(all_communities_cc, f)

with open(DATA_PATH + f"pkl/graph_layout.pkl", "wb") as f:
    pickle.dump(all_graph_layout, f)

with open(DATA_PATH + f"json/communities.json", "w") as f:
    json.dump(all_communities_kcore, f, indent=4)

with open(DATA_PATH + f"json/communities_desc_template.json", "w") as f:
    json.dump(all_communities_desc, f, indent=4)
