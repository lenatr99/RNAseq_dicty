"""
This script preprocesses the differential expression data for the top genes or gene sets based on a p-value threshold. The preprocessed data can be saved as a compressed file for later use in the network visualization script.
"""

import pandas as pd
import re
import msgpack
import gzip
import os
from CONSTANTS import *

GENESETS = input("Do you want to preprocess genes [0] or gene sets [1]? ")
if GENESETS not in ["0", "1"]:
    raise ValueError("Invalid input")

if GENESETS == "1":
    GENESETS = True
else:
    GENESETS = False


def merge_annotations(results_df, annotations_df):
    return results_df.join(annotations_df, on="Gene")

def round_to_n_significant_digits(x, n):
    if x == 0:
        return 0
    else:
        from math import log10, floor

        return round(x, -int(floor(log10(abs(x)))) + (n - 1))


def get_top_genes(p, strain):
    de_results_path = f"../data/time_series/genes/DE_{strain}.csv"
    strain_data = pd.read_csv(de_results_path)
    filtered_data = strain_data[strain_data["adj.P.Val"] < p].sort_values(
        by="adj.P.Val"
    )
    filtered_data["Rank"] = filtered_data["adj.P.Val"].rank(method="min").astype(int)
    filtered_data = filtered_data.rename(columns={"row_IDs": "Gene"})
    filtered_data = filtered_data.drop(
        columns=["b_0", "b_1", "b_2", "b_3", "d_0", "d_1", "d_2", "d_3"]
    )
    filtered_data = merge_annotations(filtered_data, annotations_df)
    filtered_data = merge_annotations(filtered_data, otherinfo_df)
    if "X" in filtered_data.columns:
        del filtered_data["X"]
    filtered_data["Symbol"] = filtered_data["Symbol"].apply(
        lambda x: "" if len(str(x)) > 10 else x
    )
    return filtered_data


def get_top_gene_sets(p, strain):
    de_results_path = f"../data/time_series/gsea/DE_{strain}_{TERM}_gsea.csv"
    strain_data = pd.read_csv(de_results_path).set_index("pathway")
    filtered_data = strain_data[strain_data["adj.pValue"] < p].sort_values(
        by="adj.pValue"
    )
    filtered_data["Rank"] = filtered_data["adj.pValue"].rank(method="min").astype(int)
    filtered_data = filtered_data.rename(columns={"description": "Description"})
    filtered_data["Description"] = filtered_data["Description"].str[1:-1]
    filtered_data = filtered_data.rename(columns={"description": "Term"})
    filtered_data = filtered_data.drop(columns=["Unnamed: 0", "overlap"])
    # save only the top 100 gene sets
    filtered_data = filtered_data.iloc[:100]
    if not os.path.exists("../data/time_series/gsea/json"):
        os.makedirs("../data/time_series/gsea/json")
    filtered_data.to_json(
        f"../data/time_series/gsea/json/{TERM}_{strain}.json", orient="records", indent=4
    )
    return filtered_data


def preprocess_data(p_threshold, genesets=False, term=""):
    strains = [strain for strain in metadata["strain"].unique() if strain != "AX4"]
    if genesets:
        results = {
            strain: get_top_gene_sets(p_threshold, strain).to_dict(orient="records")
            for strain in strains
        }
    else:
        results = {
            strain: get_top_genes(p_threshold, strain).to_dict(orient="records")
            for strain in strains
        }
    if genesets:
        save_data = input(f"Save compressed data for gene set for {term}? (y/n) ")
    else:
        save_data = input(f"Save compressed data? (y/n) ")
    if save_data not in ["y", "n"]:
        raise ValueError("Invalid input")
    if save_data == "y" and genesets == False:
        with gzip.open("../data/time_series/preprocessed_data.msgpack.gz", "wb") as f:
            msgpack.dump(results, f, use_bin_type=True)
    elif save_data == "y" and genesets == True:
        with gzip.open(
            f"../data/time_series/preprocessed_data_genesets_{TERM}.msgpack.gz", "wb"
        ) as f:
            msgpack.dump(results, f, use_bin_type=True)


def filter_by_p_value(data, p_threshold):
    """Filter the preprocessed data by a new p-value threshold."""
    filtered_data = {}
    for strain, hours_data in data.items():
        filtered_data[strain] = {}
        for hour, genes in hours_data.items():
            filtered_data[strain][hour] = [
                gene for gene in genes if gene["adj.P.Val"] < p_threshold
            ]
    return filtered_data

if GENESETS:
    for TERM in TERMS:
        annotations_df = pd.read_csv(annotations_path).set_index("Gene")
        otherinfo_df = pd.read_csv(otherinfo_path).set_index("Gene")
        df = pd.read_csv(data_path).set_index("ddb_g").transpose()
        df.index = df.index.astype(str)

        metadata = df.index.to_series().apply(
            lambda s: pd.Series(
                {
                    "replicate": (
                        int(re.search(r"_r(\d+)", s).group(1))
                        if re.search(r"_r(\d+)", s)
                        else None
                    ),
                    "hour": (
                        int(re.search(r"(\d+)h", s).group(1))
                        if re.search(r"(\d+)h", s)
                        else None
                    ),
                    "strain": (
                        re.search(r"^(.*?)_(\d+h|r\d+)", s).group(1).replace(".", "_")
                        if re.search(r"^(.*?)_(\d+h|r\d+)", s)
                        else None
                    ),
                }
            )
        )

        combined_df = pd.concat([metadata, df], axis=1)

        p_threshold = 1
        preprocess_data(p_threshold, genesets=GENESETS, term=TERM)

else:
    annotations_df = pd.read_csv(annotations_path).set_index("Gene")
    otherinfo_df = pd.read_csv(otherinfo_path).set_index("Gene")
    df = pd.read_csv(data_path).set_index("ddb_g").transpose()
    df.index = df.index.astype(str)

    metadata = df.index.to_series().apply(
        lambda s: pd.Series(
            {
                "replicate": (
                    int(re.search(r"_r(\d+)", s).group(1))
                    if re.search(r"_r(\d+)", s)
                    else None
                ),
                "hour": (
                    int(re.search(r"(\d+)h", s).group(1))
                    if re.search(r"(\d+)h", s)
                    else None
                ),
                "strain": (
                    re.search(r"^(.*?)_(\d+h|r\d+)", s).group(1).replace(".", "_")
                    if re.search(r"^(.*?)_(\d+h|r\d+)", s)
                    else None
                ),
            }
        )
    )

    combined_df = pd.concat([metadata, df], axis=1)

    p_threshold = 1
    preprocess_data(p_threshold, genesets=GENESETS)