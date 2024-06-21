"""
This script preprocesses the differential expression data for the top genes or gene sets based on a p-value threshold. The preprocessed data is saved as a compressed file for later use in the network visualization script.
"""

import pandas as pd
import re
import gzip
from scipy.stats import ttest_ind
import msgpack
import numpy as np
from differential_expression.code._constants import *


annotations_df = pd.read_csv(annotations_path)
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
combined_df.to_csv("../data/combined_df.csv")
gene_to_entrez = pd.Series(
    annotations_df["Entrez ID"].values, index=annotations_df["Gene"]
).to_dict()
combined_df_entrez = combined_df.copy()
combined_df_entrez.columns = [
    gene_to_entrez.get(x, x) for x in combined_df_entrez.columns
]
non_gene_columns = ["replicate", "hour", "strain"]
combined_df_entrez.columns = [
    gene_to_entrez.get(x, x) if x not in non_gene_columns else x
    for x in combined_df_entrez.columns
]
combined_df_entrez.to_csv("../data/combined_df_entrez.csv", index=False)

otherinfo_df = pd.read_csv(otherinfo_path).set_index("Gene")
annotations_df = annotations_df.set_index("Gene")


def merge_annotations(results_df, annotations_df):
    return results_df.join(annotations_df, on="Gene")


def get_top_genes(p, strain, time_point):
    strain_data = combined_df[combined_df["strain"] == strain]
    time_point_data = strain_data[strain_data["hour"] == time_point].drop(
        columns=["strain", "hour", "replicate"]
    )
    ax4_data = combined_df[combined_df["strain"] == "AX4"]
    ax4_time_point_data = ax4_data[ax4_data["hour"] == time_point].drop(
        columns=["strain", "hour", "replicate"]
    )

    t_stats, p_values = ttest_ind(
        time_point_data, ax4_time_point_data, equal_var=False, axis=0
    )
    fold_change = time_point_data.mean() / ax4_time_point_data.mean()
    results_df = pd.DataFrame(
        {
            "Gene": time_point_data.columns,
            "P-Value": p_values,
            "Fold Change": fold_change,
        }
    )
    results_df["Regulation"] = np.where(
        results_df["Fold Change"] > 1, "Upregulated", "Downregulated"
    )
    results_df = results_df[results_df["P-Value"] < p].sort_values(by="P-Value")
    results_df["Rank"] = results_df["P-Value"].rank(method="min").astype(int)
    results_df["P-Value"] = results_df["P-Value"].round(7)
    results_df = merge_annotations(results_df, annotations_df)
    results_df = merge_annotations(results_df, otherinfo_df)
    results_df["Symbol"] = results_df["Symbol"].apply(
        lambda x: "" if len(str(x)) > 10 else x
    )
    return results_df


def preprocess_data(p_threshold):
    strains = [strain for strain in metadata["strain"].unique() if strain != "AX4"]
    hours = [int(hour) for hour in sorted(metadata["hour"].unique())]
    results = {
        strain: {
            str(hour): get_top_genes(p_threshold, strain, hour).to_dict(
                orient="records"
            )
            for hour in hours
        }
        for strain in strains
    }
    save_data = input("Save data? (y/n) ")
    if save_data == "y":
        with gzip.open("../data/pairwise/preprocessed_data.msgpack.gz", "wb") as f:
            msgpack.dump(results, f, use_bin_type=True)


def filter_by_p_value(data, p_threshold):
    """Filter the preprocessed data by a new p-value threshold."""
    filtered_data = {}
    for strain, hours_data in data.items():
        filtered_data[strain] = {}
        for hour, genes in hours_data.items():
            filtered_data[strain][hour] = [
                gene for gene in genes if gene["P-Value"] < p_threshold
            ]
    return filtered_data


p_threshold = 1
preprocess_data(p_threshold)
