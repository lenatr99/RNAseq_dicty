import pandas as pd
from flask import Flask, render_template
from flask import request
import json
import gzip
import msgpack
from flask import send_file
import os
import math

def round_to_n_significant_digits(x, n):
    if x == 0:
        return 0
    else:
        from math import log10, floor

        return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def apply_color(row, gene_sets_to_color):
    colors = {
        "0": "background-color: #CAE1FF",  # Light blue
        "1": "background-color: #FFC0CB",  # Light red
        "2": "background-color: #90EE90",  # Light green
        "3": "background-color: #FFD700",  # Gold
    }
    default_style = [""] * len(row)

    for key, descriptions in gene_sets_to_color.items():
        if row["Description"] in descriptions:
            return [colors[key]] * len(row)
    return default_style

def clean_nan_values(data):
    """Recursively clean NaN values and replace them with None."""
    if isinstance(data, dict):
        return {k: clean_nan_values(v) for k, v in data.items()}
    elif isinstance(data, list):
        return [clean_nan_values(i) for i in data]
    elif isinstance(data, float) and math.isnan(data):
        return ""
    return data

app = Flask(
    __name__,
    static_folder="/Users/lenatrnovec/FAKS/BIOLAB/RNAseq_dicty/differential_expression/templates/static",  # Change this to your static folder
)


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


def filter_by_p_value2(data, p_threshold2):
    """Filter the preprocessed data by a new p-value threshold."""
    filtered_data = {}
    for strain, genes in data.items():
        filtered_data[strain] = [
            gene for gene in genes if gene["adj.P.Val"] < p_threshold2
        ]
    return filtered_data


with open("results/pairwise/communities/infomap/json/communities.json", "r") as file:
    gene_sets_to_color_infomap = json.load(file)

with open("results/pairwise/communities/louvain/json/communities.json", "r") as file:
    gene_sets_to_color_louvain = json.load(file)

with open("results/time_series/communities/infomap/json/communities.json", "r") as file:
    gene_sets_to_color_infomap_gs = json.load(file)

with open("results/pairwise/communities/louvain/json/communities.json", "r") as file:
    gene_sets_to_color_louvain_gs = json.load(file)

strain_names = {
    "B1": "tgrB1<sup>-</sup>",
    "C1": "tgrC1<sup>-</sup>",
    "rgB": "rapgapB<sup>-</sup>",
    "B1_rgB": "tgrB1<sup>-</sup> rapgapB<sup>-</sup>",
    "C1_rgB": "tgrC1<sup>-</sup> rapgapB<sup>-</sup>",
    "AX4L846F": "AX4 L846F<sup> </sup>",
    "B1_L846F": "tgrB1<sup>-</sup> L846F",
}


@app.route("/", methods=["GET", "POST"])
def index():
    p_threshold = 0.01

    if request.method == "POST":
        p_threshold = float(request.form["p_value"])
    with gzip.open("results/pairwise/preprocessed_data.msgpack.gz", "rb") as f:
        complete_data = msgpack.load(f, raw=False)
    filtered_data = filter_by_p_value(complete_data, p_threshold)
    gene_set_data_bp = {}
    gene_set_data_mf = {}
    gene_set_data_db = {}
    for strain in filtered_data:
        gene_set_data_bp[strain] = {}
        gene_set_data_mf[strain] = {}
        gene_set_data_db[strain] = {}
        for hour in filtered_data[strain]:
            with open(
                f"results/pairwise/gsea/BP/DE_{hour}_{strain}.json", "r"
            ) as f:
                gene_set_data_bp[strain][hour] = json.load(f)
            with open(
                f"results/pairwise/gsea/MF/DE_{hour}_{strain}.json", "r"
            ) as f:
                gene_set_data_mf[strain][hour] = json.load(f)
            with open(
                f"results/pairwise/gsea/DB/DE_{hour}_{strain}.json", "r"
            ) as f:
                gene_set_data_db[strain][hour] = json.load(f)
    strains = list(filtered_data.keys())
    hours = sorted(filtered_data[strains[0]], key=int)
    filtered_data = clean_nan_values(filtered_data)
    return render_template(
        "index.html",
        data=filtered_data,
        strain_names=strain_names,
        gene_set_data_bp=gene_set_data_bp,
        gene_set_data_mf=gene_set_data_mf,
        gene_set_data_db=gene_set_data_db,
        strains=strains,
        hours=hours,
        p_threshold=p_threshold,
        gene_sets_to_color_infomap=gene_sets_to_color_infomap,
        gene_sets_to_color_louvain=gene_sets_to_color_louvain,
    )


@app.route("/time-series", methods=["GET", "POST"])
def time_series():
    p_threshold2 = 0.0001

    if request.method == "POST":
        p_threshold2 = float(request.form["p_value"])
    with gzip.open("results/time_series/preprocessed_data.msgpack.gz", "rb") as f:
        complete_data = msgpack.load(f, raw=False)
    for strain, genes in complete_data.items():
        for gene in genes:
            gene["adj.P.Val"] = round_to_n_significant_digits(gene["adj.P.Val"], 3)
    filtered_data = filter_by_p_value2(complete_data, p_threshold2)
    gene_set_data_bp = {}
    gene_set_data_mf = {}
    gene_set_data_db = {}
    for strain in filtered_data:
        with open(f"results/time_series/gsea/json/BP_{strain}.json", "r") as f:
            gene_set_data_bp[strain] = json.load(f)
            for gene_set in gene_set_data_bp[strain]:
                gene_set["adj.pValue"] = round_to_n_significant_digits(
                    gene_set["adj.pValue"], 3
                )
        with open(f"results/time_series/gsea/json/MF_{strain}.json", "r") as f:
            gene_set_data_mf[strain] = json.load(f)
            for gene_set in gene_set_data_mf[strain]:
                gene_set["adj.pValue"] = round_to_n_significant_digits(
                    gene_set["adj.pValue"], 3
                )
        with open(f"results/time_series/gsea/json/DB_{strain}.json", "r") as f:
            gene_set_data_db[strain] = json.load(f)
            for gene_set in gene_set_data_db[strain]:
                gene_set["adj.pValue"] = round_to_n_significant_digits(
                    gene_set["adj.pValue"], 3
                )
    strains = list(filtered_data.keys())
    filtered_data = clean_nan_values(filtered_data)
    return render_template(
        "time-series.html",
        data=filtered_data,
        strain_names=strain_names,
        strains=strains,
        p_threshold2=p_threshold2,
        gene_set_data_bp=gene_set_data_bp,
        gene_set_data_mf=gene_set_data_mf,
        gene_set_data_db=gene_set_data_db,
        gene_sets_to_color_infomap=gene_sets_to_color_infomap_gs,
        gene_sets_to_color_louvain=gene_sets_to_color_louvain_gs,
    )


@app.route("/download-excel/<strain>/<p_threshold2>")
def download_excel(strain, p_threshold2=0.0001):
    p_threshold2 = float(p_threshold2)
    with gzip.open("results/time_series/preprocessed_data.msgpack.gz", "rb") as f:
        complete_data = msgpack.load(f, raw=False)
    filtered_data = filter_by_p_value2(complete_data, p_threshold2)
    data = filtered_data[strain]
    df = pd.DataFrame(data)
    filename = f"time_series_{strain}_{p_threshold2}.xlsx"
    filepath = os.path.join("temporary", filename)
    os.makedirs("temporary", exist_ok=True)
    df.to_excel(filepath, index=False)

    return send_file(filepath, as_attachment=True, download_name=filename)


@app.route("/download-all-excel")
def download_all_excel():
    p_threshold2 = 0.0001
    with gzip.open("results/time_series/preprocessed_data.msgpack.gz", "rb") as f:
        complete_data = msgpack.load(f, raw=False)
    filtered_data = filter_by_p_value2(complete_data, p_threshold2)
    writer = pd.ExcelWriter("all_data.xlsx", engine="openpyxl")
    for strain, hours_data in filtered_data.items():
        for hour, data in hours_data.items():
            df = pd.DataFrame(data)
            df.to_excel(writer, sheet_name=f"{strain}_{hour}", index=False)
    writer.save()
    writer.close()

    return send_file(
        "all_data.xlsx", as_attachment=True, attachment_filename="all_data.xlsx"
    )


@app.route("/download-gene-sets/<strain>/<set_type>/<clustering_method>")
def download_gene_sets(strain, set_type, clustering_method):
    gene_set_data_bp = {}
    gene_set_data_mf = {}
    gene_set_data_db = {}
    with open(f"results/time_series/gsea/json/BP_{strain}.json", "r") as f:
        gene_set_data_bp[strain] = json.load(f)
    with open(f"results/time_series/gsea/json/MF_{strain}.json", "r") as f:
        gene_set_data_mf[strain] = json.load(f)
    with open(f"results/time_series/gsea/json/DB_{strain}.json", "r") as f:
        gene_set_data_db[strain] = json.load(f)
    if set_type == "bp":
        gene_set_data = gene_set_data_bp
    elif set_type == "mf":
        gene_set_data = gene_set_data_mf
    elif set_type == "db":
        gene_set_data = gene_set_data_db
    else:
        return "Invalid gene set type", 404
    if clustering_method not in ["infomap", "louvain"]:
        return "Invalid clustering method", 404

    color_mapping = (
        gene_sets_to_color_infomap_gs
        if clustering_method == "infomap"
        else gene_sets_to_color_louvain_gs
    )

    data = gene_set_data[strain]
    df = pd.DataFrame(data)

    set_type = set_type.upper()

    df = df.style.apply(
        lambda x: apply_color(x, color_mapping[strain][set_type]), axis=1
    )

    filename = f"gene_sets_{set_type}_{clustering_method}_{strain}.xlsx"
    filepath = os.path.join("temporary", filename)
    os.makedirs("temporary", exist_ok=True)
    df.to_excel(filepath, engine="openpyxl", index=False)

    return send_file(filepath, as_attachment=True, download_name=filename)


if __name__ == "__main__":
    app.run(debug=True)
