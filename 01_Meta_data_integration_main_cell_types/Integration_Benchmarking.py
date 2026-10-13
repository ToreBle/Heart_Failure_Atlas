import scanpy as sc
from harmony import harmonize
import gc
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
import scib


adata = sc.read_h5ad("snRNA_seq_human_heart_atlas.h5ad")

cell_type_map = {
    # Adipocytes
    "Adipocyte": "Adipocytes",
    "Adipocytes": "Adipocytes",
    "fat cell": "Adipocytes",

    # Cardiomyocytes
    "Cardiomyocyte": "Cardiomyocytes",
    "Cardiomyocyte_I": "Cardiomyocytes",
    "Cardiomyocyte_II": "Cardiomyocytes",
    "Cardiomyocyte_III": "Cardiomyocytes",
    "Ventricular_Cardiomyocyte": "Cardiomyocytes",
    "cardiac muscle cell": "Cardiomyocytes",

    # Cycling cells
    "Cycling cells": "Cycling cells",

    # Endocardial
    "Endocardial": "Endocardial",
    "Endocardium": "Endocardial",

    # Endothelial
    "Endothelial": "Endothelial",
    "Endothelial_I": "Endothelial",
    "Endothelial_II": "Endothelial",
    "Endothelial_III": "Endothelial",
    "Endothelium": "Endothelial",
    "cardiac endothelial cell": "Endothelial",
    "endothelial cell": "Endothelial",

    # Epicardium
    "Epicardial": "Epicardium",
    "Epicardium": "Epicardium",
    "mesothelial cell": "Epicardium",

    # Fibroblast
    "Activated_fibroblast": "Fibroblast",
    "Fibroblast": "Fibroblast",
    "Fibroblast_I": "Fibroblast",
    "Fibroblast_II": "Fibroblast",
    "Fibroblasts": "Fibroblast",
    "cardiac ventricle fibroblast": "Fibroblast",
    "fibroblast of cardiac tissue": "Fibroblast",

    # Lymphatic EC
    "Lymphatic": "Lymphatic EC",
    "Lymphatic_endothelial": "Lymphatic EC",
    "endothelial cell of lymphatic vessel": "Lymphatic EC",

    # Lymphocytes
    "B_Cells": "Lymphocytes",
    "Lymphocyte": "Lymphocytes",
    "Lymphoid": "Lymphocytes",
    "T/NK_Cells": "Lymphocytes",
    "lymphocyte": "Lymphocytes",

    # Mast cells
    "Mast": "Mast cells",
    "Mast_Cells": "Mast cells",
    "Mast_cell": "Mast cells",
    "mast cell": "Mast cells",

    # Myeloid
    "Macrophage": "Myeloid",
    "Macrophages": "Myeloid",
    "Monocytes": "Myeloid",
    "Myeloid": "Myeloid",
    "myeloid cell": "Myeloid",
    "macrophage": "Myeloid",

    # Neuronal
    "Neuronal": "Neuronal",
    "Neurons": "Neuronal",
    "cardiac neuron": "Neuronal",

    # Pericytes
    "Pericyte": "Pericytes",
    "Pericyte_I": "Pericytes",
    "Pericyte_II": "Pericytes",
    "Pericytes": "Pericytes",
    "mural cell": "Pericytes",
    "pericyte cell": "Pericytes",

    # VSMC
    "Smooth_Muscle": "VSMC",
    "Smooth_muscle_cells": "VSMC",
    "VSMC": "VSMC",
    "vSMCs": "VSMC",
    "vascular associated smooth muscle cell": "VSMC",

    # Ambigous
    "NotAssigned": "NotAssigned",
    "native cell": "NotAssigned",
}

adata.obs["cell_type_harmonized"] = (
    adata.obs["Cell_type"]
    .replace(cell_type_map)
)

sc.external.pp.harmony_integrate(adata, key="Donor", adjusted_basis="X_harmony")

sc.pp.neighbors(adata, use_rep='X_harmony')

sc.tl.leiden(
    adata,
    resolution=0.1,
    flavor="igraph",
    directed=False,
    key_added='harmony_leiden_01' 
)

sc.tl.leiden(
    adata,
    resolution=0.2,
    flavor="igraph",
    directed=False,
    key_added='harmony_leiden_02' 
)

sc.tl.leiden(
    adata,
    resolution=0.3,
    flavor="igraph",
    directed=False,
    key_added='harmony_leiden_03' 
)

#Harmony Annotation
cell_dict = {'Pericytes': ['0'], 'VSMC': ['1'], 'ischemic VSMC': ['2'], 'Fibroblast': ['3', '15'],
    'Cardiomyocytes': ['4'],'Lympathic EC': ['5'],'Lymphocytes': ['6'], 
             'Endothelial':['7'], 'Endocardial': ['8'], 'Myeloid': ['9'], "Mast cells":["10"],
            'Neuronal':['11'], 'Adipocytes':['12'], 'Epicardium': ['13'], 'ischemic EC':["14"]}
adata.obs['Annotation_Harmony'] = np.nan
for i in cell_dict.keys():
    ind = pd.Series(adata.obs.harmony_leiden_02).isin(cell_dict[i])
    adata.obs.loc[ind,'Annotation_Harmony'] = i

sns.set_style("whitegrid", {'axes.grid' : False})
plt.figure(figsize = (6,6))
Harmony_Xtable = pd.crosstab(adata.obs['Annotation_Harmony'], adata.obs['Annotation_1'])
Harmony_Xtable = Harmony_Xtable/Harmony_Xtable.sum()*100
sns.heatmap(Harmony_Xtable, xticklabels=True, yticklabels=True, square=True, cbar_kws={'label': 'Percentage of cells', "shrink": 0.4}, cmap="Blues")

adata_bbknn = adata.copy()
sc.external.pp.bbknn(adata_bbknn, batch_key="Donor")

sc.tl.leiden(
    adata_bbknn,
    resolution=0.1,
    flavor="igraph",
    directed=False,
    key_added='BBKNN_leiden_01' 
)

sc.tl.leiden(
    adata_bbknn,
    resolution=0.2,
    flavor="igraph",
    directed=False,
    n_iterations = 2,
    key_added='BBKNN_leiden_02' 
)

sc.tl.leiden(
    adata_bbknn,
    resolution=0.3,
    flavor="igraph",
    directed=False,
    key_added='BBKNN_leiden_03' 
)

#BBKNN Annotation
cell_dict = {'Pericytes': ['0'], 'VSMC': ['1'], 'Fibroblast': ['2'], 'Cardiomyocytes': ['3'],
                'Lympathic EC': ['4'], 'Endothelial':['5'], 'Epicardium':['6'], 'Myeloid': ['7'], 'Endocardial': ['8'], 'Lymphocytes': ['9'],
                "Mast cells":["10"], 'Neuronal':['11'], 'Adipocytes':['12']}
adata_bbknn.obs['Annotation_BBKNN'] = np.nan
for i in cell_dict.keys():
    ind = pd.Series(adata_bbknn.obs.BBKNN_leiden_03).isin(cell_dict[i])
    adata_bbknn.obs.loc[ind,'Annotation_BBKNN'] = i

sns.set_style("whitegrid", {'axes.grid' : False})
plt.figure(figsize = (6,6))
BBKNN_Xtable = pd.crosstab(adata_bbknn.obs['Annotation_BBKNN'], adata_bbknn.obs['Annotation_1'])
BBKNN_Xtable = BBKNN_Xtable/BBKNN_Xtable.sum()*100
sns.heatmap(BBKNN_Xtable, xticklabels=True, yticklabels=True, square=True, cbar_kws={'label': 'Percentage of cells', "shrink": 0.4}, cmap="Blues")


# --------------------------------------------------
# Benchmark BBKNN, scVI and Harmony
# --------------------------------------------------

methods = {
    "BBKNN": {
        "adata": adata_bbknn,
        "type": "knn",
        "cluster_key": "Annotation_BBKNN",
        "embed": None,
    },
    "scVI": {
        "adata": adata_scvi,
        "type": "embed",
        "cluster_key": "Annotation_1",
        "embed": "X_scVI",
    },
    "Harmony": {
        "adata": adata_harmony,
        "type": "embed",
        "cluster_key": "Annotation_Harmony",
        "embed": "X_harmony",
    },
}

results = []

for method_name, config in methods.items():

    adata_int = config["adata"]
    cluster_key = config["cluster_key"]

    # Check required annotations
    if cluster_key not in adata_int.obs:
        raise KeyError(
            f"{cluster_key} missing for {method_name}"
        )

    # Use the same reference cell types
    # for all integration methods
    label_key = "cell_type_harmonized"

    # Exclude unassigned cells from ARI/NMI
    valid = (
        adata_int.obs[cluster_key].notna()
        & adata_int.obs[label_key].notna()
    )

    adata_labels = adata_int[valid].copy()

    # ARI and NMI using method-specific annotations
    ari_score = scib.metrics.ari(
        adata_labels,
        cluster_key=cluster_key,
        label_key=label_key
    )

    nmi_score = scib.metrics.nmi(
        adata_labels,
        cluster_key=cluster_key,
        label_key=label_key
    )

    # Other integration metrics
    kwargs = dict(
        batch_key="Donor",
        label_key=label_key,
        type_=config["type"],

        graph_conn_=True,
        kBET_=False,
        ilisi_=True,
        clisi_=True,

        # ARI/NMI calculated separately
        nmi_=False,
        ari_=False,
        isolated_labels_f1_=False,

        silhouette_=False,
        pcr_=False,
        hvg_score_=False,
        cell_cycle_=False,
        trajectory_=False,

        n_cores=6,
    )

    if config["type"] == "embed":
        kwargs["embed"] = config["embed"]

    metrics_df = scib.metrics.metrics(
        adata,
        adata_int,
        **kwargs
    )

    # Extract scores from scIB result
    scores = metrics_df["score"]

    def get_metric(names):
        for name in names:
            if name in scores.index:
                return float(scores[name])
        raise KeyError(
            f"Metric not found. Available: "
            f"{list(scores.index)}"
        )

    ilisi_score = get_metric([
        "iLISI_score",
        "iLISI_score/batch",
        "iLISI"
    ])

    clisi_score = get_metric([
        "cLISI_score",
        "cLISI_score/label",
        "cLISI"
    ])

    # Collect all results
    for benchmark, value in {
        "ARI": ari_score,
        "NMI": nmi_score,
        "cLISI": clisi_score,
        "iLISI": ilisi_score,
    }.items():

        results.append({
            "Benchmark": benchmark,
            "Method": method_name,
            "Value": float(value)
        })

df = pd.DataFrame(results)

benchmark_order = ["ARI", "NMI", "cLISI", "iLISI"]
method_order = ["BBKNN", "scVI", "Harmony"]

df["Benchmark"] = pd.Categorical(
    df["Benchmark"],
    categories=benchmark_order,
    ordered=True
)

df["Method"] = pd.Categorical(
    df["Method"],
    categories=method_order,
    ordered=True
)

df = df.sort_values(
    ["Benchmark", "Method"]
).reset_index(drop=True)

df.to_csv(
    "heart_atlas_integration_benchmarking.csv",
    sep=";",
    decimal=".",
    index=False
)

# Read CSV (tab-separated, German decimal commas)
df = pd.read_csv(
    "heart_atlas_integration_benchmarking.csv",
    sep=";",
    decimal="."
)

# Ensure correct order of benchmarks
benchmark_order = ["ARI", "NMI", "cLISI", "iLISI"]
df["Benchmark"] = pd.Categorical(
    df["Benchmark"],
    categories=benchmark_order,
    ordered=True
)

# Pivot table for plotting
plot_df = df.pivot(
    index="Benchmark",
    columns="Method",
    values="Value"
).loc[benchmark_order]

# Colorblind-friendly colors
colors = {
    "BBKNN": "#0072B2",   # blue
    "scVI": "#E69F00",    # orange
    "Harmony": "#009E73"  # green
}

ax = plot_df.plot(
    kind="bar",
    figsize=(8, 5),
    width=0.75,
    color=[colors[m] for m in plot_df.columns],
    linewidth=0.6
)

ax.set_ylabel("Score", fontsize=12)
ax.set_xlabel("")

plt.xticks(rotation=0, ha="right")
plt.yticks(fontsize=10)
plt.legend(title="Method", frameon=False)

# Light horizontal grid
ax.grid(axis="y", linestyle="--", alpha=0.4)
ax.set_axisbelow(True)

plt.tight_layout()
plt.show()