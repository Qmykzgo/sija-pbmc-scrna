# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: sjia-scrna
#     language: python
#     name: sjia-scrna
# ---

# %% [markdown]
# # 🧬 SJIA PBMC Single-Cell RNA-seq Analysis
#
# **Interactive walkthrough** of the full scRNA-seq analysis pipeline.
# This notebook mirrors the modular pipeline but allows step-by-step
# exploration with inline visualizations.
#
# **Dataset:** GSE207633 — PBMCs from pediatric SJIA patients
#
# **Pipeline:**
# ```
# Data Loading → QC → Preprocessing → Integration → Clustering → Annotation → Composition → Pathways
# ```

# %% [markdown]
# ## Setup

# %%
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import yaml

warnings.filterwarnings("ignore")
sc.set_figure_params(dpi=100, frameon=False, fontsize=12)
plt.rcParams["axes.grid"] = False

import sys
from pathlib import Path

# Load config (robust to running from project root or notebooks/ dir)
config_path = Path("config.yaml") if Path("config.yaml").exists() else Path("../config.yaml")
if str(Path(".").resolve()) not in sys.path:
    sys.path.insert(0, str(Path(".").resolve()))
if str(Path("..").resolve()) not in sys.path:
    sys.path.insert(0, str(Path("..").resolve()))

with open(config_path) as f:
    config = yaml.safe_load(f)

print(f"Scanpy version: {sc.__version__}")

# %% [markdown]
# ## 1. Data Loading
#
# Load 4 PBMC samples from GEO (2 Healthy controls, 2 MAS patients).

# %%
from src import data_loader

adata = data_loader.run(config)
print(adata)

# %%
# Cells per sample
adata.obs["id"].value_counts().plot.barh(figsize=(6, 2))
plt.xlabel("Number of cells")
plt.title("Cells per sample (pre-QC)")
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 2. Quality Control
#
# - Annotate mitochondrial and ribosomal genes
# - Detect doublets with Scrublet
# - Filter by count thresholds

# %%
from src import qc

adata = qc.run(adata, config)
print(f"After QC: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# %%
# QC violin plots
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
sc.pl.violin(adata, keys="n_genes_by_counts", ax=axes[0], show=False)
sc.pl.violin(adata, keys="total_counts", ax=axes[1], show=False)
sc.pl.violin(adata, keys="pct_counts_mt", ax=axes[2], show=False)
sc.pl.violin(adata, keys="pct_counts_ribo", ax=axes[3], show=False)
plt.tight_layout()
plt.show()

# %%
# Counts vs genes scatter
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", size=1)

# %%
# Highest expressed genes
sc.pl.highest_expr_genes(adata, n_top=20, palette="Blues", width=0.3)

# %% [markdown]
# ## 3. Preprocessing
#
# Normalize → log1p → HVG selection → PCA

# %%
from src import preprocessing

adata = preprocessing.run(adata, config)
print(f"After preprocessing: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# %%
# PCA variance ratio
sc.pl.pca_variance_ratio(adata, n_pcs=40, log=True)

# %%
# PCA colored by batch
sc.pl.pca(
    adata,
    color=["id", "id"],
    dimensions=[(0, 1), (2, 3)],
    ncols=2, size=2, wspace=0.3,
)

# %% [markdown]
# ## 4. Batch Integration
#
# Compare BBKNN, Harmony, and scVI. Select best method via scib metrics.

# %%
from src import integration

adata = integration.run(adata, config)
print(f"After integration: {adata.n_obs:,} cells")

# %%
# UMAP after integration
sc.pl.umap(adata, color=["id", "condition"], wspace=0.4)

# %% [markdown]
# ## 5. Clustering
#
# Leiden community detection on KNN graph.

# %%
from src import clustering

adata = clustering.run(adata, config)

# %%
# Clusters on UMAP
sc.pl.embedding(
    adata, color="leiden", basis="X_umap",
    legend_loc="on data", frameon=False, s=15,
    legend_fontsize=8, legend_fontoutline=1,
)

# %% [markdown]
# ## 6. Cell Type Annotation
#
# CellTypist automated annotation with majority voting.

# %%
from src import annotation

adata = annotation.run(adata, config)

# %%
# CellTypist annotations
if "majority_voting" in adata.obs.columns:
    sc.pl.embedding(
        adata, color="majority_voting", basis="X_umap",
        legend_loc="on data", frameon=False, s=15,
        legend_fontsize=7, legend_fontoutline=1,
    )

# %%
# Marker genes
markers = config["annotation"]["markers"]
all_markers = [g for genes in markers.values() for g in genes
               if g in adata.var_names or (adata.raw is not None and g in adata.raw.var_names)]
all_markers = list(dict.fromkeys(all_markers))
sc.pl.umap(adata, color=all_markers[:8], ncols=4)

# %%
# Dotplot
cell_type_col = "majority_voting" if "majority_voting" in adata.obs.columns else "leiden"
sc.pl.dotplot(adata, var_names=all_markers[:20], groupby=cell_type_col, standard_scale="var")

# %% [markdown]
# ## 7. Compositional Analysis
#
# scCODA and Milo for differential abundance testing.

# %%
from src import differential

adata = differential.run(adata, config)

# %% [markdown]
# ## 8. Pathway Activity
#
# Decoupler + PROGENy pathway scoring.

# %%
from src import pathway

adata = pathway.run(adata, config)

# %% [markdown]
# ## 9. Summary Figure

# %%
from src import visualization

visualization.run(adata, config)
print("✅ Pipeline complete! Check results/figures/ for all outputs.")

# %%
# Final object summary
print(adata)
