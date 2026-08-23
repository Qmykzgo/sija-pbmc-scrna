"""
Clustering Module
=================
KNN graph construction, Leiden community detection, and UMAP
visualization for cell population identification.
"""

import logging
import warnings
from pathlib import Path

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)


def build_knn_graph(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Construct a K-nearest neighbor graph on the embedding space.

    The KNN graph is the foundation for both community detection
    (clustering) and UMAP embedding.

    Parameters
    ----------
    adata : ad.AnnData
        Integrated AnnData with corrected embeddings.
    config : dict
        Pipeline config with n_neighbors and n_pcs.

    Returns
    -------
    ad.AnnData
        AnnData with KNN graph in .obsp.
    """
    n_neighbors = config["clustering"]["n_neighbors"]
    n_pcs = config["clustering"]["n_pcs"]

    # Determine the best embedding to use
    if "X_scVI" in adata.obsm:
        use_rep = "X_scVI"
        n_pcs_use = None  # scVI latent is typically 10D, use all
    elif "X_pca_harmony" in adata.obsm:
        use_rep = "X_pca_harmony"
        n_pcs_use = n_pcs
    else:
        use_rep = "X_pca"
        n_pcs_use = n_pcs

    logger.info(f"Building KNN graph (n_neighbors={n_neighbors}, rep={use_rep})...")

    if n_pcs_use:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs_use, use_rep=use_rep)
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)

    logger.info(f"  KNN graph: {adata.obsp['distances'].shape}")
    return adata


def cluster_leiden(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Perform Leiden community detection on the KNN graph.

    Leiden is a refined version of Louvain that guarantees connected
    communities and generally produces better partitions.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData with KNN graph.
    config : dict
        Pipeline config with resolution.

    Returns
    -------
    ad.AnnData
        AnnData with 'leiden' cluster assignments in .obs.
    """
    resolution = config["clustering"]["resolution"]
    logger.info(f"Running Leiden clustering (resolution={resolution})...")

    try:
        sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2, directed=False)
    except Exception as e:
        logger.warning(f"Leiden with flavor='igraph' failed ({e}), falling back to default leiden...")
        sc.tl.leiden(adata, resolution=resolution)

    n_clusters = adata.obs["leiden"].nunique()
    logger.info(f"  Found {n_clusters} clusters")
    logger.info(f"  Cells per cluster:\n{adata.obs['leiden'].value_counts().sort_index().to_string()}")

    return adata


def compute_umap(adata: ad.AnnData) -> ad.AnnData:
    """
    Compute UMAP embedding for visualization.

    UMAP provides a 2D representation that preserves local and some
    global structure of the high-dimensional data.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData with KNN graph.

    Returns
    -------
    ad.AnnData
        AnnData with UMAP coordinates in .obsm['X_umap'].
    """
    logger.info("Computing UMAP embedding...")
    sc.tl.umap(adata)
    logger.info(f"  UMAP shape: {adata.obsm['X_umap'].shape}")
    return adata


def plot_clusters(adata: ad.AnnData, config: dict) -> None:
    """
    Generate clustering visualization plots.

    Creates UMAPs colored by cluster, batch, condition, and QC metrics.

    Parameters
    ----------
    adata : ad.AnnData
        Clustered AnnData.
    config : dict
        Pipeline config with figure paths.
    """
    fig_dir = Path(config["paths"]["figures_dir"])
    sc.set_figure_params(dpi=150, frameon=False, fontsize=12)

    # Leiden clusters with labels on data
    sc.pl.embedding(
        adata,
        color="leiden",
        basis="X_umap",
        legend_loc="on data",
        frameon=False,
        s=15,
        legend_fontsize=8,
        legend_fontoutline=1,
        show=False,
    )
    plt.savefig(fig_dir / "clustering_leiden.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Batch overlay
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sc.pl.umap(adata, color="id", ax=axes[0], show=False, title="Sample")
    sc.pl.umap(adata, color="condition", ax=axes[1], show=False, title="Condition")
    plt.tight_layout()
    plt.savefig(fig_dir / "clustering_batch_condition.png", dpi=200, bbox_inches="tight")
    plt.close()

    # QC metrics on UMAP
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    sc.pl.umap(adata, color="total_counts", ax=axes[0], show=False)
    sc.pl.umap(adata, color="n_genes_by_counts", ax=axes[1], show=False)
    sc.pl.umap(adata, color="pct_counts_mt", ax=axes[2], show=False)
    plt.tight_layout()
    plt.savefig(fig_dir / "clustering_qc_metrics.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Cell type if available
    if "Cell Type" in adata.obs.columns:
        sc.pl.embedding(
            adata,
            color="Cell Type",
            basis="X_umap",
            legend_loc="on data",
            frameon=False,
            s=15,
            legend_fontsize=8,
            legend_fontoutline=1,
            show=False,
        )
        plt.savefig(fig_dir / "clustering_celltype_author.png", dpi=200, bbox_inches="tight")
        plt.close()

    logger.info(f"Clustering plots saved to {fig_dir}")


def run(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Execute full clustering pipeline.

    Steps:
        1. Build KNN graph (on corrected embedding)
        2. Compute UMAP embedding
        3. Leiden community detection
        4. Generate cluster visualizations
        5. Save checkpoint

    Parameters
    ----------
    adata : ad.AnnData
        Integrated AnnData.
    config : dict
        Pipeline configuration.

    Returns
    -------
    ad.AnnData
        Clustered AnnData.
    """
    checkpoint = Path(config["paths"]["checkpoints_dir"]) / "05_clustered.h5ad"

    if checkpoint.exists():
        logger.info(f"Loading cached clustered data from {checkpoint}")
        return sc.read_h5ad(str(checkpoint))

    adata = build_knn_graph(adata, config)
    adata = compute_umap(adata)
    adata = cluster_leiden(adata, config)
    plot_clusters(adata, config)

    # Save
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(checkpoint))
    logger.info(f"Saved clustering checkpoint: {checkpoint}")

    return adata
