"""
Preprocessing Module
====================
Normalization, log-transformation, highly variable gene (HVG)
selection, and PCA dimensionality reduction.
"""

import logging
import warnings
from pathlib import Path

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)


def normalize_and_log(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Normalize counts to a target sum per cell and apply log1p transform.

    This brings cells to comparable total count levels and stabilizes
    variance for downstream analysis.

    Parameters
    ----------
    adata : ad.AnnData
        QC-filtered AnnData with raw counts.
    config : dict
        Pipeline config with target_sum.

    Returns
    -------
    ad.AnnData
        Normalized, log-transformed AnnData.
    """
    target_sum = config["preprocessing"]["target_sum"]
    logger.info(f"Normalizing to target_sum={target_sum:,.0f} and log1p transform...")

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    # Store raw counts for differential expression later
    adata.raw = adata

    logger.info("  Normalization complete. Raw counts stored in adata.raw")
    return adata


def select_hvg(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Identify and subset to highly variable genes (HVGs).

    HVGs capture the most biologically informative variation while
    reducing noise from lowly-variable housekeeping genes.
    Uses batch-aware selection across samples.

    Parameters
    ----------
    adata : ad.AnnData
        Normalized AnnData.
    config : dict
        Pipeline config with n_top_genes.

    Returns
    -------
    ad.AnnData
        AnnData subsetted to HVGs.
    """
    n_top = config["preprocessing"]["n_top_genes"]
    logger.info(f"Selecting top {n_top:,} highly variable genes (batch-aware)...")

    batch_key = "id" if "id" in adata.obs.columns else None
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, batch_key=batch_key)
    n_hvg = adata.var["highly_variable"].sum()
    logger.info(f"  {n_hvg} HVGs identified")

    # Plot HVG selection
    fig_dir = Path(config["paths"]["figures_dir"])
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(fig_dir / "preprocessing_hvg.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Subset to HVGs
    adata = adata[:, adata.var.highly_variable].copy()
    logger.info(f"  Subsetted to {adata.n_vars} HVGs")

    return adata


def run_pca(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Perform scaling and PCA dimensionality reduction.

    Scales genes to unit variance and clips extreme values, then
    reduces high-dimensional gene space to top principal components.

    Parameters
    ----------
    adata : ad.AnnData
        HVG-subsetted AnnData.
    config : dict
        Pipeline config with n_pcs.

    Returns
    -------
    ad.AnnData
        AnnData with PCA embedding in .obsm['X_pca'].
    """
    n_pcs = config["preprocessing"]["n_pcs"]
    logger.info("Scaling HVGs to unit variance (max_value=10)...")
    sc.pp.scale(adata, max_value=10)

    logger.info(f"Running PCA (n_components={n_pcs})...")
    sc.pp.pca(adata, n_comps=n_pcs)

    # Plot variance ratio
    fig_dir = Path(config["paths"]["figures_dir"])
    sc.pl.pca_variance_ratio(adata, n_pcs=min(40, n_pcs), log=True, show=False)
    plt.savefig(fig_dir / "preprocessing_pca_variance.png", dpi=200, bbox_inches="tight")
    plt.close()

    # PCA colored by batch, condition, and QC
    available_colors = [c for c in ["id", "condition", "pct_counts_ribo", "total_counts"] if c in adata.obs.columns]
    if available_colors:
        sc.pl.pca(
            adata,
            color=available_colors,
            ncols=2,
            size=2,
            wspace=0.3,
            show=False,
        )
        plt.savefig(fig_dir / "preprocessing_pca_batch.png", dpi=200, bbox_inches="tight")
        plt.close()

    logger.info(f"  PCA shape: {adata.obsm['X_pca'].shape}")
    return adata


def run(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Execute full preprocessing pipeline.

    Steps:
        1. Normalize counts (target_sum=10000)
        2. Log1p transform
        3. Select highly variable genes
        4. PCA dimensionality reduction
        5. Save checkpoint

    Parameters
    ----------
    adata : ad.AnnData
        QC-filtered AnnData.
    config : dict
        Pipeline configuration.

    Returns
    -------
    ad.AnnData
        Preprocessed AnnData with PCA embeddings.
    """
    checkpoint = Path(config["paths"]["checkpoints_dir"]) / "03_preprocessed.h5ad"

    if checkpoint.exists():
        logger.info(f"Loading cached preprocessed data from {checkpoint}")
        return sc.read_h5ad(str(checkpoint))

    adata = normalize_and_log(adata, config)
    adata = select_hvg(adata, config)
    adata = run_pca(adata, config)

    # Save
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(checkpoint))
    logger.info(f"Saved preprocessing checkpoint: {checkpoint}")

    return adata
