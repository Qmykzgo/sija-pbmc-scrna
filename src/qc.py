"""
Quality Control Module
======================
Comprehensive QC filtering for scRNA-seq data including:
- Mitochondrial and ribosomal gene annotation
- Cell and gene count filtering
- Doublet detection with Scrublet
- QC visualization (violin plots, scatter plots)
"""

import logging
import warnings
from pathlib import Path

import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)


def annotate_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    """
    Calculate QC metrics: mitochondrial %, ribosomal %, total counts.

    Annotates genes with boolean flags for MT and ribosomal genes,
    then computes per-cell QC statistics.

    Parameters
    ----------
    adata : ad.AnnData
        Raw merged AnnData.

    Returns
    -------
    ad.AnnData
        AnnData with QC metrics in .obs and gene annotations in .var.
    """
    logger.info("Annotating QC metrics...")

    # Mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    n_mt = adata.var["mt"].sum()
    logger.info(f"  Found {n_mt} mitochondrial genes")

    # Ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    n_ribo = adata.var["ribo"].sum()
    logger.info(f"  Found {n_ribo} ribosomal genes")

    # Calculate metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # Ensure raw count layer exists
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # Total counts shorthand as 1D array
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()

    logger.info(
        f"  Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}\n"
        f"  Median counts/cell: {adata.obs['total_counts'].median():.0f}\n"
        f"  Median MT%: {adata.obs['pct_counts_mt'].median():.1f}%"
    )

    return adata


def detect_doublets(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Detect doublets using Scrublet.

    Scrublet simulates doublets by combining random pairs of cells and
    scores each real cell by its similarity to simulated doublets.
    Runs batch-aware doublet detection across samples.

    Parameters
    ----------
    adata : ad.AnnData
        QC-annotated AnnData.
    config : dict
        Pipeline config with scrublet_threshold.

    Returns
    -------
    ad.AnnData
        AnnData with 'predicted_doublet' and 'doublet_score' in .obs.
    """
    threshold = config["qc"]["scrublet_threshold"]
    logger.info(f"Running Scrublet doublet detection (threshold={threshold})...")

    batch_key = "id" if "id" in adata.obs.columns else ("batch" if "batch" in adata.obs.columns else None)
    try:
        if batch_key:
            sc.pp.scrublet(adata, batch_key=batch_key, threshold=threshold)
        else:
            sc.pp.scrublet(adata, threshold=threshold)
    except Exception as e:
        logger.warning(f"Batch-aware scrublet failed ({e}), attempting standard scrublet...")
        sc.pp.scrublet(adata, threshold=threshold)

    n_doublets = int(adata.obs["predicted_doublet"].fillna(False).sum())
    n_total = adata.n_obs
    pct = n_doublets / n_total * 100
    logger.info(f"  Detected {n_doublets:,} doublets ({pct:.1f}%) out of {n_total:,} cells")

    return adata


def filter_cells_and_genes(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Apply cell and gene filtering based on QC thresholds.

    Filtering steps:
        1. Remove predicted doublets
        2. Filter cells with too few/many genes
        3. Filter cells with too few/many UMI counts
        4. Filter cells with high mitochondrial %
        5. Filter genes expressed in too few cells

    Parameters
    ----------
    adata : ad.AnnData
        Doublet-annotated AnnData.
    config : dict
        Pipeline config with QC thresholds.

    Returns
    -------
    ad.AnnData
        Filtered AnnData.
    """
    qc = config["qc"]
    n_before = adata.n_obs
    logger.info(f"Filtering cells (starting with {n_before:,})...")

    # Remove doublets safely
    if "predicted_doublet" in adata.obs.columns:
        adata = adata[~adata.obs["predicted_doublet"].fillna(False).astype(bool)].copy()
        logger.info(f"  After doublet removal: {adata.n_obs:,}")

    # Gene and cell count filters
    sc.pp.filter_cells(adata, min_genes=qc["min_genes"])
    sc.pp.filter_genes(adata, min_cells=qc["min_cells"])

    # Count thresholds
    adata = adata[adata.obs["total_counts"] > qc["min_counts"], :].copy()
    adata = adata[adata.obs["total_counts"] < qc["max_counts"], :].copy()

    # Gene count thresholds
    adata = adata[adata.obs["n_genes_by_counts"] > qc["min_genes_by_counts"], :].copy()
    adata = adata[adata.obs["n_genes_by_counts"] < qc["max_genes_by_counts"], :].copy()

    # Mitochondrial threshold
    adata = adata[adata.obs["pct_counts_mt"] < qc["max_mt_pct"], :].copy()

    n_after = adata.n_obs
    n_removed = n_before - n_after
    logger.info(
        f"  After all filters: {n_after:,} cells "
        f"(removed {n_removed:,}, {n_removed/n_before*100:.1f}%)"
    )
    logger.info(f"  Remaining genes: {adata.n_vars:,}")

    return adata


def remove_uninformative_genes(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Remove mitochondrial, ribosomal, and MALAT1 genes.

    These genes are informative for QC but can dominate variation
    and obscure biological signals in downstream analysis.

    Parameters
    ----------
    adata : ad.AnnData
        Filtered AnnData.
    config : dict
        Pipeline config with gene prefixes to remove.

    Returns
    -------
    ad.AnnData
        AnnData with uninformative genes removed.
    """
    prefixes = config["qc"]["remove_gene_prefixes"]
    logger.info(f"Removing genes with prefixes: {prefixes}")

    keep = np.ones(adata.n_vars, dtype=bool)
    for prefix in prefixes:
        mask = adata.var_names.str.startswith(prefix)
        n_remove = mask.sum()
        keep = keep & ~mask
        logger.info(f"  {prefix}: removing {n_remove} genes")

    adata = adata[:, keep].copy()
    logger.info(f"  Remaining: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    return adata


def plot_qc(adata: ad.AnnData, config: dict, stage: str = "pre") -> None:
    """
    Generate QC visualization plots.

    Creates violin plots, scatter plots, and highest-expressed gene
    bar charts for quality assessment.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData at the current QC stage.
    config : dict
        Pipeline config with figure paths.
    stage : str
        Either "pre" or "post" filtering.
    """
    fig_dir = Path(config["paths"]["figures_dir"])
    fig_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Generating {stage}-filtering QC plots...")

    # Violin plots
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    sc.pl.violin(adata, keys="n_genes_by_counts", ax=axes[0], show=False)
    sc.pl.violin(adata, keys="total_counts", ax=axes[1], show=False)
    sc.pl.violin(adata, keys="pct_counts_mt", ax=axes[2], show=False)
    sc.pl.violin(adata, keys="pct_counts_ribo", ax=axes[3], show=False)
    plt.suptitle(f"QC Metrics ({stage}-filter)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(fig_dir / f"qc_violin_{stage}.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Scatter: counts vs genes
    fig, ax = plt.subplots(figsize=(5, 5))
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", size=1, ax=ax, show=False)
    plt.savefig(fig_dir / f"qc_scatter_{stage}.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Highest expressed genes
    sc.pl.highest_expr_genes(adata, n_top=20, show=False)
    plt.savefig(fig_dir / f"qc_highest_genes_{stage}.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Cells per sample bar chart
    fig, ax = plt.subplots(figsize=(6, 2))
    adata.obs["id"].value_counts().plot.barh(ax=ax)
    ax.set_xlabel("Number of cells")
    ax.set_title(f"Cells per sample ({stage}-filter)")
    plt.tight_layout()
    plt.savefig(fig_dir / f"qc_cells_per_sample_{stage}.png", dpi=200, bbox_inches="tight")
    plt.close()

    logger.info(f"  Saved QC plots to {fig_dir}")


def run(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Execute the full QC pipeline.

    Steps:
        1. Annotate QC metrics (MT%, ribo%, counts)
        2. Generate pre-filter QC plots
        3. Detect doublets (Scrublet)
        4. Filter cells and genes
        5. Remove uninformative genes
        6. Generate post-filter QC plots
        7. Save checkpoint

    Parameters
    ----------
    adata : ad.AnnData
        Raw merged AnnData.
    config : dict
        Pipeline configuration.

    Returns
    -------
    ad.AnnData
        QC-filtered AnnData.
    """
    checkpoint = Path(config["paths"]["checkpoints_dir"]) / "02_qc_filtered.h5ad"

    if checkpoint.exists():
        logger.info(f"Loading cached QC data from {checkpoint}")
        return sc.read_h5ad(str(checkpoint))

    # 1. QC metrics
    adata = annotate_qc_metrics(adata)

    # 2. Pre-filter plots
    plot_qc(adata, config, stage="pre")

    # 3. Doublet detection
    adata = detect_doublets(adata, config)

    # 4. Filter
    adata = filter_cells_and_genes(adata, config)

    # 5. Remove MT/ribo/MALAT1
    adata = remove_uninformative_genes(adata, config)

    # 6. Post-filter plots
    plot_qc(adata, config, stage="post")

    # 7. Checkpoint
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(checkpoint))
    logger.info(f"Saved QC checkpoint: {checkpoint}")

    return adata
