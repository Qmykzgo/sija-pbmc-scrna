"""
Visualization Module
====================
Publication-quality figure generation. Consolidates key results
into a multi-panel summary figure suitable for GitHub README.
"""

import logging
import warnings
from pathlib import Path

import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)


def set_publication_style():
    """Configure matplotlib for publication-quality figures."""
    sc.set_figure_params(
        dpi=150,
        dpi_save=300,
        frameon=False,
        vector_friendly=False,
        fontsize=12,
        format="pdf",
        transparent=False,
    )
    plt.rcParams.update({
        "axes.grid": False,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
    })


def create_summary_figure(adata: ad.AnnData, config: dict) -> None:
    """
    Create a multi-panel summary figure for the analysis.

    Layout:
        Panel A: UMAP colored by cell type (CellTypist)
        Panel B: UMAP colored by condition
        Panel C: Marker gene dotplot
        Panel D: Cells per sample barplot

    Parameters
    ----------
    adata : ad.AnnData
        Final annotated AnnData.
    config : dict
        Pipeline config with figure paths.
    """
    fig_dir = Path(config["paths"]["figures_dir"])
    set_publication_style()

    logger.info("Creating summary figure...")

    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 2, hspace=0.3, wspace=0.3)

    # Panel A: Cell types
    ax_a = fig.add_subplot(gs[0, 0])
    cell_type_col = "majority_voting" if "majority_voting" in adata.obs.columns else "leiden"
    sc.pl.umap(adata, color=cell_type_col, ax=ax_a, show=False, title="Cell Types", s=10)
    ax_a.text(-0.1, 1.05, "A", transform=ax_a.transAxes, fontsize=18, fontweight="bold")

    # Panel B: Condition
    ax_b = fig.add_subplot(gs[0, 1])
    cond_col = "clin_group" if "clin_group" in adata.obs.columns else "condition"
    if cond_col in adata.obs.columns:
        sc.pl.umap(adata, color=cond_col, ax=ax_b, show=False, title="Disease Status", s=10)
    ax_b.text(-0.1, 1.05, "B", transform=ax_b.transAxes, fontsize=18, fontweight="bold")

    # Panel C: Key marker genes
    ax_c = fig.add_subplot(gs[1, 0])
    key_markers = ["CD3E", "CD14", "CD19", "NKG7"]
    available = [g for g in key_markers if g in adata.var_names or
                 (adata.raw is not None and g in adata.raw.var_names)]
    if available:
        sc.pl.umap(adata, color=available[0], ax=ax_c, show=False, title=f"{available[0]} expression", s=10)
    ax_c.text(-0.1, 1.05, "C", transform=ax_c.transAxes, fontsize=18, fontweight="bold")

    # Panel D: Cells per sample
    ax_d = fig.add_subplot(gs[1, 1])
    cell_counts = adata.obs["id"].value_counts()
    colors = plt.cm.tab10(np.linspace(0, 1, len(cell_counts)))
    cell_counts.plot.barh(ax=ax_d, color=colors)
    ax_d.set_xlabel("Number of cells")
    ax_d.set_title("Cells per Sample")
    ax_d.text(-0.1, 1.05, "D", transform=ax_d.transAxes, fontsize=18, fontweight="bold")

    plt.savefig(fig_dir / "summary_figure.png", dpi=300, bbox_inches="tight")
    plt.savefig(fig_dir / "summary_figure.pdf", bbox_inches="tight")
    plt.close()

    logger.info(f"  Summary figure saved to {fig_dir}")


def create_marker_heatmap(adata: ad.AnnData, config: dict) -> None:
    """
    Create a comprehensive marker gene heatmap grouped by cell type.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated AnnData.
    config : dict
        Pipeline config with marker genes.
    """
    fig_dir = Path(config["paths"]["figures_dir"])
    set_publication_style()

    markers_config = config["annotation"]["markers"]
    cell_type_col = "majority_voting" if "majority_voting" in adata.obs.columns else "leiden"

    # Build ordered marker dict (only genes present in data)
    var_names = set(adata.var_names)
    if adata.raw is not None:
        var_names = var_names | set(adata.raw.var_names)

    ordered_markers = {}
    for cell_type, genes in markers_config.items():
        available = [g for g in genes if g in var_names]
        if available:
            ordered_markers[cell_type] = available

    if ordered_markers:
        sc.pl.dotplot(
            adata,
            var_names=ordered_markers,
            groupby=cell_type_col,
            standard_scale="var",
            show=False,
        )
        plt.savefig(fig_dir / "marker_heatmap.png", dpi=300, bbox_inches="tight")
        plt.close()
        logger.info("  Marker heatmap saved.")


def run(adata: ad.AnnData, config: dict) -> None:
    """
    Generate all publication-quality visualizations.

    Parameters
    ----------
    adata : ad.AnnData
        Final annotated AnnData.
    config : dict
        Pipeline configuration.
    """
    fig_dir = Path(config["paths"]["figures_dir"])
    fig_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Generating publication-quality figures...")
    logger.info("=" * 60)

    create_summary_figure(adata, config)
    create_marker_heatmap(adata, config)

    logger.info(f"All figures saved to {fig_dir}")
