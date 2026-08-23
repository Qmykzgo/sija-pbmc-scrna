"""
Cell Type Annotation Module
============================
Automated cell type annotation using CellTypist with majority voting,
supplemented by marker gene validation and rank_genes_groups analysis.
"""

import logging
import warnings
from pathlib import Path

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)


def annotate_celltypist(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Annotate cell types using CellTypist reference model.

    CellTypist uses a logistic regression classifier trained on large
    reference atlases. Majority voting refines predictions by enforcing
    consistency within over-clustered groups.

    Parameters
    ----------
    adata : ad.AnnData
        Clustered AnnData with UMAP.
    config : dict
        Pipeline config with celltypist model name.

    Returns
    -------
    ad.AnnData
        AnnData with 'predicted_labels', 'majority_voting', 'conf_score' in .obs.
    """
    import celltypist
    from celltypist import models

    model_name = config["annotation"]["celltypist_model"]
    logger.info(f"Running CellTypist annotation (model: {model_name})...")

    # Download model if needed
    models.download_models(model=model_name)
    model = models.Model.load(model=model_name)

    # Run prediction with majority voting
    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=True,
    )

    # Transfer annotations safely into adata.obs preserving all layers and obsm
    if hasattr(predictions, "predicted_labels"):
        pred_df = predictions.predicted_labels
        for col in ["predicted_labels", "majority_voting", "conf_score"]:
            if col in pred_df.columns:
                adata.obs[col] = pred_df[col].values
    else:
        adata = predictions.to_adata(insert_labels=True)

    cell_type_col = "majority_voting" if "majority_voting" in adata.obs.columns else "predicted_labels"
    n_types = adata.obs[cell_type_col].nunique()
    logger.info(f"  Annotated {n_types} cell types")
    logger.info(
        f"  Cell type distribution:\n"
        f"{adata.obs[cell_type_col].value_counts().head(15).to_string()}"
    )

    return adata


def find_marker_genes(adata: ad.AnnData) -> ad.AnnData:
    """
    Find differentially expressed marker genes for each Leiden cluster.

    Uses Wilcoxon rank-sum test to identify genes that are significantly
    upregulated in each cluster compared to all other cells.

    Parameters
    ----------
    adata : ad.AnnData
        Clustered AnnData.

    Returns
    -------
    ad.AnnData
        AnnData with rank_genes_groups results in .uns.
    """
    logger.info("Finding marker genes per Leiden cluster...")

    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    logger.info("  Marker gene analysis complete.")
    return adata


def plot_annotations(adata: ad.AnnData, config: dict) -> None:
    """
    Generate annotation visualization plots.

    Creates:
    - UMAP colored by CellTypist majority voting labels
    - Marker gene dotplot and feature plots
    - Rank genes heatmap per cluster

    Parameters
    ----------
    adata : ad.AnnData
        Annotated AnnData.
    config : dict
        Pipeline config with marker genes and figure paths.
    """
    fig_dir = Path(config["paths"]["figures_dir"])
    sc.set_figure_params(dpi=150, frameon=False, fontsize=12)
    plt.rcParams["axes.grid"] = False

    # CellTypist annotations on UMAP
    if "majority_voting" in adata.obs.columns:
        sc.pl.embedding(
            adata,
            color="majority_voting",
            basis="X_umap",
            legend_loc="on data",
            frameon=False,
            s=15,
            legend_fontsize=7,
            legend_fontoutline=1,
            show=False,
        )
        plt.savefig(fig_dir / "annotation_celltypist.png", dpi=200, bbox_inches="tight")
        plt.close()

    # Marker gene feature plots
    markers_flat = []
    markers_config = config["annotation"]["markers"]
    for cell_type, genes in markers_config.items():
        for gene in genes:
            if gene in adata.var_names or (adata.raw is not None and gene in adata.raw.var_names):
                markers_flat.append(gene)

    # Deduplicate while preserving order
    markers_flat = list(dict.fromkeys(markers_flat))

    if markers_flat:
        # Select a representative subset for feature plot
        plot_genes = markers_flat[:12]
        sc.pl.umap(adata, color=plot_genes, ncols=4, show=False)
        plt.savefig(fig_dir / "annotation_marker_features.png", dpi=200, bbox_inches="tight")
        plt.close()

        # Dotplot grouped by cell type annotation
        group_key = "majority_voting" if "majority_voting" in adata.obs.columns else "leiden"
        sc.pl.dotplot(
            adata,
            var_names=markers_flat[:20],
            groupby=group_key,
            standard_scale="var",
            show=False,
        )
        plt.savefig(fig_dir / "annotation_dotplot.png", dpi=200, bbox_inches="tight")
        plt.close()

    # Rank genes per cluster
    if "rank_genes_groups" in adata.uns:
        sc.pl.rank_genes_groups(adata, n_genes=6, sharey=False, show=False)
        plt.savefig(fig_dir / "annotation_rank_genes.png", dpi=200, bbox_inches="tight")
        plt.close()

    logger.info(f"Annotation plots saved to {fig_dir}")


def run(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Execute full annotation pipeline.

    Steps:
        1. CellTypist automated annotation with majority voting
        2. Marker gene identification (Wilcoxon rank-sum)
        3. Generate annotation visualizations
        4. Save checkpoint

    Parameters
    ----------
    adata : ad.AnnData
        Clustered AnnData.
    config : dict
        Pipeline configuration.

    Returns
    -------
    ad.AnnData
        Annotated AnnData.
    """
    checkpoint = Path(config["paths"]["checkpoints_dir"]) / "06_annotated.h5ad"

    if checkpoint.exists():
        logger.info(f"Loading cached annotated data from {checkpoint}")
        return sc.read_h5ad(str(checkpoint))

    adata = annotate_celltypist(adata, config)
    adata = find_marker_genes(adata)
    plot_annotations(adata, config)

    # Save
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(checkpoint))
    logger.info(f"Saved annotation checkpoint: {checkpoint}")

    return adata
