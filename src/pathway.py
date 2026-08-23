"""
Pathway Activity Module
=======================
Infer pathway activity scores per cell using Decoupler + PROGENy.
Visualize pathway activity across cell types and conditions.
"""

import logging
import warnings
from pathlib import Path

import scanpy as sc
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)


def score_pathway_activity(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Infer pathway activity per cell using Decoupler with PROGENy.

    PROGENy infers 14 pathway activities (JAK-STAT, TNFa, NFkB, etc.)
    from gene expression footprints. Decoupler provides a unified
    interface for running multivariate linear models.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated AnnData.
    config : dict
        Pipeline config with PROGENy parameters.

    Returns
    -------
    ad.AnnData
        AnnData with pathway scores in .obsm['progeny'].
    """
    import decoupler as dc

    pathway_cfg = config["pathway"]["progeny"]
    logger.info("Inferring pathway activity with Decoupler + PROGENy...")

    # Get PROGENy model (gene weights for 14 pathways)
    if hasattr(dc, "op") and hasattr(dc.op, "progeny"):
        progeny = dc.op.progeny(
            organism=pathway_cfg["organism"],
            top=pathway_cfg["top_genes"],
        )
    elif hasattr(dc, "get_progeny"):
        progeny = dc.get_progeny(
            organism=pathway_cfg["organism"],
            top=pathway_cfg["top_genes"],
        )
    else:
        raise AttributeError("Could not find progeny resource in decoupler.")

    logger.info(f"  PROGENy model: {progeny.shape[0]} gene-pathway associations")

    # Run multivariate linear model
    if hasattr(dc, "mt") and hasattr(dc.mt, "mlm"):
        dc.mt.mlm(
            data=adata,
            net=progeny,
            verbose=True,
        )
    elif hasattr(dc, "run_mlm"):
        dc.run_mlm(
            mat=adata,
            net=progeny,
            source="source",
            target="target",
            weight="weight",
            verbose=True,
        )

    # Results are stored in adata.obsm['mlm_estimate'] and adata.obsm['mlm_pvals']
    if "mlm_estimate" in adata.obsm:
        logger.info(f"  Pathway scores shape: {adata.obsm['mlm_estimate'].shape}")

        # Also extract as a dataframe for easy access
        pathway_scores = adata.obsm["mlm_estimate"]
        if hasattr(pathway_scores, "columns"):
            pathways = list(pathway_scores.columns)
        else:
            pathways = [f"Pathway_{i}" for i in range(pathway_scores.shape[1])]
        logger.info(f"  Pathways scored: {', '.join(pathways[:10])}...")
    else:
        logger.warning("  Pathway scoring did not produce expected output.")

    return adata


def plot_pathway_activity(adata: ad.AnnData, config: dict) -> None:
    """
    Visualize pathway activity across cell types and conditions.

    Generates:
    - Pathway activity heatmap per cell type
    - UMAPs colored by key pathway activities
    - Violin plots of pathway activity by condition

    Parameters
    ----------
    adata : ad.AnnData
        AnnData with pathway scores.
    config : dict
        Pipeline config with figure paths.
    """
    import numpy as np

    fig_dir = Path(config["paths"]["figures_dir"])

    if "mlm_estimate" not in adata.obsm:
        logger.warning("No pathway scores found. Skipping pathway plots.")
        return

    pathway_scores = adata.obsm["mlm_estimate"]

    # Determine grouping columns
    cell_type_col = "majority_voting" if "majority_voting" in adata.obs.columns else "leiden"
    condition_col = "clin_group" if "clin_group" in adata.obs.columns else "condition"

    # Build acts AnnData directly from mlm_estimate
    if hasattr(pathway_scores, "values") and hasattr(pathway_scores, "columns"):
        acts = ad.AnnData(
            X=np.asarray(pathway_scores.values, dtype=float),
            obs=adata.obs.copy(),
            var=pd.DataFrame(index=list(pathway_scores.columns)),
        )
        pathways = list(pathway_scores.columns)
    elif hasattr(pathway_scores, "shape"):
        pathways = [f"Pathway_{i}" for i in range(pathway_scores.shape[1])]
        acts = ad.AnnData(
            X=np.asarray(pathway_scores, dtype=float),
            obs=adata.obs.copy(),
            var=pd.DataFrame(index=pathways),
        )
    else:
        logger.warning("Could not convert mlm_estimate to AnnData.")
        return

    # Ensure UMAP coordinates are available on acts AnnData
    if "X_umap" in adata.obsm and "X_umap" not in acts.obsm:
        acts.obsm["X_umap"] = adata.obsm["X_umap"]

    # Mean activity per cell type heatmap
    try:
        sc.pl.matrixplot(
            acts,
            var_names=pathways,
            groupby=cell_type_col,
            standard_scale="var",
            cmap="RdBu_r",
            show=False,
        )
        plt.savefig(fig_dir / "pathway_heatmap.png", dpi=200, bbox_inches="tight")
        plt.close()
    except Exception as e:
        logger.warning(f"  Pathway heatmap failed: {e}")

    # UMAP colored by key pathways
    key_pathways = ["JAK-STAT", "TNFa", "NFkB", "Trail"]
    available = [p for p in key_pathways if p in acts.var_names]
    if not available:
        # Fallback to whatever pathways are present
        available = list(acts.var_names[:4])

    if available and "X_umap" in acts.obsm:
        try:
            sc.pl.umap(acts, color=available, ncols=2, show=False, vmin=-3, vmax=3, cmap="RdBu_r")
            plt.savefig(fig_dir / "pathway_umap.png", dpi=200, bbox_inches="tight")
            plt.close()
        except Exception as e:
            logger.warning(f"  Pathway UMAP failed: {e}")

    # Violin: pathway activity by condition
    if condition_col in acts.obs.columns and available:
        for pathway in available[:4]:
            try:
                sc.pl.violin(acts, keys=pathway, groupby=condition_col, show=False)
                plt.savefig(
                    fig_dir / f"pathway_violin_{pathway.replace('-', '_')}.png",
                    dpi=200, bbox_inches="tight",
                )
                plt.close()
            except Exception as e:
                logger.warning(f"  Pathway violin for {pathway} failed: {e}")

    logger.info(f"Pathway plots saved to {fig_dir}")


def run(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Execute full pathway analysis pipeline.

    Steps:
        1. Score pathway activity (Decoupler + PROGENy)
        2. Generate pathway visualizations
        3. Save checkpoint

    Parameters
    ----------
    adata : ad.AnnData
        Annotated AnnData.
    config : dict
        Pipeline configuration.

    Returns
    -------
    ad.AnnData
        AnnData with pathway activity scores.
    """
    checkpoint = Path(config["paths"]["checkpoints_dir"]) / "08_pathway.h5ad"

    if checkpoint.exists():
        logger.info(f"Loading cached pathway data from {checkpoint}")
        return sc.read_h5ad(str(checkpoint))

    try:
        adata = score_pathway_activity(adata, config)
        plot_pathway_activity(adata, config)
    except Exception as e:
        logger.error(f"Pathway analysis failed: {e}")
        logger.error("This is non-critical. Continuing pipeline.")

    # Save
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(checkpoint))
    logger.info(f"Saved pathway checkpoint: {checkpoint}")

    return adata
