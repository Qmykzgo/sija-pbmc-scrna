"""
Batch Integration Module
========================
Compare three batch correction methods (BBKNN, Harmony, scVI) and
benchmark with scib metrics to select the best integration.
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


def integrate_bbknn(adata: ad.AnnData, batch_key: str) -> ad.AnnData:
    """
    Batch correction using BBKNN (graph-based approach).

    BBKNN modifies the KNN graph so that each cell's neighbors are
    balanced across batches, without altering expression values.

    Parameters
    ----------
    adata : ad.AnnData
        Preprocessed AnnData with PCA.
    batch_key : str
        Column in .obs identifying batches.

    Returns
    -------
    ad.AnnData
        AnnData with corrected KNN graph and UMAP.
    """
    import bbknn

    logger.info("Running BBKNN integration...")
    adata_corr = adata.copy()
    if "X_pca" not in adata_corr.obsm:
        sc.pp.pca(adata_corr)
    bbknn.bbknn(adata_corr, batch_key=batch_key)
    sc.tl.umap(adata_corr)
    logger.info("  BBKNN integration complete.")
    return adata_corr


def integrate_harmony(adata: ad.AnnData, batch_key: str) -> ad.AnnData:
    """
    Batch correction using Harmony (iterative PCA correction).

    Harmony iteratively adjusts PCA embeddings to mix batches while
    preserving biological variation. It modifies the embedding space
    rather than the expression matrix.

    Parameters
    ----------
    adata : ad.AnnData
        Preprocessed AnnData with PCA.
    batch_key : str
        Column in .obs identifying batches.

    Returns
    -------
    ad.AnnData
        AnnData with corrected PCA and UMAP embeddings.
    """
    logger.info("Running Harmony integration...")
    adata_corr = adata.copy()

    sc.external.pp.harmony_integrate(adata_corr, key=batch_key)

    # Use corrected embeddings for downstream
    adata_corr.obsm["X_pca"] = adata_corr.obsm["X_pca_harmony"]
    sc.pp.neighbors(adata_corr, n_pcs=30)
    sc.tl.umap(adata_corr)

    logger.info("  Harmony integration complete.")
    return adata_corr


def integrate_scvi(adata: ad.AnnData, batch_key: str, config: dict) -> ad.AnnData:
    """
    Batch correction using scVI (variational autoencoder).

    scVI learns a probabilistic latent space that accounts for batch
    effects through a deep generative model. Requires raw unnormalized counts.

    Parameters
    ----------
    adata : ad.AnnData
        Preprocessed AnnData.
    batch_key : str
        Column in .obs identifying batches.
    config : dict
        Pipeline config with scVI hyperparameters.

    Returns
    -------
    ad.AnnData
        AnnData with scVI latent representation and UMAP.
    """
    import scvi

    logger.info("Running scVI integration...")
    adata_corr = adata.copy()

    scvi_cfg = config["integration"]["scvi"]

    # scVI requires raw integer counts (from layers['counts'] if available)
    if "counts" in adata_corr.layers:
        scvi.model.SCVI.setup_anndata(adata_corr, layer="counts", batch_key=batch_key)
    else:
        logger.warning("No 'counts' layer found for scVI. Using .X")
        scvi.model.SCVI.setup_anndata(adata_corr, batch_key=batch_key)

    model = scvi.model.SCVI(
        adata_corr,
        n_latent=scvi_cfg["n_latent"],
        n_hidden=scvi_cfg["n_hidden"],
        n_layers=scvi_cfg["n_layers"],
        dropout_rate=scvi_cfg["dropout_rate"],
    )

    model.train(
        max_epochs=scvi_cfg["max_epochs"],
        early_stopping=scvi_cfg["early_stopping"],
    )

    # Get latent representation
    adata_corr.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata_corr, use_rep="X_scVI")
    sc.tl.umap(adata_corr)

    logger.info("  scVI integration complete.")
    return adata_corr


def benchmark_integrations(
    adata_uncorr: ad.AnnData,
    results: dict,
    batch_key: str,
    label_key: str = "Cell Type",
) -> pd.DataFrame:
    """
    Benchmark integration methods using scib metrics.

    Evaluates both batch mixing (how well batches are integrated) and
    biological conservation (how well cell types are preserved).

    Parameters
    ----------
    adata_uncorr : ad.AnnData
        Uncorrected AnnData for comparison.
    results : dict
        Mapping method_name -> corrected AnnData.
    batch_key : str
        Batch column name.
    label_key : str
        Cell type column name.

    Returns
    -------
    pd.DataFrame
        Benchmark scores for each method.
    """
    try:
        import scib
    except ImportError:
        logger.warning("scib is not installed. Skipping benchmark metrics.")
        return pd.DataFrame()

    logger.info("Benchmarking integration methods with scib...")

    scores = {}
    for method_name, adata_corr in results.items():
        logger.info(f"  Evaluating {method_name}...")

        # Determine embedding key
        if method_name == "scVI":
            embed_key = "X_scVI"
        elif method_name == "Harmony":
            embed_key = "X_pca_harmony" if "X_pca_harmony" in adata_corr.obsm else "X_pca"
        else:
            embed_key = "X_pca"

        try:
            score = scib.metrics.metrics(
                adata_uncorr,
                adata_corr,
                batch_key=batch_key,
                label_key=label_key,
                embed=embed_key,
                isolated_labels_asw_=True,
                silhouette_=True,
                graph_conn_=True,
                nmi_=True,
                ari_=True,
                isolated_labels_f1_=False,
                trajectory_=False,
                pcr_=False,
                kBET_=False,
                lisi_graph_=False,
                clisi_=False,
                ilisi_=False,
                cell_cycle_=False,
                hvg_score_=False,
            )
            scores[method_name] = score
        except Exception as e:
            logger.warning(f"  scib evaluation failed for {method_name}: {e}")
            continue

    if scores:
        benchmark_df = pd.DataFrame(scores)
        logger.info(f"\nBenchmark Results:\n{benchmark_df.to_string()}")
        return benchmark_df

    logger.warning("No benchmark scores computed.")
    return pd.DataFrame()


def plot_integrations(results: dict, config: dict) -> None:
    """
    Plot UMAP embeddings for each integration method side by side.

    Parameters
    ----------
    results : dict
        Mapping method_name -> corrected AnnData.
    config : dict
        Pipeline config with figure paths.
    """
    fig_dir = Path(config["paths"]["figures_dir"])
    methods = list(results.keys())
    n = len(methods)
    if n == 0:
        logger.warning("No integration results to plot.")
        return

    # Plot batch mixing
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    if n == 1:
        axes = [axes]
    for ax, (name, adata_corr) in zip(axes, results.items()):
        sc.pl.umap(adata_corr, color="id", ax=ax, show=False, title=f"{name} (batch)")
    plt.tight_layout()
    plt.savefig(fig_dir / "integration_batch_comparison.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Plot cell type preservation if available
    first_adata = list(results.values())[0]
    cell_type_col = "Cell Type" if "Cell Type" in first_adata.obs.columns else ("cell_type" if "cell_type" in first_adata.obs.columns else None)
    if cell_type_col:
        fig, axes = plt.subplots(1, n, figsize=(6 * n, 5))
        if n == 1:
            axes = [axes]
        for ax, (name, adata_corr) in zip(axes, results.items()):
            sc.pl.umap(
                adata_corr, color=cell_type_col, ax=ax, show=False,
                title=f"{name} (cell type)", legend_loc="right margin",
            )
        plt.tight_layout()
        plt.savefig(fig_dir / "integration_celltype_comparison.png", dpi=200, bbox_inches="tight")
        plt.close()

    logger.info(f"Integration comparison plots saved to {fig_dir}")


def run(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Execute full integration pipeline.

    Steps:
        1. Run all configured integration methods
        2. Benchmark with scib metrics
        3. Generate comparison plots
        4. Select best method and save

    Parameters
    ----------
    adata : ad.AnnData
        Preprocessed AnnData with PCA.
    config : dict
        Pipeline configuration.

    Returns
    -------
    ad.AnnData
        Best-integrated AnnData (scVI by default).
    """
    checkpoint = Path(config["paths"]["checkpoints_dir"]) / "04_integrated.h5ad"

    if checkpoint.exists():
        logger.info(f"Loading cached integrated data from {checkpoint}")
        return sc.read_h5ad(str(checkpoint))

    batch_key = config["integration"]["batch_key"]
    methods = config["integration"]["methods"]
    results = {}

    for method in methods:
        try:
            if method == "bbknn":
                results["BBKNN"] = integrate_bbknn(adata, batch_key)
            elif method == "harmony":
                results["Harmony"] = integrate_harmony(adata, batch_key)
            elif method == "scvi":
                results["scVI"] = integrate_scvi(adata, batch_key, config)
        except Exception as e:
            logger.error(f"Integration method '{method}' failed: {e}")

    if not results:
        logger.warning("All integration methods failed. Using uncorrected data.")
        sc.pp.neighbors(adata, n_pcs=30)
        sc.tl.umap(adata)
        return adata

    # Plot comparisons
    plot_integrations(results, config)

    # Benchmark if cell type labels exist
    label_key = "Cell Type" if "Cell Type" in adata.obs.columns else None
    if label_key:
        benchmark_df = benchmark_integrations(adata, results, batch_key, label_key)
        if not benchmark_df.empty:
            fig_dir = Path(config["paths"]["figures_dir"])
            benchmark_df.to_csv(fig_dir / "integration_benchmark.csv")

    # Select scVI if available, else Harmony, else first available
    for preferred in ["scVI", "Harmony", "BBKNN"]:
        if preferred in results:
            best_adata = results[preferred]
            logger.info(f"Selected integration method: {preferred}")
            break

    # Save
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    best_adata.write_h5ad(str(checkpoint))
    logger.info(f"Saved integration checkpoint: {checkpoint}")

    return best_adata
