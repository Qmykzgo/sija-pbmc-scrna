"""
Data Loading Module
===================
Download scRNA-seq data from GEO (GSE207633) and prepare AnnData objects.
Handles .h5 matrix loading, sample merging, and clinical metadata attachment.
"""

import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import logging
import warnings
from pathlib import Path

import scanpy as sc
import pandas as pd
import anndata as ad

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)


def download_geo_data(config: dict) -> Path:
    """
    Download supplementary files from GEO accession GSE207633.

    The dataset contains .h5 count matrices for each sample and
    clinical metadata in .xlsx format.

    Parameters
    ----------
    config : dict
        Pipeline configuration with paths and dataset info.

    Returns
    -------
    Path
        Path to the data directory containing downloaded files.
    """
    import GEOparse

    data_dir = Path(config["paths"]["data_dir"])
    data_dir.mkdir(parents=True, exist_ok=True)

    geo_acc = config["dataset"]["geo_accession"]
    logger.info(f"Downloading GEO dataset: {geo_acc}")

    # Download GEO series
    gse = GEOparse.get_GEO(geo=geo_acc, destdir=str(data_dir), silent=True)

    # Download supplementary files (h5 matrices, metadata xlsx)
    for gsm_name, gsm in gse.gsms.items():
        for url in gsm.metadata.get("supplementary_file", []):
            # Prefer HTTPS over FTP for firewall resilience
            url = url.replace("ftp://", "https://")
            filename = url.split("/")[-1]
            filepath = data_dir / filename
            if not filepath.exists():
                logger.info(f"  Downloading {filename}...")
                import urllib.request
                urllib.request.urlretrieve(url, filepath)
            else:
                logger.info(f"  {filename} already exists, skipping.")

    # Also download series-level supplementary files (metadata xlsx)
    for url in gse.metadata.get("supplementary_file", []):
        url = url.replace("ftp://", "https://")
        filename = url.split("/")[-1]
        filepath = data_dir / filename
        if not filepath.exists():
            logger.info(f"  Downloading {filename}...")
            import urllib.request
            urllib.request.urlretrieve(url, filepath)

    logger.info(f"Data downloaded to {data_dir}")
    return data_dir


def load_h5_matrices(config: dict) -> dict:
    """
    Load .h5 count matrices for each sample.

    Parameters
    ----------
    config : dict
        Pipeline configuration.

    Returns
    -------
    dict
        Mapping of sample_id -> AnnData object.
    """
    data_dir = Path(config["paths"]["data_dir"])
    adatas = {}

    for sample_info in config["dataset"]["samples"]:
        sample_id = sample_info["id"]
        # Find matching .h5 file
        h5_files = list(data_dir.glob(f"*{sample_info['geo_id']}*.h5"))

        if not h5_files:
            logger.warning(f"No .h5 file found for {sample_id}, skipping.")
            continue

        h5_path = h5_files[0]
        logger.info(f"Loading {h5_path.name} -> {sample_id}")

        adata = sc.read_10x_h5(str(h5_path))
        adata.var_names_make_unique()

        # Preserve raw integer counts in layers
        adata.layers["counts"] = adata.X.copy()

        # Attach sample-level metadata
        adata.obs["barcode"] = adata.obs_names.astype(str)
        adata.obs["raw_id"] = sample_id
        adata.obs["batch"] = sample_info["geo_id"]
        adata.obs["condition"] = sample_info["condition"]
        adata.obs["patient"] = sample_info["patient"]
        patient_clean = sample_info["patient"].replace("-", "")
        adata.obs["author_cell_id"] = [f"{bc}.{patient_clean}" for bc in adata.obs_names]

        adatas[sample_id] = adata

    logger.info(f"Loaded {len(adatas)} samples.")
    return adatas


def merge_samples(adatas: dict) -> ad.AnnData:
    """
    Merge multiple AnnData objects into a single combined object.

    Handles gene intersection across samples and concatenation with
    proper batch tracking in .obs['id'].

    Parameters
    ----------
    adatas : dict
        Mapping of sample_id -> AnnData.

    Returns
    -------
    ad.AnnData
        Merged AnnData with all samples concatenated.
    """
    adata_list = list(adatas.values())
    sample_ids = list(adatas.keys())

    logger.info(f"Merging {len(adata_list)} samples...")
    merged = ad.concat(
        adata_list,
        join="inner",
        label="id",
        keys=sample_ids,
        index_unique="-",
    )
    merged.obs_names_make_unique()

    # Create a clean sample identifier
    merged.obs["id"] = merged.obs["id"].astype("category")

    # Extract patient ID for downstream use
    merged.obs["pid"] = merged.obs["patient"].astype(str)

    n_cells = merged.n_obs
    n_genes = merged.n_vars
    logger.info(f"Merged object: {n_cells:,} cells x {n_genes:,} genes")
    logger.info(f"  Cells per sample:\n{merged.obs['id'].value_counts().to_string()}")

    return merged


def attach_clinical_metadata(adata: ad.AnnData, config: dict) -> ad.AnnData:
    """
    Attach clinical metadata from Excel files provided by the study authors.

    Parameters
    ----------
    adata : ad.AnnData
        Merged AnnData object.
    config : dict
        Pipeline configuration with metadata file paths.

    Returns
    -------
    ad.AnnData
        AnnData with clinical metadata columns added to .obs.
    """
    data_dir = Path(config["paths"]["data_dir"])

    # Load cell-level metadata if available
    cell_meta_file = data_dir / config["dataset"]["cell_meta"]
    if cell_meta_file.exists():
        logger.info(f"Loading cell metadata from {cell_meta_file.name}")
        cell_meta = pd.read_excel(cell_meta_file)
        cell_meta = cell_meta.set_index("Cell ID")

        # Match via author_cell_id or index
        match_col = "author_cell_id" if "author_cell_id" in adata.obs.columns else None
        if match_col:
            matched_meta = cell_meta.reindex(adata.obs[match_col].values)
            n_matched = (~matched_meta["Cell Type"].isna()).sum()
            logger.info(f"  Matched {n_matched:,} / {adata.n_obs:,} cells with author annotations")
            for col in cell_meta.columns:
                if col not in adata.obs.columns:
                    adata.obs[col] = matched_meta[col].values
        else:
            common_cells = list(set(cell_meta.index) & set(adata.obs.index))
            if common_cells:
                logger.info(f"  Matched {len(common_cells):,} cells with metadata")
                for col in cell_meta.columns:
                    if col not in adata.obs.columns:
                        adata.obs[col] = cell_meta.reindex(adata.obs.index)[col]
            else:
                logger.warning("  No cell barcodes matched. Barcode format may differ.")

    # Load patient-level clinical metadata
    clinical_meta_file = data_dir / config["dataset"]["clinical_meta"]
    if clinical_meta_file.exists():
        logger.info(f"Loading clinical metadata from {clinical_meta_file.name}")
        clinical_meta = pd.read_excel(clinical_meta_file)
        logger.info(f"  {len(clinical_meta)} patients in clinical table")

    return adata


def run(config: dict) -> ad.AnnData:
    """
    Execute the full data loading pipeline.

    Steps:
        1. Download data from GEO (if not cached)
        2. Load .h5 count matrices
        3. Merge samples into single AnnData
        4. Attach clinical metadata

    Parameters
    ----------
    config : dict
        Full pipeline configuration dictionary.

    Returns
    -------
    ad.AnnData
        Merged, metadata-annotated AnnData ready for QC.
    """
    data_dir = Path(config["paths"]["data_dir"])
    checkpoint = Path(config["paths"]["checkpoints_dir"]) / "01_raw_merged.h5ad"

    if checkpoint.exists():
        logger.info(f"Loading cached merged data from {checkpoint}")
        return sc.read_h5ad(str(checkpoint))

    # Step 1: Download
    if not any(data_dir.glob("*.h5")):
        download_geo_data(config)
    else:
        logger.info("Data files already present, skipping download.")

    # Step 2: Load
    adatas = load_h5_matrices(config)
    if not adatas:
        raise FileNotFoundError(
            f"No .h5 files found in {data_dir}. "
            "Please download the data manually from GEO: GSE207633"
        )

    # Step 3: Merge
    adata = merge_samples(adatas)

    # Step 4: Metadata
    adata = attach_clinical_metadata(adata, config)

    # Save checkpoint
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(checkpoint))
    logger.info(f"Saved checkpoint: {checkpoint}")

    return adata
