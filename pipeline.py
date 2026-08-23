#!/usr/bin/env python3
"""
SJIA PBMC scRNA-seq Analysis Pipeline
======================================
Main entry point for the single-cell RNA-seq analysis pipeline.

Analyzes PBMC samples from pediatric SJIA patients (GSE207633) through
a complete workflow: QC → Integration → Clustering → Annotation →
Compositional Analysis → Pathway Activity.

Usage:
    python pipeline.py                    # Run all steps
    python pipeline.py --step qc          # Run specific step
    python pipeline.py --step integration # Run up to integration
    python pipeline.py --resume           # Resume from last checkpoint

Author: Yerkanat
"""

import argparse
import logging
import os
import sys
import time
from pathlib import Path

import yaml

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"


if sys.platform == "win32":
    try:
        if hasattr(sys.stdout, "reconfigure"):
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        if hasattr(sys.stderr, "reconfigure"):
            sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


# ─── Configure logging ───────────────────────────────────────────────
def setup_logging(log_dir: Path) -> None:
    """Configure logging to both console and file."""
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "pipeline.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)-7s | %(name)-20s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file, mode="a", encoding="utf-8"),
        ],
    )


# ─── Load config ─────────────────────────────────────────────────────
def load_config(config_path: str = "config.yaml") -> dict:
    """Load pipeline configuration from YAML file."""
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    return config


# ─── Pipeline steps ──────────────────────────────────────────────────
STEPS = [
    "data_loading",
    "qc",
    "preprocessing",
    "integration",
    "clustering",
    "annotation",
    "differential",
    "pathway",
    "visualization",
]


def run_pipeline(config: dict, start_step: str = None, single_step: bool = False):
    """
    Execute the scRNA-seq analysis pipeline.

    Parameters
    ----------
    config : dict
        Pipeline configuration.
    start_step : str, optional
        Step to start from (or run individually).
    single_step : bool
        If True, run only the specified step.
    """
    logger = logging.getLogger("pipeline")

    # Determine which steps to run
    if start_step:
        if start_step not in STEPS:
            logger.error(f"Unknown step: {start_step}. Available: {STEPS}")
            sys.exit(1)
        start_idx = STEPS.index(start_step)
        if single_step:
            steps_to_run = [start_step]
        else:
            steps_to_run = STEPS[start_idx:]
    else:
        steps_to_run = STEPS

    logger.info("=" * 70)
    logger.info("  SJIA PBMC scRNA-seq Analysis Pipeline")
    logger.info("=" * 70)
    logger.info(f"  Steps to run: {' -> '.join(steps_to_run)}")
    logger.info("=" * 70)

    # Create output directories
    for dir_key in ["results_dir", "figures_dir", "checkpoints_dir"]:
        Path(config["paths"][dir_key]).mkdir(parents=True, exist_ok=True)

    adata = None
    total_start = time.time()

    for step in steps_to_run:
        step_start = time.time()
        logger.info(f"\n{'-' * 60}")
        logger.info(f"  Step: {step.upper()}")
        logger.info(f"{'-' * 60}")

        if step == "data_loading":
            from src import data_loader
            adata = data_loader.run(config)

        elif step == "qc":
            if adata is None:
                from src import data_loader
                adata = data_loader.run(config)
            from src import qc
            adata = qc.run(adata, config)

        elif step == "preprocessing":
            if adata is None:
                adata = _load_latest_checkpoint(config, step)
            from src import preprocessing
            adata = preprocessing.run(adata, config)

        elif step == "integration":
            if adata is None:
                adata = _load_latest_checkpoint(config, step)
            from src import integration
            adata = integration.run(adata, config)

        elif step == "clustering":
            if adata is None:
                adata = _load_latest_checkpoint(config, step)
            from src import clustering
            adata = clustering.run(adata, config)

        elif step == "annotation":
            if adata is None:
                adata = _load_latest_checkpoint(config, step)
            from src import annotation
            adata = annotation.run(adata, config)

        elif step == "differential":
            if adata is None:
                adata = _load_latest_checkpoint(config, step)
            from src import differential
            adata = differential.run(adata, config)

        elif step == "pathway":
            if adata is None:
                adata = _load_latest_checkpoint(config, step)
            from src import pathway
            adata = pathway.run(adata, config)

        elif step == "visualization":
            if adata is None:
                adata = _load_latest_checkpoint(config, step)
            from src import visualization
            visualization.run(adata, config)

        elapsed = time.time() - step_start
        logger.info(f"  Step '{step}' completed in {elapsed:.1f}s")

    # Final summary
    total_elapsed = time.time() - total_start
    logger.info(f"\n{'=' * 70}")
    logger.info(f"  Pipeline complete! Total time: {total_elapsed/60:.1f} minutes")
    logger.info(f"  Results: {config['paths']['results_dir']}")
    logger.info(f"  Figures: {config['paths']['figures_dir']}")
    logger.info(f"{'=' * 70}")


def _load_latest_checkpoint(config: dict, target_step: str):
    """Load the required prerequisite checkpoint for the target step."""
    import scanpy as sc

    checkpoint_dir = Path(config["paths"]["checkpoints_dir"])

    # Checkpoint mapping (step -> prerequisite checkpoint file)
    checkpoints = {
        "qc": "01_raw_merged.h5ad",
        "preprocessing": "02_qc_filtered.h5ad",
        "integration": "03_preprocessed.h5ad",
        "clustering": "04_integrated.h5ad",
        "annotation": "05_clustered.h5ad",
        "differential": "06_annotated.h5ad",
        "pathway": "07_differential.h5ad",
        "visualization": "08_pathway.h5ad",
    }

    logger = logging.getLogger("pipeline")

    if target_step in checkpoints:
        cp_file = checkpoint_dir / checkpoints[target_step]
        if cp_file.exists():
            logger.info(f"Resuming from checkpoint: {cp_file}")
            return sc.read_h5ad(str(cp_file))
        else:
            raise FileNotFoundError(
                f"Prerequisite checkpoint '{checkpoints[target_step]}' for step '{target_step}' was not found in {checkpoint_dir}. "
                f"Please run earlier steps first to generate this checkpoint."
            )

    raise FileNotFoundError(
        f"No checkpoint definition found for step '{target_step}'."
    )


# ─── CLI ──────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="SJIA PBMC scRNA-seq Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python pipeline.py                        # Run full pipeline
    python pipeline.py --step qc              # Run only QC step
    python pipeline.py --step integration     # Run from integration onward
    python pipeline.py --step annotation --only  # Run annotation step only
    python pipeline.py --config my_config.yaml   # Use custom config
        """,
    )
    parser.add_argument(
        "--step",
        choices=STEPS,
        help="Start from (or run only) this pipeline step.",
    )
    parser.add_argument(
        "--only",
        action="store_true",
        help="Run only the specified --step, not subsequent steps.",
    )
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to configuration file (default: config.yaml).",
    )

    args = parser.parse_args()

    config = load_config(args.config)
    setup_logging(Path(config["paths"]["results_dir"]))

    run_pipeline(config, start_step=args.step, single_step=args.only)


if __name__ == "__main__":
    main()
