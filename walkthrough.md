# Walkthrough: SJIA PBMC scRNA-seq Portfolio Project

## What Was Built

A complete, GitHub-ready single-cell RNA-seq analysis pipeline at **`D:\sjia-pbmc-scrna\`** analyzing PBMCs from pediatric SJIA patients (GEO: GSE207633).

## Project Structure

```
D:\sjia-pbmc-scrna\
├── README.md                  ← Professional README with badges & mermaid diagram
├── LICENSE                    ← MIT
├── .gitignore                 ← Ignores data/, checkpoints, __pycache__
├── environment.yml            ← Conda env with all dependencies
├── config.yaml                ← All pipeline parameters centralized
├── pipeline.py                ← CLI entry point with step selection
├── src/
│   ├── __init__.py
│   ├── data_loader.py         ← GEO download, .h5 loading, merging
│   ├── qc.py                  ← MT%, Scrublet, filtering, QC plots
│   ├── preprocessing.py       ← normalize → log1p → HVG → PCA
│   ├── integration.py         ← BBKNN vs Harmony vs scVI + scib benchmark
│   ├── clustering.py          ← KNN → Leiden → UMAP
│   ├── annotation.py          ← CellTypist + marker genes + dotplots
│   ├── differential.py        ← scCODA + Milo compositional analysis
│   ├── pathway.py             ← Decoupler + PROGENy pathway scoring
│   └── visualization.py       ← Publication-quality summary figure
├── notebooks/
│   └── 01_full_analysis.py    ← Interactive jupytext notebook
├── data/                      ← Downloaded data (gitignored)
└── results/
    ├── figures/               ← Output plots
    └── checkpoints/           ← .h5ad intermediates
```

## How to Run

### 1. Create the conda environment
```bash
cd D:\sjia-pbmc-scrna
conda env create -f environment.yml
conda activate sjia-scrna
```

### 2. Run the full pipeline
```bash
python pipeline.py
```

### 3. Or run individual steps
```bash
python pipeline.py --step qc --only         # QC only
python pipeline.py --step integration        # From integration onward
python pipeline.py --step annotation --only  # Annotation only
```

### 4. Interactive exploration
```bash
jupyter notebook notebooks/01_full_analysis.py
```

## Pipeline Steps (8-step workflow)

| # | Step | Module | What It Does |
|---|------|--------|-------------|
| 1 | Data Loading | `data_loader.py` | Downloads from GEO, loads .h5 matrices, merges 4 samples |
| 2 | QC | `qc.py` | MT% annotation, Scrublet doublets, threshold filtering |
| 3 | Preprocessing | `preprocessing.py` | Normalize (10k), log1p, 3000 HVGs, PCA (50 PCs) |
| 4 | Integration | `integration.py` | BBKNN vs Harmony vs scVI, scib benchmarking |
| 5 | Clustering | `clustering.py` | KNN graph → Leiden → UMAP embedding |
| 6 | Annotation | `annotation.py` | CellTypist (Immune_All_Low), marker gene validation |
| 7 | Composition | `differential.py` | scCODA Bayesian + Milo KNN differential abundance |
| 8 | Pathways | `pathway.py` | Decoupler + PROGENy (14 pathways per cell) |

## Key Design Decisions

- **Checkpointing**: Each step saves an `.h5ad` checkpoint so you can resume from any point
- **Config-driven**: All thresholds in `config.yaml` — change once, apply everywhere
- **Modular**: Each `src/*.py` module has a `run()` entry point and works independently
- **Portfolio-ready**: Professional README with badges, mermaid workflow diagram, methods summary, and references
