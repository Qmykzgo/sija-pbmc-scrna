# 🧬 Single-Cell RNA-seq Analysis of SJIA PBMCs

[![Python](https://img.shields.io/badge/Python-3.10+-3776AB?logo=python&logoColor=white)](https://python.org)
[![Scanpy](https://img.shields.io/badge/Scanpy-1.10+-blue)](https://scanpy.readthedocs.io)
[![scVI](https://img.shields.io/badge/scVI--tools-1.1+-green)](https://scvi-tools.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Dissecting immune cell heterogeneity in Systemic Juvenile Idiopathic Arthritis using single-cell transcriptomics**

A reproducible scRNA-seq analysis pipeline for PBMC samples from pediatric SJIA patients (GEO: [GSE207633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207633)). Demonstrates a complete workflow from raw count matrices to biological interpretation, including batch integration benchmarking, automated cell annotation, compositional analysis, and pathway activity inference.

---

## 📊 Pipeline Overview

```mermaid
graph LR
    A[📥 Data Loading<br/>GEO Download] --> B[🔍 Quality Control<br/>MT%, Doublets]
    B --> C[⚙️ Preprocessing<br/>Normalize, HVG, PCA]
    C --> D[🔗 Integration<br/>BBKNN · Harmony · scVI]
    D --> E[📍 Clustering<br/>Leiden · UMAP]
    E --> F[🏷️ Annotation<br/>CellTypist · Markers]
    F --> G[📈 Composition<br/>scCODA · Milo]
    G --> H[🧪 Pathways<br/>Decoupler · PROGENy]

    style A fill:#1a1a2e,color:#e94560,stroke:#e94560
    style B fill:#1a1a2e,color:#0f3460,stroke:#0f3460
    style C fill:#1a1a2e,color:#16213e,stroke:#16213e
    style D fill:#1a1a2e,color:#533483,stroke:#533483
    style E fill:#1a1a2e,color:#e94560,stroke:#e94560
    style F fill:#1a1a2e,color:#0f3460,stroke:#0f3460
    style G fill:#1a1a2e,color:#16213e,stroke:#16213e
    style H fill:#1a1a2e,color:#533483,stroke:#533483
```

## 🗂️ Project Structure

```
sjia-pbmc-scrna/
├── README.md                 # This file
├── LICENSE                   # MIT License
├── environment.yml           # Conda environment specification
├── config.yaml               # Pipeline parameters & thresholds
├── pipeline.py               # Main entry point (run all steps)
├── src/
│   ├── __init__.py
│   ├── data_loader.py        # GEO data download & loading
│   ├── qc.py                 # Quality control & doublet detection
│   ├── preprocessing.py      # Normalization, HVG selection, scaling
│   ├── integration.py        # Batch correction benchmarking
│   ├── clustering.py         # Dimensionality reduction & clustering
│   ├── annotation.py         # Automated cell type annotation
│   ├── differential.py       # Compositional & differential abundance
│   ├── pathway.py            # Pathway activity scoring
│   └── visualization.py      # Publication-quality figure generation
├── notebooks/
│   └── 01_full_analysis.py   # Interactive walkthrough (Jupyter-ready)
├── results/
│   └── figures/              # Generated plots
└── data/                     # Downloaded data (gitignored)
```

## 🚀 Quick Start

### 1. Clone & Setup Environment

```bash
git clone https://github.com/YOUR_USERNAME/sjia-pbmc-scrna.git
cd sjia-pbmc-scrna

conda env create -f environment.yml
conda activate sjia-scrna
```

### 2. Run the Full Pipeline

```bash
python pipeline.py
```

Or run individual steps:

```bash
python pipeline.py --step qc
python pipeline.py --step integration
python pipeline.py --step annotation
```

### 3. Explore Interactively

```bash
jupyter notebook notebooks/01_full_analysis.py
```

## 📋 Dataset

| Property | Details |
|----------|---------|
| **GEO Accession** | [GSE207633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207633) |
| **Disease** | Systemic Juvenile Idiopathic Arthritis (SJIA) |
| **Tissue** | Peripheral Blood Mononuclear Cells (PBMC) |
| **Patients** | 21 pediatric (4 subtypes + 5 controls) |
| **Technology** | 10x Genomics Chromium |
| **Subset Used** | 4 samples (2 Healthy, 2 MAS) |

## 🔬 Methods Summary

### Quality Control
- Mitochondrial gene percentage filtering (<25%)
- Ribosomal gene annotation (RPS/RPL)
- Doublet detection via **Scrublet** (threshold=0.25)
- Cell/gene count thresholds (500-25000 UMI, 700-4000 genes)

### Batch Integration Benchmarking
Three methods compared using **scib** metrics:

| Method | Approach | Key Metric |
|--------|----------|------------|
| **BBKNN** | Graph-based | Batch ASW |
| **Harmony** | Linear PCA correction | Bio conservation |
| **scVI** | Variational autoencoder | Overall score |

### Cell Type Annotation
- **CellTypist** automated annotation (`Immune_All_Low.pkl` model)
- Majority voting across over-clustering for robustness
- Marker gene validation (CD3E, CD8A, CD14, CD19, FOXP3, FCGR3A)

### Compositional Analysis
- **scCODA**: Bayesian Dirichlet-Multinomial model for cell type proportion changes
- **Milo**: KNN-based differential abundance testing across disease subtypes

### Pathway Activity
- **Decoupler + PROGENy**: Infer pathway activity scores per cell
- Differential pathway activity across SJIA clinical subtypes

## ⚙️ Configuration

All pipeline parameters are centralized in [`config.yaml`](config.yaml):

```yaml
qc:
  min_genes: 200
  min_cells: 3
  max_mt_pct: 25
  min_counts: 500
  max_counts: 25000

integration:
  methods: ["bbknn", "harmony", "scvi"]
  n_hvg: 3000
  n_pcs: 50

clustering:
  resolution: 1.0
  n_neighbors: 15
```

## 📄 License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

## 📚 References

- Schulert, G.S. et al. (2022). *Single-cell transcriptomic analysis of SJIA PBMCs*. GEO: GSE207633
- Wolf, F.A. et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*
- Gayoso, A. et al. (2022). A Python library for probabilistic analysis of single-cell omics data. *Nature Biotechnology*
- Domínguez Conde, C. et al. (2022). Cross-tissue immune cell analysis reveals tissue-specific features in humans. *Science*
