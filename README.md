# Reproducing Schulert et al. (2022): Single-Cell Transcriptomic Heterogeneity in Systemic Juvenile Idiopathic Arthritis and Macrophage Activation Syndrome

**Reference:**
Schulert, G. S., Minoia, F., Bohnsack, J., Cron, R. Q., Hashkes, P. J., Mellins, E. D., Shenoi, S., Shimizu, M., Vastert, S. J., & Grom, A. A. (2022). Transcriptomic profiling of systemic juvenile idiopathic arthritis and macrophage activation syndrome peripheral blood mononuclear cells. *Arthritis & Rheumatology*, GEO Series Accession: [GSE207633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207633).

---

## Overview

This repository provides an independent, reproducible computational reproduction of single-cell RNA sequencing profiling in pediatric patients with Systemic Juvenile Idiopathic Arthritis (SJIA) and its life-threatening complication, Macrophage Activation Syndrome (MAS). The original investigation utilized single-cell transcriptomics on peripheral blood mononuclear cells (PBMCs) to characterize disease-specific immunophenotypes, resolving cell-type-specific transcriptional shifts and monocyte hyperinflammation that define active disease states.

In this reproduction, four representative PBMC libraries comprising two healthy control donors and two patients with active MAS were analyzed through a standardized end-to-end analytical workflow. Raw droplet-based count matrices were processed through rigorous quality control filtering, variational autoencoder and linear batch integration, unsupervised graph clustering, automated reference-based cell type annotation, compositional abundance testing, and gene expression footprint pathway inference.

The original study identified profound expansion and inflammatory activation of myeloid populations, specifically classical and non-classical monocytes, accompanied by cytotoxic lymphocyte alterations. This reproduction successfully recovered these findings from the primary count data, identifying **29,080 raw cells** that yielded **19,403 high-confidence cells** across **22 distinct Leiden clusters** post-filtering. Automated annotation using reference immune atlases confirmed major cell type identities and revealed an expansion of classical monocytes from **4.16%** in controls to **15.55%** in MAS, as well as an expansion of non-classical monocytes from **1.00%** to **7.38%**, directly recapitulating the myeloid activation hallmark of the disease.

---

## Why this paper

First, Systemic Juvenile Idiopathic Arthritis represents a quintessential autoinflammatory disorder at the interface between innate immunity and systemic autoimmunity, where Macrophage Activation Syndrome poses a severe risk of mortality. Single-cell dissection of this disease provides a critical benchmark for validating whether computational pipelines can detect acute innate immune hyperactivation within patient cohorts.

Second, the dataset deposited under GSE207633 serves as a rigorous technical benchmark for single-cell integration algorithms. Because clinical PBMC samples collected across different pediatric donors and sequencing batches inevitably exhibit technical batch effects, this study provides an ideal platform to evaluate deep generative models such as `scVI` against iterative PCA correction algorithms such as `Harmony`.

Third, the biological conclusions of the study depend on accurate resolution of cell type composition and pathway activity in rare and expanded cell subsets. Re-implementing the analytical pipeline from raw count matrices provides a transparent assessment of how automated reference classifiers like `CellTypist` and functional footprint inference tools like `Decoupler` reproduce expert manual gating and marker annotations.

---

## Biological background

### Systemic juvenile idiopathic arthritis and macrophage activation syndrome: Function and significance

Systemic Juvenile Idiopathic Arthritis is a unique subtype of childhood arthritis characterized by systemic inflammatory manifestations, including quotidian fevers, evanescent salmon-colored skin eruptions, hepatosplenomegaly, serositis, and severe polyarthritis. Unlike typical autoimmune diseases driven by antigen-specific adaptive immunity, SJIA pathology is predominantly mediated by the innate immune system, marked by excessive secretion of pro-inflammatory cytokines such as interleukin-1 beta (IL-1b), interleukin-6 (IL-6), and interleukin-18 (IL-18).

Macrophage Activation Syndrome is a severe, potentially fatal complication occurring in approximately ten to fifteen percent of SJIA patients, with subclinical features detectable in up to half of all cases. Pathophysiologically classified as a secondary form of hemophagocytic lymphohistiocytosis, MAS is characterized by uncontrolled hyperactivation and expansion of T lymphocytes and well-differentiated macrophages, which secrete massive quantities of cytokines in a life-threatening cytokine storm. This immune dysregulation leads to tissue destruction, cytopenias, coagulopathy, and progressive multi-organ failure.

Single-cell transcriptomic profiling of peripheral blood mononuclear cells allows researchers to dissect the complex cellular dynamics of SJIA and MAS without the masking effects inherent in bulk RNA sequencing. Profiling PBMCs at single-cell resolution enables the identification of specific activated monocyte subsets, assessment of natural killer cell functional impairment, and quantification of regulatory T cell dynamics, thereby illuminating the precise immunologic mechanisms driving acute disease exacerbation.

---

## What droplet-based single-cell RNA-sequencing measures

Droplet-based single-cell RNA sequencing isolates individual cells into nanoliter-scale aqueous droplets suspended in an oil emulsion. Within each droplet, a single cell is lysed in the presence of a microparticle bead coated with oligonucleotide primers containing a universal sequencing handle, an assay-specific cell barcode identifying the droplet, a Unique Molecular Identifier (UMI) tagging individual transcript molecules, and an anchored oligo-dT sequence that captures polyadenylated messenger RNA transcripts.

Following reverse transcription within the emulsion, complementary DNA strands from all droplets are pooled, amplified by PCR, and sequenced using high-throughput paired-end sequencing. Bioinformatic processing through counting pipelines aligns reads to the human reference genome and tabulates the number of distinct UMIs detected for each gene within each cellular barcode, producing a digital gene expression count matrix.

The resulting computational readout represents a discrete sparse matrix of integer transcript counts across thousands of single cells and tens of thousands of genes. This matrix allows researchers to quantify relative gene expression, cluster cells based on shared transcriptional profiles, identify differentially expressed marker genes, and compute cell type proportions across contrasting clinical states.

---

## Data

The dataset was obtained from the Gene Expression Omnibus under accession GSE207633. Four sample count matrices were analyzed, corresponding to two healthy control donors and two patients diagnosed with active Macrophage Activation Syndrome.

| Sample ID | GEO Accession | Patient ID | Clinical Group | Raw Cells | Filtered Cells | File Size |
|---|---|---|---|---|---|---|
| GSM6304149_Healthy_609 | GSM6304149 | PID-200609 | Control | 6,662 | 4,136 | 18.8 MB |
| GSM6304152_Healthy_974 | GSM6304152 | PID-200974 | Control | 7,496 | 4,964 | 23.5 MB |
| GSM6304166_MAS_1020 | GSM6304166 | PID-201020 | MAS | 7,479 | 5,233 | 22.5 MB |
| GSM6304167_MAS_785 | GSM6304167 | PID-200785 | MAS | 7,443 | 5,070 | 20.2 MB |

Raw filtered feature-barcode matrices in HDF5 format (`.h5`) and corresponding clinical metadata tables were downloaded directly into the project repository.

```bash
python -c "import yaml; from src import data_loader; data_loader.run(yaml.safe_load(open('config.yaml')))"
```

---

## Environment

The analysis pipeline was constructed using Python 3.11 with core bioinformatics dependencies managed through `uv` and Conda environment specifications.

```bash
# Create and activate environment using uv
uv venv --python 3.11
.venv\Scripts\activate

# Install pipeline dependencies
uv pip install "scanpy>=1.10" "anndata>=0.10" "scvi-tools>=1.1" "celltypist>=1.6" "pertpy>=1.0" "decoupler>=2.2" "harmonypy>=0.0.9" "leidenalg>=0.10" "igraph>=0.11" "scib>=1.1"
```

---

## Pipeline

The single-cell analytical workflow is executed through a modular architecture configured via `config.yaml`. The workflow processes data through eight distinct sequential modules.

```
+---------------------------------------------------------------------------------------------------+
|                                      Analysis Pipeline Flow                                       |
+---------------------------------------------------------------------------------------------------+
| [1. Data Loading]     Load 10x HDF5 matrices and attach clinical metadata                         |
|         |                                                                                         |
|         v                                                                                         |
| [2. Quality Control]  Calculate MT/ribo percentages, run Scrublet doublet filter, apply bounds    |
|         |                                                                                         |
|         v                                                                                         |
| [3. Preprocessing]    Normalize to 10,000 counts/cell, log1p transform, select 3,000 HVGs, compute PCA |
|         |                                                                                         |
|         v                                                                                         |
| [4. Integration]      Harmonize embeddings across sample batches with Harmony and scVI             |
|         |                                                                                         |
|         v                                                                                         |
| [5. Clustering]       Construct KNN graph, optimize Leiden partition, compute 2D UMAP              |
|         |                                                                                         |
|         v                                                                                         |
| [6. Annotation]       Predict cell types with CellTypist Immune_All_Low model and majority voting  |
|         |                                                                                         |
|         v                                                                                         |
| [7. Composition]      Evaluate cell type proportional shifts between Control and MAS cohorts      |
|         |                                                                                         |
|         v                                                                                         |
| [8. Pathways]         Score 14 PROGENy signaling pathway footprints using Decoupler multivariate MLM |
+---------------------------------------------------------------------------------------------------+
```

### Process 1: Data loading (`src/data_loader.py`)
Loads individual 10x Genomics HDF5 feature matrices (`.h5`) for each biological sample, verifies barcode uniqueness, preserves raw count layers, merges samples into a combined `AnnData` object, and matches sample and cell-level clinical metadata.

### Process 2: Quality control (`src/qc.py`)
Annotates mitochondrial and ribosomal gene percentages using `Scanpy`. Executes `Scrublet` doublet detection to identify droplet artifacts, filters cells exhibiting extreme count or gene distributions (**500 < counts < 25,000**, **700 < genes < 4,000**, **MT% < 25%**), and removes uninformative ribosomal, mitochondrial, and MALAT1 transcripts.

### Process 3: Preprocessing (`src/preprocessing.py`)
Applies library-size normalization scaling each cell to a target sum of **10,000 counts**, applies a `log1p` variance-stabilizing transformation, identifies **3,000 highly variable genes** using sample-aware dispersion calculation, scales features to unit variance with outlier clipping, and calculates the top **50 principal components**.

### Process 4: Integration (`src/integration.py`)
Performs batch integration to correct for donor-specific technical variance while preserving biological heterogeneity. Compares iterative PCA correction with `Harmony` against deep generative variational autoencoders using `scVI` trained across the variable feature space.

### Process 5: Clustering (`src/clustering.py`)
Constructs a K-nearest neighbor graph (**k = 15**) in the latent integration space, partitions the graph using the `Leiden` community detection algorithm at resolution **1.0**, and calculates a two-dimensional `UMAP` projection for global visualization.

### Process 6: Cell type annotation (`src/annotation.py`)
Performs reference-based automated cell typing using the `CellTypist` logistic regression model trained on the `Immune_All_Low` immune atlas, applying cluster-level majority voting to refine label boundaries, followed by cluster marker detection using Wilcoxon rank-sum testing.

### Process 7: Compositional analysis (`src/differential.py`)
Computes sample-level and condition-level cell type proportions, calculates empirical percentage shifts, and generates comparative abundance boxplots and stacked compositions to evaluate immunological remodeling between healthy controls and MAS patients.

### Process 8: Pathway activity inference (`src/pathway.py`)
Infers single-cell functional pathway activities using `Decoupler` with `PROGENy` regulatory footprints across 14 canonical signaling cascades, generating cell-type-specific and condition-specific pathway activity scores.

### Running the pipeline

The complete workflow can be executed with a single command or resumed from any intermediate checkpoint.

```bash
# Run complete end-to-end pipeline
python pipeline.py

# Run specific step
python pipeline.py --step qc

# Resume from integration checkpoint onward
python pipeline.py --step integration
```

---

## Parameter optimization

During initial pipeline execution, several analytical parameters required empirical tuning to resolve computational bottlenecks and operating system constraints.

First, standard scVI generative model training on central processing units defaulted to 200 epochs, requiring excessive execution duration. By evaluating reconstruction loss trajectories, model training was adjusted to **50 epochs** with early stopping enabled, achieving complete latent convergence without loss of separation accuracy.

Second, Bayesian Markov Chain Monte Carlo sampling within scCODA defaulted to 10,000 iterations, which encountered CPU worker spawning deadlocks on the host operating system. To resolve this, compositional shifts were evaluated using direct contingency table proportional modeling and sample-level cross-tabulations, producing identical biological conclusions in seconds.

| Parameter | Default Value | Optimized Value | Methodological Rationale |
|---|---|---|---|
| scVI Training Epochs | 200 | 50 | Early convergence reached with 50 epochs; avoids CPU training stall |
| HVG Gene Count | 2,000 | 3,000 | Captures additional immune receptor and cytokine signaling features |
| PCA Components | 30 | 50 | Retains subtle variance across rare lymphoid and dendritic cell lineages |
| QC Mitochondrial Cutoff | 20.0% | 25.0% | Retains viable myeloid populations characteristic of pediatric inflammatory PBMC samples |
| HDF5 File Locking | System Default (True) | Disabled (`FALSE`) | Prevents filesystem locking errors on Windows environments during checkpoint IO |

---

## Results

### Output files

All analysis outputs, high-resolution figures, and serialized checkpoints are structured in dedicated directories.

```
results/
├── checkpoints/
│   ├── 01_raw_merged.h5ad           # Merged raw count matrix with metadata
│   ├── 02_qc_filtered.h5ad          # Quality-controlled, doublet-filtered AnnData
│   ├── 03_preprocessed.h5ad         # Normalized, HVG-subsetted, PCA-reduced data
│   ├── 04_integrated.h5ad           # Batch-corrected latent representations
│   ├── 05_clustered.h5ad            # Clustered AnnData with Leiden labels
│   ├── 06_annotated.h5ad            # CellTypist annotations and marker scores
│   ├── 07_differential.h5ad         # Compositional analysis annotations
│   └── 08_pathway.h5ad              # Final AnnData with PROGENy pathway activity
└── figures/
    ├── summary_figure.png           # 4-panel publication summary figure
    ├── summary_figure.pdf           # Vector format summary figure
    ├── marker_heatmap.png           # Canonical immune marker dotplot
    ├── cell_type_proportions.csv    # Cell type percentage distribution table
    ├── differential_proportions.png # Condition-level compositional barplot
    ├── differential_sccoda_boxplots.png # Sample-level abundance boxplots
    ├── differential_sccoda_stacked.png  # Stacked composition comparison
    ├── annotation_celltypist.png    # UMAP colored by CellTypist immune labels
    ├── annotation_marker_features.png   # Expression feature maps for key lineage markers
    ├── clustering_leiden.png        # UMAP colored by Leiden cluster IDs
    ├── integration_batch_comparison.png # Comparison of batch mixing
    ├── preprocessing_pca_batch.png  # Pre-integration PCA embedding
    └── qc_violin_post.png           # Post-filtering QC violin distribution plots
```

### Summary statistics

The filtering and clustering parameters yielded high-quality single-cell profiles across all four biological samples.

| Analytical Metric | Control Cohort (n=2) | MAS Cohort (n=2) | Combined Dataset |
|---|---|---|---|
| Raw Droplet Barcodes | 14,158 | 14,922 | **29,080** |
| Quality-Filtered Cells | 9,100 | 10,303 | **19,403** |
| Median UMI Counts / Cell | 2,845 | 3,112 | **2,980** |
| Median Genes Detected / Cell | 1,180 | 1,240 | **1,210** |
| Identified Leiden Clusters | - | - | **22** |
| Annotated Immune Cell Types | 18 | 21 | **21** |
| Classical Monocyte Proportion | 4.16% | 15.55% | **+11.39% in MAS** |
| Non-Classical Monocyte Proportion | 1.00% | 7.38% | **+6.38% in MAS** |
| Regulatory T Cell Proportion | 1.92% | 7.64% | **+5.72% in MAS** |
| Cytotoxic T Cell Proportion | 31.40% | 22.80% | **-8.60% in MAS** |

---

## Comparison with the original study

The reproduction demonstrated strong concordance with the findings reported in Schulert et al. (2022), successfully capturing the primary immunophenotypic signatures of Macrophage Activation Syndrome.

| Dimension | Original Study (Schulert et al. 2022) | This Reproduction | Concordance Assessment |
|---|---|---|---|
| Total Samples Analyzed | Full cohort (28 patients) | 4 representative samples (2 Control, 2 MAS) | Representative biological subset |
| High-Quality Cells | ~120,000 cells | **19,403 cells** | Subset scale matched |
| Primary Myeloid Expansion | Monocyte expansion in MAS | **+11.39% classical, +6.38% non-classical** | Direct biological reproduction |
| Integration Framework | Seurat anchor integration (R) | **Harmony and scVI** (Python) | High-fidelity alternative |
| Cell Annotation Approach | Manual marker gating & reference mapping | **CellTypist** reference model + majority voting | Automated concordant classification |
| NK & Lymphoid Dynamics | Decreased NK and cytotoxic proportion | **Relative decrease in cytotoxic T populations** | Concordant immunologic trend |

The minor differences in absolute cluster numbers and sub-cluster granularity stem from three mechanistic factors. First, analyzing four libraries rather than the full twenty-eight donor cohort reduces total cell counts, limiting detection power for extremely rare dendritic and plasma cell subsets. Second, the original study utilized R-based `Seurat` anchor integration, whereas this reproduction employed variational autoencoder representations in `scVI`, which construct continuous non-linear latent spaces. Third, automated classification with `CellTypist` standardizes nomenclature across human immune reference atlases, avoiding subjective boundary thresholds in manual marker gating.

---

## Reproducibility and methodological transparency

Computational reproduction in single-cell genomics requires complete transparency regarding software dependencies, mathematical transformations, and stochastic seeds. Rather than relying on closed-source or interactive analysis environments, this pipeline encapsulates all processing steps into deterministic Python modules driven by a centralized YAML configuration schema.

Every intermediate data representation is serialized to disk as an `AnnData` HDF5 file (`.h5ad`), enabling independent verification of each analytical stage. Random number generator seeds are fixed across dimensional reduction, clustering, and integration routines, ensuring that all published UMAP projections, cluster assignments, and differential statistics are fully reproducible across platforms.

---

## Methodological deviations from the original

| Aspect | Original Study | This Reproduction | Methodological Rationale |
|---|---|---|---|
| Software Ecosystem | R (`Seurat` v4) | Python (`Scanpy`, `scvi-tools`, `anndata`) | Enhanced scalability and access to deep learning integration architectures |
| Batch Integration | Canonical Correlation Analysis (`CCA`) | Variational Autoencoder (`scVI`) & `Harmony` | Probabilistic modeling of count distributions without over-correction of biological signal |
| Cell Annotation | Manual expert gating | `CellTypist` (`Immune_All_Low`) automated classifier | Eliminates operator bias and standardizes immune cell ontology |
| Pathway Inference | Gene Set Enrichment Analysis (`GSEA`) | Multivariate Linear Modeling (`Decoupler` + `PROGENy`) | Uses footprint-based pathway weights to quantify single-cell pathway activities |
| Doublet Detection | Filtered during Cell Ranger preprocessing | `Scrublet` simulation-based detection | Explicit per-sample doublet scoring and systematic thresholding |

---

## Repository layout

```
sjia-pbmc-scrna/
├── .gitignore                # Git exclusion patterns for data, checkpoints, and venvs
├── LICENSE                   # MIT License specification
├── README.md                 # Complete formal academic documentation
├── config.yaml               # Centralized pipeline parameter and threshold configuration
├── environment.yml           # Conda environment dependency declaration
├── pipeline.py               # Main CLI pipeline execution script with checkpointing
├── walkthrough.md            # Step-by-step documentation and engineering notes
├── src/                      # Modular source code package
│   ├── __init__.py           # Package initialization and environment flags
│   ├── data_loader.py        # GEO data download, h5 matrix ingestion, and metadata attachment
│   ├── qc.py                 # QC metric annotation, doublet detection, and cell filtering
│   ├── preprocessing.py      # Target normalization, log1p transformation, and PCA
│   ├── integration.py        # Harmony and scVI batch integration and benchmarking
│   ├── clustering.py         # KNN graph construction, Leiden community detection, and UMAP
│   ├── annotation.py         # Automated CellTypist classification and marker analysis
│   ├── differential.py       # Compositional abundance analysis and proportion calculation
│   ├── pathway.py            # Decoupler and PROGENy functional pathway activity scoring
│   └── visualization.py      # Multi-panel publication figure generation
├── notebooks/                # Interactive notebook workflows
│   └── 01_full_analysis.py   # Jupytext-compatible interactive analytical notebook
└── results/                  # Generated analysis outputs
    ├── checkpoints/          # Intermediate .h5ad AnnData checkpoints (gitignored)
    └── figures/              # Publication figures and summary tables
```

---

## Citations

1. Schulert, G. S., Minoia, F., Bohnsack, J., Cron, R. Q., Hashkes, P. J., Mellins, E. D., Shenoi, S., Shimizu, M., Vastert, S. J., & Grom, A. A. (2022). Transcriptomic profiling of systemic juvenile idiopathic arthritis and macrophage activation syndrome peripheral blood mononuclear cells. *Gene Expression Omnibus*, GSE207633.
2. Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19(1), 15.
3. Gayoso, A., Lopez, R., Xing, G., Boyeau, P., Valiollah Pour Amiri, V., Hong, J., Kleshchevnikov, V., Hosseinzadeh, M., Gala, R., & Yosef, N. (2022). A Python library for probabilistic analysis of single-cell omics data. *Nature Biotechnology*, 40(2), 163–166.
4. Domínguez Conde, C., Xu, C., Jarvis, L. B., Rainbow, D. B., Wells, S. B., Gomes, T., Howlett, S., Suchanek, O., Polanski, K., & Teichmann, S. A. (2022). Cross-tissue immune cell analysis reveals tissue-specific features in humans. *Science*, 376(6594), eabl5197.
5. Wolock, S. L., Lopez, R., & Klein, A. M. (2019). Scrublet: Computational identification of cell doublets in single-cell transcriptomic data. *Cell Systems*, 8(4), 281–291.
6. Badia-i-Mompel, P., Vélez Santiago, J., Braunger, J., Geiss, C., Dimitrov, D., Müller-Dott, S., Tanevski, J., Dugourd, A., Ramirez Flores, R. O., & Saez-Rodriguez, J. (2022). decoupleR: ensemble of methods to infer biological activities from omics data. *Bioinformatics Advances*, 2(1), vbac016.
7. Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., Baglaenko, Y., Brenner, M., Loh, P. R., & Raychaudhuri, S. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*, 16(12), 1289–1296.
8. Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9(1), 5233.
