# epitome_tools: ML-Powered Cell Typing and Doublet Detection for Pituitary Single-Cell Data

A Python package for automated cell type annotation and doublet detection in single-cell RNA-seq and ATAC-seq data from pituitary gland studies, powered by XGBoost models trained on the Consensus Pituitary Atlas.

[![DOI](https://zenodo.org/badge/930878390.svg)](https://doi.org/10.5281/zenodo.17154160)

## Overview

Epitome Tools provides:
- **Automated cell type annotation** using pre-trained XGBoost classifiers achieving 90-95% accuracy
- **Pituitary-specific doublet detection** with 86-92% accuracy, outperforming generic methods
- **Support for both RNA and ATAC modalities** with concordant predictions
- **Quality control workflows** optimized for pituitary datasets
- **End-to-end processing pipelines** from raw data to annotated cell atlases

Trained on 1.1+ million cells from 256 biological replicates, these models enable consistent, expert-level analysis across the pituitary research community.

## Installation

```bash
pip install epitome-tools
```

## Quick Start

```python
import scanpy as sc
from epitome_tools.workflow import celltype_doublet_workflow

# Load your preprocessed AnnData (normalized, log-transformed)
adata = sc.read_h5ad("your_data.h5ad")

# Run complete annotation and doublet detection in one line
adata = celltype_doublet_workflow(
    adata, 
    active_assay="sc",  # Options: "sc", "sn", "multi_rna"
    modality="rna",     # Options: "rna", "atac"
    smoothing=True      # KNN-based label smoothing
)

# Access results
print(adata.obs['cell_type_final'])           # Smoothed cell type predictions
print(adata.obs['thresholded_doublet_epitome']) # Boolean: True = singlet
```

## Workflow Module

The `workflow` module provides the main interface for cell typing and doublet detection.

### Cell Type Prediction

Automated annotation using models trained on the Consensus Pituitary Atlas:

```python
from epitome_tools.workflow import cell_type_workflow

adata = cell_type_workflow(
    adata_to_use=adata,
    active_assay="sc",      # Assay type: "sc", "sn", or "multi_rna"
    modality="rna",         # Modality: "rna" or "atac"
    in_place=True,          # Modify adata directly
    nan_or_zero='nan',      # Handle missing features
    smoothing=True          # Apply neighbor-based smoothing
)
```

**Key Features:**
- **90-95% accuracy** on unseen datasets across 5-fold cross-validation
- **11 pituitary cell types**: Stem cells, corticotrophs, melanotrophs, gonadotrophs, somatotrophs, lactotrophs, thyrotrophs, plus endothelial, mesenchymal, immune cells, and pituicytes
- **Optional smoothing**: Uses 10 nearest neighbors in PCA space for robust predictions
- **Automatic validation**: Checks feature compatibility (≥70% overlap) and normalization

**Output columns in `adata.obs`:**
- `predicted_cell_type`: Raw model predictions
- `predicted_cell_type_proba`: Maximum prediction probability
- `cell_type_final`: Smoothed predictions (if `smoothing=True`)

### Doublet Detection

Pituitary-optimized doublet detection outperforming scrublet:

```python
from epitome_tools.workflow import doublet_workflow

adata = doublet_workflow(
    adata_to_use=adata,
    active_assay="sc",
    modality="rna",
    in_place=True
)
```

**Key Features:**
- **86-92% accuracy** with calibrated thresholds (95% TPR at 10% FPR)
- **Cell type-aware**: Trained on all pairwise doublet combinations
- **Superior performance**: Consistently outperforms scrublet, especially on challenging datasets

**Output columns in `adata.obs`:**
- `init_predicted_doublet_epitome`: Raw predictions
- `thresholded_doublet_epitome`: Boolean calls (True = singlet, False = doublet)
- `doublet_score_epitome`: Confidence score (higher = more likely singlet)

**Filtering doublets:**
```python
# Remove predicted doublets
adata_clean = adata[adata.obs['thresholded_doublet_epitome']].copy()
```

### Combined Workflow

Run both analyses sequentially:

```python
from epitome_tools.workflow import celltype_doublet_workflow

adata = celltype_doublet_workflow(
    adata,
    active_assay="sc",
    modality="rna",
    smoothing=True
)
```

## Data Requirements

### Input Format

Your AnnData object should be:

1. **Quality filtered**: Remove low-quality cells and genes
2. **Normalized**: Total counts per cell = 10,000
3. **Log-transformed**: Log1p applied

**Example preprocessing:**

```python
import scanpy as sc

# Basic QC
sc.pp.filter_cells(adata, min_genes=200)

# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Ready for Epitome Tools
adata = celltype_doublet_workflow(adata)
```

### Assay Types

- **`sc`**: Single-cell RNA-seq
- **`sn`**: Single-nucleus RNA-seq
- **`multi_rna`**: Multiome RNA component

### Compatibility Checks

Models automatically validate:
- **Feature coverage**: Warns if <70% of model features present
- **Normalization**: Detects improper log-normalization
- **Cell count**: Warns if >50,000 cells (check filtering)

## Atlas Module

The `atlas` module provides tools for large-scale dataset processing and atlas construction, particularly for pituitary studies.

### Data Acquisition

```python
from epitome_tools.atlas import sra_downloader, fastq_downloader

# Download SRA datasets
sra_downloader(raw_path="/data/raw/", df=metadata_df)

# Extract FASTQ files
fastq_downloader(raw_path="/data/raw/", df=metadata_df)
```

### Quantification

```python
from epitome_tools.atlas import run_kb_count_nac

# kallisto|bustools nascent RNA quantification
run_kb_count_nac(
    nac_path="/data/processed/nac/",
    raw_path="/data/raw/",
    df=metadata_df,
    species="mouse"  # or "human"
)
```

### Quality Control

```python
from epitome_tools.atlas import filtering_cells, qc1, qc2

# Multi-step filtering for atlas-scale data
filtering_cells(df=metadata_df, raw_path="/data/")

# Standard QC workflows
adata = qc1(adata, species="mouse", min_genes=800, min_counts=1000)
adata, qc_fig = qc2(adata, species="mouse", experiment="sc")
```

### Visualization

```python
from epitome_tools.atlas import kneeplot, create_assignment_barchart

# Barcode ranking knee plot
fig, data = kneeplot(adata_filtered_final, adata_filtered_init, adata_raw)

# Cell type distribution
fig = create_assignment_barchart(adata.obs['predicted_cell_type'])
```

### Utility Functions

```python
from epitome_tools.atlas import get_sparsity, shrink_high_values_normalized_df

# Calculate matrix sparsity
sparsity = get_sparsity(adata)

# Adjust ambient RNA profiles
adjusted_df = shrink_high_values_normalized_df(
    df=ambient_profile,
    percentile=50,
    shrinkage_factor=0.5
)
```

## Advanced Usage

### Manual Model Loading

```python
from epitome_tools.celltyping import load_celltype_model
from epitome_tools.doublets import load_doublet_model

# Load cell typing model
model, encoder, features = load_celltype_model(
    model_path="path/to/model.json",
    label_encoder_path="path/to/encoder.pkl"
)

# Load doublet model
model, encoder, threshold, features = load_doublet_model(
    model_path="path/to/model.json",
    label_encoder_path="path/to/encoder.pkl",
    threshold_path="path/to/threshold.pkl"
)
```

### Custom Matrix Preparation

```python
from epitome_tools.celltyping import prepare_matrix_celltype
from epitome_tools.doublets import prepare_matrix_doublet

# Prepare for cell typing
X = prepare_matrix_celltype(
    adata, 
    feature_names, 
    active_assay="sn",
    nan_or_zero='zero'
)

# Prepare for doublet detection
X = prepare_matrix_doublet(
    adata,
    feature_names,
    active_assay="sc"
)
```

## Model Performance

### Cell Type Model
- **Cross-validation accuracy**: 90-95% on unseen studies
- **Most cell types**: ~95% recall
- **Thyrotrophs**: 76% recall (lower due to rarity)
- **RNA-ATAC concordance**: 89.2% in multiome data

### Doublet Model
- **Cross-validation accuracy**: 86-92% on unseen studies
- **Calibrated threshold**: 95% TPR at 10% FPR
- **Benchmark vs scrublet**: Superior on all test datasets
- **Weakest performance**: Somatotroph-lactotroph doublets

## Citation

If you use Epitome Tools in your research, please cite:

**For the models and package:**
```
Kövér et al. (2025). Consensus Pituitary Atlas, a scalable resource for annotation, 
novel marker discovery and age-sex analysis. [Manuscript in preparation]
```

**For the Epitome platform:**
```
Kövér, B., Kaufman-Cook, J., Sherwin, O., Vazquez Segoviano, M., Kemkem, Y., 
Lu, H.-C., & Andoniadou, C. (2025). Electronic Pituitary Omics (epitome) platform. 
Zenodo. https://doi.org/10.5281/zenodo.17154160
```

## Related Resources

- **Epitome Platform**: https://epitome.sites.er.kcl.ac.uk/ - Interactive web interface for exploring the Consensus Pituitary Atlas
- **Consensus Pituitary Atlas**: 1.1M+ cells, uniformly processed

## Support

For issues, questions, or feature requests:
- **GitHub Issues**: https://github.com/Andoniadou-Lab/epitome_tools/issues
- **Email**: epitome@kcl.ac.uk


## Acknowledgments
This work was supported by the Wellcome Trust Advanced Therapies for Regenerative Medicine PhD Programme (218461/Z/19/Z). Special thanks to the Andoniadou Lab at King's College London and all contributors to the Consensus Pituitary Atlas.

---

**Developer and Lead Curator**: Bence Kövér  
**Lab**: Andoniadou Lab, King's College London  
**Contact**: bence.kover@kcl.ac.uk
