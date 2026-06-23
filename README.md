# End-to-End Single-Cell RNA-Seq Preprocessing & Clustering Pipeline

An operational, reproducible computational pipeline built to preprocess, clean, and cluster a raw high-dimensional single-cell dataset (10k Human Peripheral Blood Mononuclear Cells) using standard quantitative data-science strategies.

---

## 📂 Project Architecture
```text
pbmc10k-pipeline/
├── data/                       # Initial genomics matrices (git-ignored)
├── src/
│   ├── 01_download_qc.sh       # Bash automation stream 
│   └── 02_analysis.R           # Seurat normalization & dimensionality reduction
└── results/
    └── figures/
        └── umap_clusters.png   # Final cluster visualization# pbmc10k-pipeline
