# GDSCtools

Utility functions for downstream single-cell RNA-seq analysis, providing a
pipeline from differential expression through GO enrichment to network
visualization.

## Installation

```r
# Install from local source
install.packages("path/to/GDSC_Utils", repos = NULL, type = "source")

# Or with devtools
devtools::install_local("path/to/GDSC_Utils")
```

## Dependencies

**CRAN:** igraph, dplyr, ggplot2, ggwordcloud, patchwork, scales, stringr

**Bioconductor:** Seurat, clusterProfiler, org.Mm.eg.db (or org.Hs.eg.db)

Install Bioconductor dependencies first if needed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db"))
```

## Functions

### Differential Expression

| Function | Description |
|---|---|
| `run_cluster_de()` | Cluster-specific DE via `Seurat::FindMarkers` with fold-change and p-value filtering |
| `write_de_results()` | Export DE or GO result lists to per-cluster CSV files |

### GO Enrichment

| Function | Description |
|---|---|
| `run_go_ora()` | GO over-representation analysis per cluster via `clusterProfiler::enrichGO` with QC summary |

### Network Analysis

| Function | Description |
|---|---|
| `run_comparison_networks()` | Top-level driver: build overlap networks for all clusters in a comparison |
| `build_network_inputs()` | Extract significant GO terms and pair with fold-change data |
| `generate_overlap_network()` | Build Jaccard-overlap network, Leiden clustering, plots and CSV export |
| `write_cluster_results()` | Extract all GO terms from a specific Leiden cluster |

### Visualization

| Function | Description |
|---|---|
| `write_cluster_wordclouds_patchwork()` | Multi-panel word clouds summarizing GO term themes per Leiden cluster |

## Typical Workflow

```r
library(GDSCUtils)
library(org.Mm.eg.db)

# 1. Differential expression
de_results <- run_cluster_de(
  obj        = seurat_obj,
  clusters   = c("Fibroblast 1", "Fibroblast 2"),
  ident1     = "Car",
  ident2     = "Bleo"
)

write_de_results(de_results, outSubDir = "results/DE/", type = "DE")

# 2. GO enrichment
universe <- rownames(seurat_obj)

go_results <- run_go_ora(
  de_list   = de_results,
  org_db    = org.Mm.eg.db,
  universe  = universe,
  ont       = "BP",
  outputDir = "results/GO/"
)

write_de_results(go_results, outSubDir = "results/GO/", type = "GO",
                 split_direction = FALSE)

# 3. Network analysis
net_results <- run_comparison_networks(
  go_list     = go_results,
  de_list     = de_results,
  comp_label  = "Car_vs_Bleo",
  base_outdir = "results/"
)

# 4. Word clouds
write_cluster_wordclouds_patchwork(
  net_list   = net_results,
  outputDir  = "results/",
  comparison = "Car_vs_Bleo"
)

# 5. Inspect a specific cluster
ecm_terms <- write_cluster_results(
  net_result = net_results[["Fibroblast_1"]],
  go_term    = "extracellular matrix organization",
  outdir     = "results/Car_vs_Bleo/networks/"
)
```

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `overlap_threshold` | 0.1 | Minimum Jaccard overlap to draw a network edge |
| `leiden_resolution` | 0.03 | Leiden clustering resolution (lower = fewer clusters) |
| `merge_similarity` | 0.6 | Word-level Jaccard threshold for merging redundant GO labels |
| `padj_thresh` | 0.05 | Adjusted p-value cutoff for DE filtering |
| `logfc_thresh` | 0.5 | Minimum absolute log2FC for DE filtering |
| `min_terms` | 5 | Minimum phrases to render a word cloud panel |

## Attribution

Network analysis approach adapted from `pathway-analysis.R` by Owen Wilkins.

## License

MIT
