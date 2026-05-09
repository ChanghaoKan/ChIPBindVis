```
   _____ _     _____ _____  ____  _           _ __     ___     
  / ____| |   |_   _|  __ \|  _ \(_)         | |\ \   / (_)    
 | |    | |__   | | | |__) | |_) |_ _ __   __| | \ \ / / _ ___ 
 | |    | '_ \  | | |  ___/|  _ <| | '_ \ / _` |  \ V / | / __|
 | |____| | | |_| |_| |    | |_) | | | | | (_| |   | |  | \__ \
  \_____|_| |_|_____|_|    |____/|_|_| |_|\__,_|   |_|  |_|___/
```

*One-command publication-ready ChIP-seq binding visualization at target gene promoters*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19613845.svg)](https://doi.org/10.5281/zenodo.19613845)
![R >= 4.3](https://img.shields.io/badge/R->=4.3-blue?logo=r)
![Bioconductor](https://img.shields.io/badge/Bioconductor-powered-green)
![License: MIT](https://img.shields.io/badge/license-MIT-yellow)
![Version](https://img.shields.io/badge/version-0.1.0-orange)

---

## What it does

You have a **bigWig** signal file and a **narrowPeak** file from a ChIP-seq experiment. You want to know: *does my transcription factor bind the promoter of gene X, and how does that compare to all its other targets?*

ChIPBindVis answers this in **one function call** and **two figures**:

| | |
|---|---|
| **Figure A — Track Plot** <br> *Gviz genome browser view centered on TSS* | **Figure B — Enrichment Heatmap** <br> *EnrichedHeatmap across all TF target gene TSSs* |

Built with a consistent **Morandi color palette** out of the box — muted, print-friendly tones that look good in manuscripts without any tweaking.

## Installation

```r
# 1. Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "Gviz", "rtracklayer", "GenomicRanges", "IRanges", "GenomeInfoDb",
  "GenomicFeatures", "ChIPseeker", "EnrichedHeatmap", "ComplexHeatmap",
  "circlize", "BiocParallel", "AnnotationDbi",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"
))

# 2. Install ChIPBindVis from GitHub
# install.packages("devtools")
devtools::install_github("ChanghaoKan/ChIPBindVis")
```

Or install locally from a downloaded source:

```r
devtools::install_local("ChIPBindVis-0.1.0.tar.gz")
# or, if you have the unpacked source folder:
devtools::install("path/to/ChIPBindVis")
```

## Quick Start

```r
library(ChIPBindVis)

# One command -> two publication-ready figures
chip_bindingVis(
  bigwig_file = "ENCFF354YZN.bigWig",
  peaks_file  = "ENCFF049BWK.bed",
  gene_symbol = "KIF18A",
  tf_name     = "E2F1"
)

# Save directly to PDF
chip_bindingVis(
  bigwig_file = "ENCFF354YZN.bigWig",
  peaks_file  = "ENCFF049BWK.bed",
  gene_symbol = "KIF18A",
  tf_name     = "E2F1",
  save_pdf    = TRUE
)
```

## Usage

### Individual functions

```r
# Track plot only (with wider window)
plot_chip_track("signal.bw", "peaks.bed", "KIF18A", "E2F1", extend = 10000)

# Heatmap only (returns ranking info)
res <- plot_chip_heatmap("signal.bw", "peaks.bed", "KIF18A", "E2F1")
res$target_rank
```

### Custom palette

```r
chip_bindingVis(
  bigwig_file = "signal.bw",
  peaks_file  = "peaks.bed",
  gene_symbol = "MYC",
  tf_name     = "MAX",
  palette     = morandi_palette(
    signal      = "#2E4057",
    signal_fill = "#4A6FA5",
    peak        = "#D64933"
  )
)
```

## API Reference

| Function | Description |
|---|---|
| `chip_bindingVis()` | High-level wrapper — generates both figures, optional PDF export |
| `plot_chip_track()` | Gviz genome browser track at a single gene promoter |
| `plot_chip_heatmap()` | EnrichedHeatmap across all target gene TSSs |
| `get_gene_info()` | Resolve HGNC symbol → chr, start, end, strand, TSS |
| `morandi_palette()` | Customizable Morandi color palette with named overrides |

## Input files

| File | Format | Source |
|---|---|---|
| Signal track | `.bigWig` | ENCODE, GEO, or your own pipeline |
| Called peaks | `.narrowPeak` / `.bed` | MACS2, ENCODE |
| Cytoband (optional) | `cytoBand.txt.gz` | UCSC (auto-downloaded on first use) |

The hg38 cytoBand file is fetched once from UCSC and cached at
`tools::R_user_dir("ChIPBindVis", "cache")`. If you are offline, you can
supply your own copy via `cytoband_file = "/path/to/cytoBand.txt.gz"`.

## Genome

Currently hard-coded to **human hg38**, using
`TxDb.Hsapiens.UCSC.hg38.knownGene` and `org.Hs.eg.db`. Other genomes are
not yet supported.

## Acknowledgments

This R package was developed with assistance from **Claude** (Anthropic).

## Citation

If you use ChIPBindVis in your research, please cite:

> Kan, C. (2025). ChIPBindVis: One-command publication-ready ChIP-seq binding visualization at target gene promoters (v0.1.0). Zenodo. <https://doi.org/10.5281/zenodo.19613845>

```bibtex
@software{kan2025chipbindvis,
  author    = {Kan, Changhao},
  title     = {{ChIPBindVis: One-command publication-ready ChIP-seq 
               binding visualization at target gene promoters}},
  year      = {2025},
  version   = {v0.1.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.19613845},
  url       = {https://doi.org/10.5281/zenodo.19613845}
}
```

## License

[MIT](LICENSE)
