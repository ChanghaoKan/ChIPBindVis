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
![Version](https://img.shields.io/badge/version-0.2.0-orange)

---

## What it does

You have a **bigWig** signal file and a **narrowPeak** (or **broadPeak**) file from a ChIP-seq experiment. You want to know: *does my transcription factor bind the promoter of gene X, and how does that compare to all its other targets?*

ChIPBindVis answers this in **one function call** and **two figures**:

| Figure | Preview |
|---|---|
| **Track Plot & Enrichment Heatmap** | <img width="1007" height="503" alt="image" src="https://github.com/user-attachments/assets/46f276dc-1351-4d3c-8868-ce3a703bb5d6" /> |
| Description | *Gviz genome browser tracks centered on transcription start sites (TSSs) together with EnrichedHeatmap signals across TF target gene promoters.* |

Built with a consistent **Morandi color palette** — muted, publication-ready tones designed for clean manuscript integration without additional styling.

---

## Supported genomes

| `genome =` | Species | Annotation packages required |
|---|---|---|
| `"hg38"` *(default)* | Human GRCh38 | `TxDb.Hsapiens.UCSC.hg38.knownGene`, `org.Hs.eg.db` |
| `"hg19"` | Human GRCh37 | `TxDb.Hsapiens.UCSC.hg19.knownGene`, `org.Hs.eg.db` |
| `"mm10"` | Mouse GRCm38 | `TxDb.Mmusculus.UCSC.mm10.knownGene`, `org.Mm.eg.db` |
| `"mm39"` | Mouse GRCm39 | `TxDb.Mmusculus.UCSC.mm39.knownGene`, `org.Mm.eg.db` |

> **Always pass `genome` explicitly.** The default is `"hg38"`, but relying
> on it silently for mouse data will cause errors or wrong annotations.
> Annotation packages are loaded on demand — only install the one(s) you use.

---

## Installation

```r
# 1. Core Bioconductor dependencies (required for all genomes)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "Gviz", "rtracklayer", "GenomicRanges", "IRanges", "GenomeInfoDb",
  "GenomicFeatures", "ChIPseeker", "EnrichedHeatmap", "ComplexHeatmap",
  "circlize", "AnnotationDbi"
))

# 2. Genome-specific annotation packages — install only what you need
# Human hg38
BiocManager::install(c("TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))

# Human hg19
BiocManager::install(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"))

# Mouse mm10
BiocManager::install(c("TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"))

# Mouse mm39
BiocManager::install(c("TxDb.Mmusculus.UCSC.mm39.knownGene", "org.Mm.eg.db"))

# 3. Install ChIPBindVis from GitHub
# install.packages("devtools")
devtools::install_github("ChanghaoKan/ChIPBindVis")
```

Or install locally:

```r
devtools::install_local("ChIPBindVis-0.2.0.tar.gz")
# or from an unpacked source directory:
devtools::install("path/to/ChIPBindVis")
```

---

## Quick start

### Human hg38

```r
library(ChIPBindVis)

chip_bindingVis(
  bigwig_file = "ENCFF354YZN.bigWig",
  peaks_file  = "ENCFF049BWK.narrowPeak",
  gene_symbol = "KIF18A",
  tf_name     = "E2F8",
  genome      = "hg38"        # explicit — do not omit
)
```

### Human hg19

```r
chip_bindingVis(
  bigwig_file = "signal_hg19.bigWig",
  peaks_file  = "peaks_hg19.narrowPeak",
  gene_symbol = "KIF18A",
  tf_name     = "E2F8",
  genome      = "hg19"        # explicit
)
```

### Mouse mm10

```r
# Mouse gene symbols use title-case (e.g. "Kif18a", not "KIF18A")
chip_bindingVis(
  bigwig_file = "E2f8_mm10.bigWig",
  peaks_file  = "E2f8_mm10.narrowPeak",
  gene_symbol = "Kif18a",     # mouse capitalisation
  tf_name     = "E2f8",
  genome      = "mm10"        # explicit
)
```

### Mouse mm39

```r
chip_bindingVis(
  bigwig_file = "E2f8_mm39.bigWig",
  peaks_file  = "E2f8_mm39.narrowPeak",
  gene_symbol = "Kif18a",
  tf_name     = "E2f8",
  genome      = "mm39"        # explicit
)
```

### Save to PDF

```r
chip_bindingVis(
  bigwig_file = "ENCFF354YZN.bigWig",
  peaks_file  = "ENCFF049BWK.narrowPeak",
  gene_symbol = "KIF18A",
  tf_name     = "E2F8",
  genome      = "hg38",
  save_pdf    = TRUE           # output: E2F8_KIF18A_ChIPBindVis.pdf
)
```

---

## Usage

### Individual functions

```r
# Figure A only — Gviz track, human hg38
plot_chip_track(
  bigwig_file = "signal.bw",
  peaks_file  = "peaks.narrowPeak",
  gene_symbol = "KIF18A",
  tf_name     = "E2F8",
  genome      = "hg38",
  extend      = 10000         # ±10 kb around TSS
)

# Figure A only — mouse mm10
plot_chip_track(
  bigwig_file = "signal.bw",
  peaks_file  = "peaks.narrowPeak",
  gene_symbol = "Kif18a",
  tf_name     = "E2f8",
  genome      = "mm10"
)

# Figure B only — heatmap, returns ranking info
res <- plot_chip_heatmap(
  bigwig_file = "signal.bw",
  peaks_file  = "peaks.narrowPeak",
  gene_symbol = "KIF18A",
  tf_name     = "E2F8",
  genome      = "hg38"
)
res$target_rank   # integer rank among all TF targets
res$n_targets     # total number of target TSSs in the heatmap

# Gene info lookup
get_gene_info("KIF18A",  genome = "hg38")  # human
get_gene_info("Kif18a",  genome = "mm10")  # mouse mm10
get_gene_info("Kif18a",  genome = "mm39")  # mouse mm39
```

### Custom palette

```r
chip_bindingVis(
  bigwig_file = "signal.bw",
  peaks_file  = "peaks.narrowPeak",
  gene_symbol = "MYC",
  tf_name     = "MAX",
  genome      = "hg38",
  palette     = morandi_palette(
    signal      = "#2E4057",
    signal_fill = "#4A6FA5",
    peak        = "#D64933"
  )
)
```

### Offline / air-gapped use

The cytoBand file for the selected genome is downloaded from UCSC on first
use and cached at `tools::R_user_dir("ChIPBindVis", "cache")`. To avoid the
download, supply a pre-downloaded copy:

```r
chip_bindingVis(
  bigwig_file   = "signal.bw",
  peaks_file    = "peaks.narrowPeak",
  gene_symbol   = "KIF18A",
  tf_name       = "E2F8",
  genome        = "hg38",
  cytoband_file = "/path/to/cytoBand_hg38.txt.gz"
)
```

---

## API reference

| Function | Description |
|---|---|
| `chip_bindingVis()` | High-level wrapper — both figures, optional PDF export |
| `plot_chip_track()` | Figure A: Gviz genome browser track |
| `plot_chip_heatmap()` | Figure B: TSS enrichment heatmap |
| `get_gene_info()` | Gene symbol → chr / start / end / strand / TSS |
| `morandi_palette()` | Customizable Morandi color palette |

---

## Input files

| File | Format | Source |
|---|---|---|
| Signal track | `.bigWig` / `.bw` | ENCODE, GEO, or your own pipeline |
| Called peaks | `.narrowPeak` or `.broadPeak` | MACS2, ENCODE |
| Cytoband *(optional)* | `cytoBand.txt.gz` | UCSC (auto-downloaded on first use) |

> **broadPeak** files (e.g. from histone ChIP-seq) are automatically detected
> from the filename — no extra argument needed.

---

## Strand handling

TSS coordinates are strand-corrected automatically: positive-strand (`+`)
genes use `start`, negative-strand (`-`) genes use `end`. In Figure B,
`normalizeToMatrix` reverses the signal window for negative-strand targets
so that the left side of every heatmap row consistently represents upstream
sequence — the `-Xkb / TSS / +Xkb` axis labels are biologically accurate
for all genes regardless of strand.

---

## Common mistakes

| Mistake | Result | Fix |
|---|---|---|
| Using mouse data without `genome = "mm10"` | `Gene not found` error | Always pass `genome` explicitly |
| Human-style symbol for mouse (`"KIF18A"`) | `Gene not found` error | Use mouse capitalisation: `"Kif18a"` |
| Missing annotation package for the genome | `Required package not installed` error | Run the genome-specific `BiocManager::install()` above |

---

## Acknowledgments

This R package was developed with assistance from **Claude** (Anthropic).

---

## Citation

If you use ChIPBindVis in your research, please cite:

> Kan, C. (2025). ChIPBindVis: One-command publication-ready ChIP-seq
> binding visualization at target gene promoters (v0.2.0). Zenodo.
> <https://doi.org/10.5281/zenodo.19613845>

```bibtex
@software{kan2025chipbindvis,
  author    = {Kan, Changhao},
  title     = {{ChIPBindVis: One-command publication-ready ChIP-seq
               binding visualization at target gene promoters}},
  year      = {2025},
  version   = {v0.2.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.19613845},
  url       = {https://doi.org/10.5281/zenodo.19613845}
}
```

## License

[MIT](LICENSE)
