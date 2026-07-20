```
   _____ _     _____ _____  ____  _           _ __     ___     
  / ____| |   |_   _|  __ \|  _ \(_)         | |\ \   / (_)    
 | |    | |__   | | | |__) | |_) |_ _ __   __| | \ \ / / _ ___ 
 | |    | '_ \  | | |  ___/|  _ <| | '_ \ / _` |  \ V / | / __|
 | |____| | | |_| |_| |    | |_) | | | | | (_| |   | |  | \__ \
  \_____|_| |_|_____|_|    |____/|_|_| |_|\__,_|   |_|  |_|___/
```

*ChIP-seq signal and promoter-associated peak visualization around gene TSSs*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19613844.svg)](https://doi.org/10.5281/zenodo.19613844)
[![R-CMD-check](https://github.com/ChanghaoKan/ChIPBindVis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ChanghaoKan/ChIPBindVis/actions/workflows/R-CMD-check.yaml)
![R >= 4.3](https://img.shields.io/badge/R->=4.3-blue?logo=r)
![Bioconductor](https://img.shields.io/badge/Bioconductor-powered-green)
![License: MIT](https://img.shields.io/badge/license-MIT-yellow)
![Version](https://img.shields.io/badge/version-0.2.1-orange)

---

## What it does

You have a **bigWig** signal file and a **bed** (or **broadPeak**) file from a ChIP-seq experiment. You want to inspect signal near gene X and ask whether its configured promoter window has an associated called peak, then compare TSS-centered signal across other genes with promoter-associated peaks.

ChIPBindVis produces **two complementary views in one function call**:

| Output | Description |
|---|---|
| **Selected-gene track** | Gviz signal, called peaks, TSS marker, and nearby gene models around the selected TSS. |
| **TSS enrichment heatmap** | ChIP-seq signal across gene-level TSSs that have promoter-associated peaks, with conditional query highlighting. |

Built with a consistent **Morandi color palette** — muted tones intended as a configurable starting point for manuscript figures.

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
devtools::install_local("ChIPBindVis-0.2.1.tar.gz")
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
  peaks_file  = "ENCFF049BWK.bed.gz",
  gene_symbol = "KIF18A",
  tf_name     = "E2F1",
  genome      = "hg38"        # explicit — do not omit
)
```

These two ENCODE files belong to the released
[E2F1 TF ChIP-seq experiment ENCSR717ZZW](https://www.encodeproject.org/experiments/ENCSR717ZZW/)
in HepG2: [ENCFF354YZN](https://www.encodeproject.org/files/ENCFF354YZN/)
is the GRCh38 signal p-value bigWig and
[ENCFF049BWK](https://www.encodeproject.org/files/ENCFF049BWK/) is the
IDR-thresholded narrowPeak file. The `tf_name` is therefore `"E2F1"`, not
`"E2F8"`.

### Human hg19

```r
chip_bindingVis(
  bigwig_file = "signal_hg19.bigWig",
  peaks_file  = "peaks_hg19.bed",
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
  peaks_file  = "E2f8_mm10.bed",
  gene_symbol = "Kif18a",     # mouse capitalisation
  tf_name     = "E2f8",
  genome      = "mm10"        # explicit
)
```

### Mouse mm39

```r
chip_bindingVis(
  bigwig_file = "E2f8_mm39.bigWig",
  peaks_file  = "E2f8_mm39.bed",
  gene_symbol = "Kif18a",
  tf_name     = "E2f8",
  genome      = "mm39"        # explicit
)
```

### Save to PDF

```r
chip_bindingVis(
  bigwig_file = "ENCFF354YZN.bigWig",
  peaks_file  = "ENCFF049BWK.bed.gz",
  gene_symbol = "KIF18A",
  tf_name     = "E2F1",
  genome      = "hg38",
  save_pdf    = TRUE           # output: E2F1_KIF18A_ChIPBindVis.pdf
)
```

---

## Usage

### Individual functions

```r
# Figure A only — Gviz track, human hg38
plot_chip_track(
  bigwig_file = "signal.bw",
  peaks_file  = "peaks.bed",
  gene_symbol = "KIF18A",
  tf_name     = "E2F8",
  genome      = "hg38",
  extend      = 10000         # ±10 kb around TSS
)

# Figure A only — mouse mm10
plot_chip_track(
  bigwig_file = "signal.bw",
  peaks_file  = "peaks.bed",
  gene_symbol = "Kif18a",
  tf_name     = "E2f8",
  genome      = "mm10"
)

# Figure B only — heatmap, returns ranking info
res <- plot_chip_heatmap(
  bigwig_file = "signal.bw",
  peaks_file  = "peaks.bed",
  gene_symbol = "KIF18A",
  tf_name     = "E2F8",
  genome      = "hg38"
)
res$rank                          # rank among promoter-associated gene TSSs
res$n_targets                     # number of TSS rows in the heatmap
res$query_promoter_peak_detected  # TRUE only if the query has a promoter peak
res$query_status                  # explicit detected/ranked status

# Gene info lookup
get_gene_info("KIF18A",  genome = "hg38")  # human
get_gene_info("Kif18a",  genome = "mm10")  # mouse mm10
get_gene_info("Kif18a",  genome = "mm39")  # mouse mm39
```

### Custom palette

```r
chip_bindingVis(
  bigwig_file = "signal.bw",
  peaks_file  = "peaks.bed",
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
  peaks_file    = "peaks.bed",
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
| Called peaks | `.bed` or `.broadPeak` | MACS2, ENCODE |
| Cytoband *(optional)* | `cytoBand.txt.gz` | UCSC (auto-downloaded on first use) |

> **broadPeak** files (e.g. from histone ChIP-seq) are automatically detected
> from the filename — no extra argument needed.

## Interpretation

`plot_chip_heatmap()` uses `ChIPseeker::annotatePeak()` and retains only
annotations classified as promoter peaks within `tss_region`. The query gene
is highlighted and ranked only when it is in that promoter-associated set. If
it is absent, ChIPBindVis does **not** insert it into the heatmap; `rank` is
`NA`, `query_promoter_peak_detected` is `FALSE`, and `query_status` explains
the result.

Each row uses a gene-level TSS derived from the TxDb gene range. It is not a
claim about a canonical transcript or a specific isoform. A promoter-associated
ChIP-seq peak supports an association in the supplied experiment, but by itself
does not establish direct TF–DNA binding or transcriptional regulation.

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

For reproducibility, cite the archived version used. The v0.2.1 archive is
available at <https://doi.org/10.5281/zenodo.21461267>. The concept DOI
<https://doi.org/10.5281/zenodo.19613844> always resolves to the latest
archived version.

> Kan, C. (2026). ChIPBindVis: ChIP-seq Signal and Peak Visualization at Gene
> Promoters (Version 0.2.1). Zenodo.
> <https://doi.org/10.5281/zenodo.21461267>

```bibtex
@software{kan2026chipbindvis,
  author    = {Kan, Changhao},
  title     = {{ChIPBindVis: ChIP-seq Signal and Peak Visualization at Gene
               Promoters}},
  version   = {0.2.1},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.21461267},
  url       = {https://doi.org/10.5281/zenodo.21461267}
}
```

## License

[MIT](LICENSE.md)
