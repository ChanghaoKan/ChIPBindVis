# ChIPBindVis 0.2.1

## Correctness and reproducibility updates

* Corrected the bundled ENCODE example metadata: `ENCFF354YZN` and
  `ENCFF049BWK` are E2F1 files from experiment `ENCSR717ZZW`, not E2F8.
* Heatmap rows are now restricted to genes with annotations classified as
  promoter peaks within `tss_region`; distal nearest-gene assignments are not
  described as promoter-associated targets.
* A query gene without a promoter-associated peak is no longer inserted into
  the heatmap. `rank` is `NA`, and `query_promoter_peak_detected` plus
  `query_status` make the outcome explicit.
* Heatmap TSSs are derived from TxDb gene ranges and documented as gene-level
  proxies rather than canonical-transcript TSSs.

# ChIPBindVis 0.2.0

## New features

* **Multi-genome support** — added `hg19`, `mm10`, and `mm39` alongside the
  original `hg38`. Annotation packages are loaded on demand, so users only
  install what they need.
* **broadPeak support** — histone ChIP-seq peak files are auto-detected
  from the filename; no extra argument required.
* **Automatic strand correction** — TSS coordinates and `normalizeToMatrix`
  windows are now reversed for negative-strand genes, so the
  `-Xkb / TSS / +Xkb` axis is biologically accurate for every gene.
* **Cached cytoBand download** — the UCSC cytoBand file is downloaded once
  and cached at `tools::R_user_dir("ChIPBindVis", "cache")`. Offline use
  is supported via the `cytoband_file =` argument.
* **Local install** — `devtools::install_local("ChIPBindVis-0.2.0.tar.gz")`
  now works out of the box.

## Breaking changes

* All plotting functions now take a `genome =` argument. Default is
  `"hg38"`; passing it explicitly is strongly recommended (especially for
  mouse data — silent defaults will produce wrong annotations).

## API changes

* `plot_chip_heatmap()` return value includes `n_targets` (total
  number of target TSSs in the heatmap) alongside `rank`.
* Mouse gene symbols must follow title case (e.g. `"Kif18a"`, not
  `"KIF18A"`).

## Bug fixes

* Corrected the install path in the README: `ChanghaoKan/ChIPBindVis`
  (previously pointed to a non-existent `Changhao-Kan/ChIPBindVis`).

# ChIPBindVis 0.1.0

* First public release.
* `chip_bindingVis()` wrapper produces a Gviz track plot and an
  EnrichedHeatmap from a bigWig + narrowPeak pair, centered on a target
  gene promoter.
* Human (hg38) only.
* Built-in Morandi color palette, customizable via `morandi_palette()`.
* Optional one-line PDF export.
