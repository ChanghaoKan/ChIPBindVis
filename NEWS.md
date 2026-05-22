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

* `plot_chip_heatmap()` return value now includes `n_targets` (total
  number of target TSSs in the heatmap) alongside `target_rank`.
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
