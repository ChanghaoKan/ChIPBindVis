#' ChIPBindVis: ChIP-seq Binding Visualization at Target Gene Promoters
#'
#' One-command, publication-ready ChIP-seq binding visualization at target gene
#' promoters. Given a bigWig signal track and a narrowPeak file, the package
#' produces (A) a Gviz genome-browser view centered on a target gene's TSS, and
#' (B) an EnrichedHeatmap across all peak-annotated target gene TSSs, with the
#' query gene highlighted.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{chip_bindingVis}} — high-level wrapper, both figures
#'   \item \code{\link{plot_chip_track}} — Figure A only (single-gene track)
#'   \item \code{\link{plot_chip_heatmap}} — Figure B only (TSS heatmap)
#'   \item \code{\link{get_gene_info}} — symbol -> chr / start / end / strand / TSS
#'   \item \code{\link{morandi_palette}} — built-in muted color palette
#' }
#'
#' @section Genome:
#' Hard-coded to human hg38 via
#' \pkg{TxDb.Hsapiens.UCSC.hg38.knownGene} and \pkg{org.Hs.eg.db}.
#'
#' @keywords internal
"_PACKAGE"

## Silence R CMD check NOTEs about pipe-style symbol references
utils::globalVariables(c("score", "SYMBOL", "tx_name"))
