#' ChIPBindVis: ChIP-seq Binding Visualization at Target Gene Promoters
#'
#' One-command, publication-ready ChIP-seq binding visualization at target
#' gene promoters. Supports human (hg38, hg19) and mouse (mm10, mm39)
#' assemblies.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{chip_bindingVis}} — high-level wrapper, both figures
#'   \item \code{\link{plot_chip_track}} — Figure A: single-gene Gviz track
#'   \item \code{\link{plot_chip_heatmap}} — Figure B: TSS enrichment heatmap
#'   \item \code{\link{get_gene_info}} — gene symbol → chr/start/end/strand/TSS
#'   \item \code{\link{morandi_palette}} — built-in muted color palette
#' }
#'
#' @section Supported genomes:
#' \describe{
#'   \item{hg38}{Human GRCh38 (default)}
#'   \item{hg19}{Human GRCh37}
#'   \item{mm10}{Mouse GRCm38}
#'   \item{mm39}{Mouse GRCm39}
#' }
#' Required annotation packages are loaded on demand and only need to be
#' installed for the genome(s) actually used.
#'
#' @keywords internal
"_PACKAGE"

## Silence R CMD check NOTEs
utils::globalVariables(c("score", "SYMBOL", "tx_name"))
