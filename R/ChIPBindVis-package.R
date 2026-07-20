#' ChIPBindVis: ChIP-seq Signal and Peak Visualization at Gene Promoters
#'
#' Visualizes ChIP-seq signal and called peaks around a selected gene and
#' across gene-level TSSs with promoter-associated peaks. The query gene is
#' highlighted only when a promoter-associated peak is detected within the
#' configured window. Supports human (hg38, hg19) and mouse (mm10, mm39)
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
