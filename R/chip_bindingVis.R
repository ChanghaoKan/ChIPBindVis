#' One-command ChIP-seq binding visualization (Figures A + B)
#'
#' High-level wrapper that calls \code{\link{plot_chip_track}} and
#' \code{\link{plot_chip_heatmap}} in sequence. Both figures are drawn on
#' the active graphics device, or saved to a multi-page PDF when
#' \code{save_pdf = TRUE}.
#'
#' @inheritParams plot_chip_track
#' @inheritParams plot_chip_heatmap
#' @param save_pdf Logical. If \code{TRUE}, save both figures to a PDF
#'   (default \code{FALSE}).
#' @param pdf_file Output PDF path. Defaults to
#'   \code{<tf_name>_<gene_symbol>_ChIPBindVis.pdf}.
#' @param pdf_width,pdf_height Track-page dimensions (inches). Default 8 × 4.
#' @param pdf_width_heatmap,pdf_height_heatmap Heatmap-page dimensions
#'   (inches). Default 6 × 8.
#'
#' @return Invisibly, a list with elements \code{track} and \code{heatmap}
#'   (return values of the two underlying functions).
#' @export
#'
#' @examples
#' \dontrun{
#' # Human (default hg38)
#' chip_bindingVis(
#'   bigwig_file = "E2F8.bigWig",
#'   peaks_file  = "E2F8.narrowPeak",
#'   gene_symbol = "KIF18A",
#'   tf_name     = "E2F8"
#' )
#'
#' # Mouse mm10, save to PDF
#' chip_bindingVis(
#'   bigwig_file = "E2f8_mm10.bigWig",
#'   peaks_file  = "E2f8_mm10.narrowPeak",
#'   gene_symbol = "Kif18a",
#'   tf_name     = "E2f8",
#'   genome      = "mm10",
#'   save_pdf    = TRUE
#' )
#' }
chip_bindingVis <- function(bigwig_file,
                            peaks_file,
                            gene_symbol,
                            tf_name,
                            genome             = "hg38",
                            extend             = 5000,
                            tss_window         = 2000,
                            tss_region         = c(-2000, 2000),
                            palette            = morandi_palette(),
                            cytoband_file      = NULL,
                            save_pdf           = FALSE,
                            pdf_file           = NULL,
                            pdf_width          = 8,
                            pdf_height         = 4,
                            pdf_width_heatmap  = 6,
                            pdf_height_heatmap = 8) {

  if (isTRUE(save_pdf)) {
    if (is.null(pdf_file))
      pdf_file <- paste0(tf_name, "_", gene_symbol, "_ChIPBindVis.pdf")
    grDevices::pdf(pdf_file,
                   width  = max(pdf_width, pdf_width_heatmap),
                   height = max(pdf_height, pdf_height_heatmap))
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  track_res <- plot_chip_track(
    bigwig_file   = bigwig_file,
    peaks_file    = peaks_file,
    gene_symbol   = gene_symbol,
    tf_name       = tf_name,
    genome        = genome,
    extend        = extend,
    palette       = palette,
    cytoband_file = cytoband_file,
    plot          = TRUE
  )

  heatmap_res <- plot_chip_heatmap(
    bigwig_file = bigwig_file,
    peaks_file  = peaks_file,
    gene_symbol = gene_symbol,
    tf_name     = tf_name,
    genome      = genome,
    tss_window  = tss_window,
    tss_region  = tss_region,
    palette     = palette,
    plot        = TRUE
  )

  if (isTRUE(save_pdf))
    message("Saved to: ", normalizePath(pdf_file, mustWork = FALSE))

  invisible(list(track = track_res, heatmap = heatmap_res))
}
