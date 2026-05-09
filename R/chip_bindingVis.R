#' One-command ChIP-seq binding visualization (Figures A + B)
#'
#' High-level wrapper that calls \code{\link{plot_chip_track}} and
#' \code{\link{plot_chip_heatmap}} in sequence. Both figures are drawn on
#' the active graphics device, or saved to a single multi-page PDF when
#' \code{save_pdf = TRUE}.
#'
#' @inheritParams plot_chip_track
#' @inheritParams plot_chip_heatmap
#' @param save_pdf Logical. If \code{TRUE}, save the two figures to a PDF
#'   instead of drawing on screen (default \code{FALSE}).
#' @param pdf_file Output PDF path. Used only when \code{save_pdf = TRUE}.
#'   Default: \code{paste0(tf_name, "_", gene_symbol, "_ChIPBindVis.pdf")}.
#' @param pdf_width,pdf_height PDF dimensions in inches.
#'   Defaults: \code{8 x 4} for the track, \code{6 x 8} for the heatmap.
#' @param pdf_width_heatmap,pdf_height_heatmap Override the heatmap-page size.
#'
#' @return Invisibly, a list with elements \code{track} and \code{heatmap}
#'   (the return values of the two underlying functions).
#' @export
#'
#' @examples
#' \dontrun{
#' chip_bindingVis(
#'   bigwig_file = "ENCFF354YZN.bigWig",
#'   peaks_file  = "ENCFF049BWK.bed",
#'   gene_symbol = "KIF18A",
#'   tf_name     = "E2F1"
#' )
#'
#' chip_bindingVis(
#'   bigwig_file = "ENCFF354YZN.bigWig",
#'   peaks_file  = "ENCFF049BWK.bed",
#'   gene_symbol = "KIF18A",
#'   tf_name     = "E2F1",
#'   save_pdf    = TRUE
#' )
#' }
chip_bindingVis <- function(bigwig_file,
                            peaks_file,
                            gene_symbol,
                            tf_name,
                            extend            = 5000,
                            tss_window        = 2000,
                            tss_region        = c(-2000, 2000),
                            palette           = morandi_palette(),
                            cytoband_file     = NULL,
                            save_pdf          = FALSE,
                            pdf_file          = NULL,
                            pdf_width         = 8,
                            pdf_height        = 4,
                            pdf_width_heatmap = 6,
                            pdf_height_heatmap = 8) {

  if (isTRUE(save_pdf)) {
    if (is.null(pdf_file)) {
      pdf_file <- paste0(tf_name, "_", gene_symbol, "_ChIPBindVis.pdf")
    }

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
    tss_window  = tss_window,
    tss_region  = tss_region,
    palette     = palette,
    plot        = TRUE
  )

  if (isTRUE(save_pdf)) {
    message("Saved figures to: ", normalizePath(pdf_file, mustWork = FALSE))
  }

  invisible(list(track = track_res, heatmap = heatmap_res))
}
