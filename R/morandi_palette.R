#' Morandi color palette for ChIPBindVis
#'
#' Returns a named list of muted, print-friendly Morandi tones used as the
#' default palette across all ChIPBindVis plotting functions. Pass any subset
#' of arguments to override individual colors.
#'
#' @param signal Color of the signal track outline (default \code{"#5B7C99"}).
#' @param signal_fill Fill color of the signal track polygon
#'   (default \code{"#7D9BB3"}).
#' @param peak Color of peak annotations and target highlight
#'   (default \code{"#B87B7B"}).
#' @param gene Color of gene-region track (default \code{"#6B8E7B"}).
#' @param tss Color of TSS marker (default \code{"#C4A77D"}).
#' @param heatmap_low,heatmap_mid,heatmap_high Heatmap color ramp endpoints
#'   (defaults \code{"#F5F2EB"}, \code{"#B8C5D0"}, \code{"#5B7C99"}).
#' @param axis Color of axis labels and titles (default \code{"#4A4A4A"}).
#' @param bg Background color (default \code{"white"}).
#'
#' @return A named list of hex color strings.
#' @export
#'
#' @examples
#' pal <- morandi_palette()
#' pal$signal
#'
#' # Override just the peak color
#' pal2 <- morandi_palette(peak = "#D64933")
morandi_palette <- function(signal       = "#5B7C99",
                            signal_fill  = "#7D9BB3",
                            peak         = "#B87B7B",
                            gene         = "#6B8E7B",
                            tss          = "#C4A77D",
                            heatmap_low  = "#F5F2EB",
                            heatmap_mid  = "#B8C5D0",
                            heatmap_high = "#5B7C99",
                            axis         = "#4A4A4A",
                            bg           = "white") {
  list(
    signal       = signal,
    signal_fill  = signal_fill,
    peak         = peak,
    gene         = gene,
    tss          = tss,
    heatmap_low  = heatmap_low,
    heatmap_mid  = heatmap_mid,
    heatmap_high = heatmap_high,
    axis         = axis,
    bg           = bg
  )
}
