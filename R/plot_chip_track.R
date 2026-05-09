#' Plot a ChIP-seq Gviz track at a target gene's promoter (Figure A)
#'
#' Builds a multi-track Gviz figure (ideogram, axis, signal, peaks, TSS,
#' gene model) centered on the TSS of \code{gene_symbol}, with signal taken
#' from \code{bigwig_file} and peak rectangles from \code{peaks_file}.
#'
#' @param bigwig_file Path to a bigWig signal file (\code{.bw} / \code{.bigWig}).
#' @param peaks_file Path to a narrowPeak / BED peak file.
#' @param gene_symbol HGNC gene symbol of the target gene (e.g. \code{"KIF18A"}).
#' @param tf_name Display name of the transcription factor (used in the
#'   track title and main title).
#' @param extend Half-window in bp around the TSS (default \code{5000}, i.e.
#'   TSS \eqn{\pm} 5 kb).
#' @param palette A list of colors, typically from \code{\link{morandi_palette}}.
#' @param cytoband_file Optional path to a user-supplied hg38
#'   \code{cytoBand.txt(.gz)}. If \code{NULL}, the file is auto-downloaded
#'   from UCSC and cached on first use.
#' @param plot Logical. If \code{TRUE} (default) the figure is drawn on the
#'   active graphics device.
#'
#' @return Invisibly, a list with \code{tracks} (the list of Gviz track
#'   objects), \code{sizes} (relative track sizes), \code{from}, \code{to},
#'   and \code{target} (the resolved gene info).
#' @export
plot_chip_track <- function(bigwig_file,
                            peaks_file,
                            gene_symbol,
                            tf_name,
                            extend        = 5000,
                            palette       = morandi_palette(),
                            cytoband_file = NULL,
                            plot          = TRUE) {

  TARGET <- get_gene_info(gene_symbol)

  message("Transcription factor: ", tf_name)
  message("Target gene: ", TARGET$symbol, " ",
          TARGET$chr, ":",
          format(TARGET$start, big.mark = ","), "-",
          format(TARGET$end,   big.mark = ","),
          " (", TARGET$strand, ")")
  message("TSS: ", format(TARGET$tss, big.mark = ","))

  ## Peaks ------------------------------------------------------------------
  peaks <- .read_peaks(peaks_file)

  ## Y-axis range from the signal in the view window -----------------------
  view_start <- TARGET$tss - extend
  view_end   <- TARGET$tss + extend

  bw_region <- rtracklayer::import(
    bigwig_file,
    format    = "BigWig",
    selection = GenomicRanges::GRanges(
      TARGET$chr,
      IRanges::IRanges(view_start, view_end)
    )
  )
  y_max <- ceiling(max(bw_region$score, na.rm = TRUE) * 1.1)
  if (is.na(y_max) || y_max <= 0) y_max <- 30

  ## Cytoband + ideogram ---------------------------------------------------
  cytoband <- .get_cytoband_hg38(cytoband_file)
  colnames(cytoband) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

  itrack <- Gviz::IdeogramTrack(
    genome     = "hg38",
    chromosome = TARGET$chr,
    bands      = cytoband,
    fontcolor  = palette$axis,
    fontsize   = 9
  )

  gtrack <- Gviz::GenomeAxisTrack(
    col       = palette$axis,
    fontcolor = palette$axis,
    fontsize  = 10
  )

  ## Signal track -----------------------------------------------------------
  signal_track <- Gviz::DataTrack(
    range            = bigwig_file,
    genome           = "hg38",
    chromosome       = TARGET$chr,
    name             = tf_name,
    type             = "polygon",
    fill.mountain    = palette$signal_fill,
    col.mountain     = palette$signal,
    ylim             = c(0, y_max),
    col.axis         = palette$axis,
    fontcolor.title  = palette$axis,
    background.title = palette$bg,
    cex.title        = 1.1,
    rotation.title   = 0
  )

  ## Peaks in view ----------------------------------------------------------
  peak_region <- GenomicRanges::GRanges(
    TARGET$chr, IRanges::IRanges(view_start, view_end)
  )
  peaks_in_view <- peaks[IRanges::overlapsAny(peaks, peak_region)]

  peak_track <- NULL
  if (length(peaks_in_view) > 0L) {
    peak_track <- Gviz::AnnotationTrack(
      peaks_in_view,
      name             = "Peaks",
      fill             = palette$peak,
      col              = palette$peak,
      stacking         = "dense",
      fontcolor.title  = palette$axis,
      background.title = palette$bg,
      rotation.title   = 0
    )
  }

  ## TSS marker -------------------------------------------------------------
  tss_track <- Gviz::AnnotationTrack(
    GenomicRanges::GRanges(
      TARGET$chr,
      IRanges::IRanges(TARGET$tss - 50, TARGET$tss + 50)
    ),
    name             = "TSS",
    fill             = palette$tss,
    col              = palette$tss,
    fontcolor.title  = palette$axis,
    background.title = palette$bg,
    rotation.title   = 0
  )

  ## Gene model -------------------------------------------------------------
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

  gene_track <- Gviz::GeneRegionTrack(
    txdb,
    genome               = "hg38",
    chromosome           = TARGET$chr,
    name                 = "Genes",
    transcriptAnnotation = "symbol",
    collapseTranscripts  = "meta",
    fill                 = palette$gene,
    col                  = palette$gene,
    fontcolor.group      = palette$axis,
    fontsize.group       = 11,
    fontface.group       = "bold.italic",
    fontcolor.title      = palette$axis,
    background.title     = palette$bg,
    rotation.title       = 0
  )

  ## Assemble ---------------------------------------------------------------
  if (!is.null(peak_track)) {
    track_list  <- list(itrack, gtrack, signal_track, peak_track,
                        tss_track, gene_track)
    track_sizes <- c(0.5, 1, 4, 0.5, 0.3, 1.5)
  } else {
    track_list  <- list(itrack, gtrack, signal_track,
                        tss_track, gene_track)
    track_sizes <- c(0.5, 1, 4, 0.3, 1.5)
  }

  if (isTRUE(plot)) {
    Gviz::plotTracks(
      track_list,
      from              = view_start,
      to                = view_end,
      sizes             = track_sizes,
      background.panel  = palette$bg,
      background.title  = palette$bg,
      col.border.title  = NA,
      main              = paste0(tf_name, " ChIP-seq at ",
                                 TARGET$symbol, " Promoter"),
      cex.main          = 1.4,
      fontface.main     = 2,
      col.main          = palette$axis
    )
  }

  invisible(list(
    tracks = track_list,
    sizes  = track_sizes,
    from   = view_start,
    to     = view_end,
    target = TARGET
  ))
}
