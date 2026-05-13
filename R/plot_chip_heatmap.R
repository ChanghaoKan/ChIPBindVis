#' Plot a ChIP-seq Gviz track at a target gene's promoter (Figure A)
#'
#' Builds a multi-track Gviz figure (ideogram, axis, signal, peaks, TSS,
#' gene model) centered on the TSS of \code{gene_symbol}. Supports human
#' (hg38 / hg19) and mouse (mm10 / mm39) assemblies.
#'
#' @param bigwig_file Path to a bigWig signal file (\code{.bw} /
#'   \code{.bigWig}).
#' @param peaks_file Path to a narrowPeak or broadPeak file.
#' @param gene_symbol Gene symbol of the target gene (e.g. \code{"KIF18A"}
#'   for human or \code{"Kif18a"} for mouse).
#' @param tf_name Display name of the transcription factor (used in track
#'   and main titles).
#' @param genome Reference genome assembly: \code{"hg38"} (default),
#'   \code{"hg19"}, \code{"mm10"}, or \code{"mm39"}.
#' @param extend Half-window in bp around the TSS (default \code{5000}).
#' @param palette A named list of colors, typically from
#'   \code{\link{morandi_palette}}.
#' @param cytoband_file Optional path to a pre-downloaded
#'   \code{cytoBand.txt(.gz)} for the chosen genome. If \code{NULL}, the file
#'   is auto-downloaded from UCSC and cached on first use.
#' @param plot Logical. If \code{TRUE} (default) the figure is drawn on
#'   the active graphics device.
#'
#' @return Invisibly, a list with \code{tracks}, \code{sizes}, \code{from},
#'   \code{to}, and \code{target} (the resolved gene info).
#' @export
#'
#' @examples
#' \dontrun{
#' # Human (hg38, default)
#' plot_chip_track("E2F8.bigWig", "E2F8.narrowPeak", "KIF18A", "E2F8")
#'
#' # Mouse (mm10)
#' plot_chip_track("E2f8.bigWig", "E2f8.narrowPeak", "Kif18a", "E2f8",
#'                 genome = "mm10")
#' }
plot_chip_track <- function(bigwig_file,
                            peaks_file,
                            gene_symbol,
                            tf_name,
                            genome        = "hg38",
                            extend        = 5000,
                            palette       = morandi_palette(),
                            cytoband_file = NULL,
                            plot          = TRUE) {

  TARGET <- get_gene_info(gene_symbol, genome = genome)
  res    <- .genome_resources(genome)
  txdb   <- .load_txdb(res$txdb)

  message("Genome: ", genome)
  message("Transcription factor: ", tf_name)
  message("Target gene: ", TARGET$symbol, "  ", TARGET$chr, ":",
          format(TARGET$start, big.mark = ","), "-",
          format(TARGET$end,   big.mark = ","),
          "  (", TARGET$strand, "-strand)")
  message("TSS: ", format(TARGET$tss, big.mark = ","))

  peaks      <- .read_peaks(peaks_file)
  view_start <- TARGET$tss - extend
  view_end   <- TARGET$tss + extend

  # Read only the view window from bigWig to determine y-axis ceiling
  bw_region <- rtracklayer::import(
    bigwig_file, format = "BigWig",
    selection = GenomicRanges::GRanges(
      TARGET$chr, IRanges::IRanges(view_start, view_end))
  )
  y_max <- ceiling(max(bw_region$score, na.rm = TRUE) * 1.1)
  if (is.na(y_max) || y_max <= 0) y_max <- 30

  # ── IdeogramTrack ───────────────────────────────────────────
  # .get_cytoband filters to TARGET$chr + resets rownames, which prevents
  # Gviz's "breaks not unique" error on narrow windows.
  # showBandId = FALSE avoids a separate label-overflow bug.
  cyto_chr <- .get_cytoband(genome, TARGET$chr, cytoband_file)

  itrack <- tryCatch({
    Gviz::IdeogramTrack(
      genome     = genome,
      chromosome = TARGET$chr,
      bands      = cyto_chr,
      showBandId = FALSE,
      fontcolor  = palette$axis,
      fontsize   = 9
    )
  }, error = function(e) {
    message("IdeogramTrack failed (skipped): ", conditionMessage(e))
    NULL
  })

  gtrack <- Gviz::GenomeAxisTrack(
    col = palette$axis, fontcolor = palette$axis, fontsize = 10
  )

  # ── Signal track ────────────────────────────────────────────
  signal_track <- Gviz::DataTrack(
    range            = bigwig_file,
    genome           = genome,
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

  # ── Peaks in view ───────────────────────────────────────────
  peak_region   <- GenomicRanges::GRanges(
    TARGET$chr, IRanges::IRanges(view_start, view_end)
  )
  peaks_in_view <- peaks[IRanges::overlapsAny(peaks, peak_region)]

  peak_track <- NULL
  if (length(peaks_in_view) > 0L) {
    peak_track <- Gviz::AnnotationTrack(
      peaks_in_view, name = "Peaks",
      fill = palette$peak, col = palette$peak, stacking = "dense",
      fontcolor.title = palette$axis, background.title = palette$bg,
      rotation.title = 0
    )
  }

  # ── TSS marker (±50 bp) ─────────────────────────────────────
  tss_track <- Gviz::AnnotationTrack(
    GenomicRanges::GRanges(
      TARGET$chr, IRanges::IRanges(TARGET$tss - 50, TARGET$tss + 50)
    ),
    name = "TSS", fill = palette$tss, col = palette$tss,
    fontcolor.title = palette$axis, background.title = palette$bg,
    rotation.title = 0
  )

  # ── Gene model ──────────────────────────────────────────────
  orgdb <- .load_orgdb(res$orgdb)

  gene_track <- Gviz::GeneRegionTrack(
    txdb,
    genome               = genome,
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

  # GeneRegionTrack from TxDb has no symbol column — Gviz falls back to
  # transcript ID (ENST...). Inject gene symbols by mapping Entrez -> symbol.
  entrez_ids <- as.character(gene_track@range$gene)
  if (length(entrez_ids) > 0 && any(!is.na(entrez_ids))) {
    sym_map <- suppressMessages(AnnotationDbi::select(
      orgdb,
      keys    = unique(stats::na.omit(entrez_ids)),
      columns = "SYMBOL",
      keytype = "ENTREZID"
    ))
    symbols <- sym_map$SYMBOL[match(entrez_ids, sym_map$ENTREZID)]
    symbols[is.na(symbols)] <- entrez_ids[is.na(symbols)]
    gene_track@range$symbol <- symbols
  }

  # ── Assemble track list ─────────────────────────────────────
  # core_tracks excludes itrack (added separately if available)
  if (!is.null(peak_track)) {
    core_tracks <- list(gtrack, signal_track, peak_track, tss_track, gene_track)
    core_sizes  <- c(1, 4, 0.5, 0.3, 1.5)
  } else {
    core_tracks <- list(gtrack, signal_track, tss_track, gene_track)
    core_sizes  <- c(1, 4, 0.3, 1.5)
  }

  if (!is.null(itrack)) {
    track_list  <- c(list(itrack), core_tracks)
    track_sizes <- c(0.5, core_sizes)
  } else {
    track_list  <- core_tracks
    track_sizes <- core_sizes
  }

  if (isTRUE(plot)) {
    # chromosome must be explicit to prevent Gviz from guessing wrong seqlevel
    Gviz::plotTracks(
      track_list,
      from             = view_start,
      to               = view_end,
      chromosome       = TARGET$chr,
      sizes            = track_sizes,
      background.panel = palette$bg,
      background.title = palette$bg,
      col.border.title = NA,
      main             = paste0(tf_name, " ChIP-seq at ",
                                TARGET$symbol, " Promoter"),
      cex.main         = 1.4,
      fontface.main    = 2,
      col.main         = palette$axis
    )
  }

  invisible(list(tracks = track_list, sizes = track_sizes,
                 from = view_start, to = view_end, target = TARGET))
}
