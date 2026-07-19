#' Plot a ChIP-seq Gviz track around a selected gene's TSS (Figure A)
#'
#' Builds a multi-track Gviz figure (ideogram, axis, signal, peaks, TSS,
#' gene model) centered on the TSS of \code{gene_symbol}. Supports human
#' (hg38 / hg19) and mouse (mm10 / mm39) assemblies.
#'
#' Gene symbol labels are placed on the upstream side of each TSS:
#' \code{+}-strand genes get labels on the left of the gene body, and
#' \code{-}-strand genes get labels on the right. Internally this is done
#' by splitting the gene model into two sub-tracks with different
#' \code{just.group} settings, since Gviz's \code{just.group} is a
#' track-level (not per-feature) parameter.
#'
#' @param bigwig_file Path to a bigWig signal file (\code{.bw} /
#'   \code{.bigWig}).
#' @param peaks_file Path to a narrowPeak or broadPeak file.
#' @param gene_symbol Gene symbol of the selected gene (e.g. \code{"KIF18A"}
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
  orgdb  <- .load_orgdb(res$orgdb)

  message("Genome: ", genome)
  message("Transcription factor: ", tf_name)
  message("Selected gene: ", TARGET$symbol, "  ", TARGET$chr, ":",
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

  # ── Inject gene symbols ─────────────────────────────────────
  # GeneRegionTrack built from a TxDb stores only Entrez IDs in @range$gene,
  # so transcriptAnnotation = "symbol" falls back to displaying transcript
  # IDs (ENST.../ENSMUST...). Map Entrez -> SYMBOL via OrgDb and write the
  # result to @range$symbol so labels render as e.g. "KIF18A" / "METTL15".
  entrez_ids <- as.character(gene_track@range$gene)
  if (length(entrez_ids) > 0L && any(!is.na(entrez_ids))) {
    sym_lookup <- tryCatch(
      AnnotationDbi::mapIds(
        orgdb,
        keys      = unique(stats::na.omit(entrez_ids)),
        column    = "SYMBOL",
        keytype   = "ENTREZID",
        multiVals = "first"
      ),
      error = function(e) {
        message("Symbol lookup failed (using Entrez IDs instead): ",
                conditionMessage(e))
        NULL
      }
    )
    if (!is.null(sym_lookup)) {
      symbols <- unname(sym_lookup[entrez_ids])
      # Fall back to Entrez ID when a symbol is missing
      symbols[is.na(symbols)] <- entrez_ids[is.na(symbols)]
      gene_track@range$symbol <- symbols
    }
  }

  # ── Strand-aware label placement ────────────────────────────
  # Gviz's `just.group` is a track-level parameter and can't be set per
  # gene. To place each gene's symbol on its upstream side, we split
  # gene_track into two sub-tracks by strand:
  #   + strand → TSS at left of gene body  → just.group = "left"
  #   - strand → TSS at right of gene body → just.group = "right"
  strand_per_exon <- as.character(
    GenomicRanges::strand(gene_track@range)
  )

  gene_track_plus <- gene_track
  gene_track_plus@range <- gene_track@range[strand_per_exon == "+"]
  Gviz::displayPars(gene_track_plus)$just.group <- "left"
  gene_track_plus@name <- ""        # blank title to avoid duplicate "Genes"

  gene_track_minus <- gene_track
  gene_track_minus@range <- gene_track@range[strand_per_exon == "-"]
  Gviz::displayPars(gene_track_minus)$just.group <- "right"
  gene_track_minus@name <- "Genes"

  # Drop empty strand tracks so we don't waste a row when only one strand
  # has genes in the view window.
  gene_tracks <- list()
  gene_sizes  <- numeric(0)
  if (length(gene_track_minus@range) > 0L) {
    gene_tracks <- c(gene_tracks, list(gene_track_minus))
    gene_sizes  <- c(gene_sizes, 0.8)
  }
  if (length(gene_track_plus@range) > 0L) {
    gene_tracks <- c(gene_tracks, list(gene_track_plus))
    gene_sizes  <- c(gene_sizes, 0.8)
  }
  # Fallback: no genes in view → keep combined track so plotTracks
  # still has a "Genes" row
  if (length(gene_tracks) == 0L) {
    gene_tracks <- list(gene_track)
    gene_sizes  <- 1.5
  }

  # ── Assemble track list ─────────────────────────────────────
  if (!is.null(peak_track)) {
    core_tracks <- c(list(gtrack, signal_track, peak_track, tss_track),
                     gene_tracks)
    core_sizes  <- c(1, 4, 0.5, 0.3, gene_sizes)
  } else {
    core_tracks <- c(list(gtrack, signal_track, tss_track),
                     gene_tracks)
    core_sizes  <- c(1, 4, 0.3, gene_sizes)
  }

  if (!is.null(itrack)) {
    track_list  <- c(list(itrack), core_tracks)
    track_sizes <- c(0.5, core_sizes)
  } else {
    track_list  <- core_tracks
    track_sizes <- core_sizes
  }

  if (isTRUE(plot)) {
    Gviz::plotTracks(
      track_list,
      from             = view_start,
      to               = view_end,
      chromosome       = TARGET$chr,
      sizes            = track_sizes,
      background.panel = palette$bg,
      background.title = palette$bg,
      col.border.title = NA,
      main             = paste0(tf_name, " ChIP-seq near ",
                                TARGET$symbol, " TSS"),
      cex.main         = 1.4,
      fontface.main    = 2,
      col.main         = palette$axis
    )
  }

  invisible(list(tracks = track_list, sizes = track_sizes,
                 from = view_start, to = view_end, target = TARGET))
}
