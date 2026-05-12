#' Plot a ChIP-seq EnrichedHeatmap across all target gene TSSs (Figure B)
#'
#' Annotates peaks to nearest genes via \pkg{ChIPseeker}, builds a
#' TSS-centered \pkg{EnrichedHeatmap} of the bigWig signal across all target
#' genes, ranks them by total signal, and highlights the query gene.
#' Supports human (hg38 / hg19) and mouse (mm10 / mm39) assemblies.
#'
#' @inheritParams plot_chip_track
#' @param tss_window Half-window (bp) around each TSS for the heatmap matrix
#'   (default \code{2000}). Axis labels ("-Xkb" / "+Xkb") are generated
#'   automatically from this value.
#' @param tss_region 2-element integer vector defining the promoter window
#'   used by \pkg{ChIPseeker} when assigning peaks to genes
#'   (default \code{c(-2000, 2000)}).
#' @param plot Logical. If \code{TRUE} (default) the heatmap is drawn.
#'
#' @return Invisibly, a list with components:
#' \describe{
#'   \item{matrix}{Signal matrix (rows ordered by total signal, decreasing).}
#'   \item{target_rank}{Integer rank of the query gene (\code{NA} if absent).}
#'   \item{n_targets}{Number of TSSs in the heatmap.}
#'   \item{tss_gr}{The \code{GRanges} of TSSs used.}
#' }
#'
#' @details
#' **Strand handling**: \code{normalizeToMatrix} respects the \code{strand}
#' slot of the \code{target} \code{GRanges}. For negative-strand genes the
#' signal window is automatically reversed so that the left side of every
#' heatmap row corresponds to upstream sequence, making the "-Xkb" / "+Xkb"
#' axis labels biologically consistent across all rows regardless of strand.
#'
#' @export
plot_chip_heatmap <- function(bigwig_file,
                              peaks_file,
                              gene_symbol,
                              tf_name,
                              genome     = "hg38",
                              tss_window = 2000,
                              tss_region = c(-2000, 2000),
                              palette    = morandi_palette(),
                              plot       = TRUE) {

  TARGET <- get_gene_info(gene_symbol, genome = genome)
  res    <- .genome_resources(genome)
  txdb   <- .load_txdb(res$txdb)
  orgdb  <- .load_orgdb(res$orgdb)
  peaks  <- .read_peaks(peaks_file)

  # ── Peak → gene annotation ──────────────────────────────────
  message("Annotating peaks to genes (ChIPseeker)...")
  peakAnno <- ChIPseeker::annotatePeak(
    peaks, tssRegion = tss_region,
    TxDb = txdb, annoDb = res$orgdb, verbose = FALSE
  )

  target_genes <- unique(stats::na.omit(as.data.frame(peakAnno)$SYMBOL))
  # Always include the query gene even if it has no local peak
  if (!TARGET$symbol %in% target_genes)
    target_genes <- c(target_genes, TARGET$symbol)

  # ── TSS GRanges: one per gene, standard chromosomes only ────
  all_tss <- GenomicFeatures::promoters(txdb, upstream = 0, downstream = 1)
  all_tss <- all_tss[!duplicated(all_tss$tx_name)]
  all_tss <- all_tss[as.character(GenomeInfoDb::seqnames(all_tss)) %in% res$std_chr]
  all_tss <- GenomeInfoDb::keepSeqlevels(all_tss, res$std_chr,
                                         pruning.mode = "coarse")

  # TX → Entrez → Symbol mapping
  tx2gene <- AnnotationDbi::select(
    txdb, keys = all_tss$tx_name,
    columns = c("TXNAME", "GENEID"), keytype = "TXNAME"
  )
  gene2symbol <- AnnotationDbi::select(
    orgdb, keys = unique(stats::na.omit(tx2gene$GENEID)),
    columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID"
  )
  tx2gene        <- merge(tx2gene, gene2symbol,
                          by.x = "GENEID", by.y = "ENTREZID", all.x = TRUE)
  all_tss$SYMBOL <- tx2gene$SYMBOL[match(all_tss$tx_name, tx2gene$TXNAME)]

  all_tss_filtered <- all_tss[all_tss$SYMBOL %in% target_genes]
  all_tss_filtered <- all_tss_filtered[!duplicated(all_tss_filtered$SYMBOL)]

  message(tf_name, " target genes: ", length(target_genes))
  message("TSSs in heatmap: ", length(all_tss_filtered))

  # ── Signal matrix ───────────────────────────────────────────
  # normalizeToMatrix flips the window for negative-strand targets,
  # so upstream is always on the left of the resulting matrix
  message("Computing signal matrix (strand-aware)...")
  bw  <- rtracklayer::import(bigwig_file, format = "BigWig")
  mat <- EnrichedHeatmap::normalizeToMatrix(
    signal       = bw,
    target       = all_tss_filtered,
    extend       = tss_window,
    mean_mode    = "w0",
    value_column = "score",
    background   = 0,
    smooth       = TRUE
  )

  # Sort rows by total signal (highest → lowest)
  order_idx   <- order(rowSums(mat, na.rm = TRUE), decreasing = TRUE)
  mat_ordered <- mat[order_idx, ]

  target_idx  <- which(all_tss_filtered$SYMBOL == TARGET$symbol)
  target_rank <- if (length(target_idx) > 0L) {
    rk <- which(order_idx == target_idx[1])
    message(TARGET$symbol, " rank: ", rk, "/", nrow(mat),
            " (top ", round(rk / nrow(mat) * 100, 1), "%)")
    rk
  } else {
    NA_integer_
  }

  # ── Color ramp: 0 → median → 95th percentile ────────────────
  col_fun <- circlize::colorRamp2(
    c(0,
      stats::quantile(mat_ordered, 0.5,  na.rm = TRUE),
      stats::quantile(mat_ordered, 0.95, na.rm = TRUE)),
    c(palette$heatmap_low, palette$heatmap_mid, palette$heatmap_high)
  )

  # ── Axis labels derived from tss_window ─────────────────────
  lbl_left  <- paste0("-", .bp_label(tss_window))
  lbl_right <- paste0("+", .bp_label(tss_window))

  # ── Heatmap ─────────────────────────────────────────────────
  # axis_name_gp col="transparent": hide EnrichedHeatmap's default axis
  # labels ("-2000 / start / 2000") while keeping their layout space,
  # then draw our own labels via decorate_heatmap_body below.
  ht <- EnrichedHeatmap::EnrichedHeatmap(
    mat_ordered,
    name            = tf_name,
    col             = col_fun,
    axis_name_gp    = grid::gpar(col = "transparent", fontsize = 0.001),
    top_annotation  = ComplexHeatmap::HeatmapAnnotation(
      enriched = EnrichedHeatmap::anno_enriched(
        gp     = grid::gpar(col  = palette$signal,
                            lwd  = 2,
                            fill = palette$signal_fill),
        height = grid::unit(3, "cm")
      ),
      show_annotation_name = FALSE
    ),
    column_title    = paste0(tf_name, " Binding at Target Gene TSS (",
                             TARGET$symbol, " highlighted)"),
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold",
                                 col = palette$axis),
    show_row_names  = FALSE,
    heatmap_legend_param = list(
      title         = tf_name,
      title_gp      = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp     = grid::gpar(fontsize = 9),
      legend_height = grid::unit(4, "cm")
    ),
    border         = FALSE,
    use_raster     = TRUE,
    raster_quality = 3,
    pos_line       = FALSE
  )

  # Right-side annotation marking the query gene's row
  if (!is.na(target_rank)) {
    ht <- ht + ComplexHeatmap::rowAnnotation(
      mark = ComplexHeatmap::anno_mark(
        at        = target_rank,
        labels    = TARGET$symbol,
        labels_gp = grid::gpar(fontsize = 10, fontface = "bold",
                               col = palette$peak),
        link_gp   = grid::gpar(col = palette$peak, lwd = 1.5)
      ),
      show_annotation_name = FALSE
    )
  }

  if (isTRUE(plot)) {
    ComplexHeatmap::draw(ht, padding = grid::unit(c(5, 5, 5, 5), "mm"))

    # Draw custom x-axis labels (replaces the suppressed defaults)
    ComplexHeatmap::decorate_heatmap_body(tf_name, {
      for (pos in list(c(0, lbl_left, "plain"),
                       c(0.5, "TSS",     "bold"),
                       c(1,   lbl_right, "plain"))) {
        grid::grid.text(
          pos[[2]],
          x    = grid::unit(as.numeric(pos[[1]]), "npc"),
          y    = grid::unit(-2, "mm"),
          just = "top",
          gp   = grid::gpar(fontsize = 9, col = palette$axis,
                            fontface = pos[[3]])
        )
      }
    })

    # Top-left annotation: query gene rank
    if (!is.na(target_rank)) {
      grid::grid.text(
        paste0(TARGET$symbol, ": #", target_rank, "/", nrow(mat),
               " (top ", round(target_rank / nrow(mat) * 100, 1), "%)"),
        x    = grid::unit(0.02, "npc"),
        y    = grid::unit(0.98, "npc"),
        just = c("left", "top"),
        gp   = grid::gpar(fontsize = 9, col = palette$peak, fontface = "bold")
      )
    }
  }

  invisible(list(matrix      = mat_ordered,
                 target_rank = target_rank,
                 n_targets   = nrow(mat),
                 tss_gr      = all_tss_filtered))
}
