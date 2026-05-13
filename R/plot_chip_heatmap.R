#' Plot a TSS enrichment heatmap of TF ChIP-seq signal (Figure B)
#'
#' Builds an \code{\link[EnrichedHeatmap]{EnrichedHeatmap}} of TF ChIP-seq
#' signal centered on the TSS of every peak-annotated target gene, sorted
#' by total signal. The user-specified \code{gene_symbol} is highlighted
#' with an external label, its rank is annotated in the top-left corner,
#' and the x-axis is labelled \code{-Nkb / TSS / +Nkb}.
#'
#' Target genes are obtained via \code{ChIPseeker::annotatePeak} on the
#' peak file. Each target gene contributes one canonical TSS (duplicates
#' collapsed by SYMBOL).
#'
#' @inheritParams plot_chip_track
#' @param tss_window Half-window in bp around each TSS for the signal
#'   matrix (default \code{2000}).
#' @param tss_region Two-element numeric vector passed to
#'   \code{ChIPseeker::annotatePeak} for peak-to-gene assignment
#'   (default \code{c(-2000, 2000)}).
#' @param plot Logical. If \code{TRUE} (default) the heatmap is drawn on
#'   the active graphics device.
#'
#' @return Invisibly, a list with \code{matrix} (sorted normalized
#'   matrix), \code{order} (row order indices), \code{rank} (rank of the
#'   target gene), \code{n_targets} (total rows), and \code{target}
#'   (resolved gene info).
#' @export
#'
#' @examples
#' \dontrun{
#' plot_chip_heatmap("E2F8.bigWig", "E2F8.narrowPeak", "KIF18A", "E2F8")
#' }
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

  peaks <- .read_peaks(peaks_file)

  message("Genome: ", genome,
          " | TF: ", tf_name,
          " | Target: ", TARGET$symbol)

  # ── Step 1: peak → gene annotation ───────────────────────────
  message("Annotating peaks to extract target gene set...")
  peakAnno <- ChIPseeker::annotatePeak(
    peaks,
    tssRegion = tss_region,
    TxDb      = txdb,
    annoDb    = res$orgdb,
    verbose   = FALSE
  )
  target_genes <- unique(stats::na.omit(as.data.frame(peakAnno)$SYMBOL))
  # Ensure the highlighted target is in the set even if it has no direct peak
  if (!TARGET$symbol %in% target_genes)
    target_genes <- c(target_genes, TARGET$symbol)

  # ── Step 2: all TSS → filter to target genes → de-duplicate ──
  all_tss <- GenomicFeatures::promoters(txdb, upstream = 0, downstream = 1)
  all_tss <- all_tss[!duplicated(all_tss$tx_name)]
  # Keep only standard chromosomes for the chosen genome
  all_tss <- all_tss[GenomeInfoDb::seqnames(all_tss) %in% res$std_chr]
  all_tss <- GenomeInfoDb::keepSeqlevels(
    all_tss, res$std_chr, pruning.mode = "coarse"
  )

  # Map transcript -> gene -> symbol
  tx2gene <- AnnotationDbi::select(
    txdb,
    keys    = all_tss$tx_name,
    columns = c("TXNAME", "GENEID"),
    keytype = "TXNAME"
  )
  gene2sym <- AnnotationDbi::select(
    orgdb,
    keys    = unique(stats::na.omit(tx2gene$GENEID)),
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "ENTREZID"
  )
  tx2gene <- merge(tx2gene, gene2sym,
                   by.x = "GENEID", by.y = "ENTREZID", all.x = TRUE)
  all_tss$SYMBOL <- tx2gene$SYMBOL[match(all_tss$tx_name, tx2gene$TXNAME)]

  tss_use <- all_tss[all_tss$SYMBOL %in% target_genes]
  tss_use <- tss_use[!duplicated(tss_use$SYMBOL)]

  message("Targets: ", length(target_genes),
          " | TSS rows used: ", length(tss_use))

  if (length(tss_use) == 0L)
    stop("No TSS rows passed filtering. Check that peaks_file contains ",
         "binding sites overlapping known TSSs.", call. = FALSE)

  # ── Step 3: build signal matrix around each TSS ──────────────
  bw  <- rtracklayer::import(bigwig_file, format = "BigWig")
  mat <- EnrichedHeatmap::normalizeToMatrix(
    signal       = bw,
    target       = tss_use,
    extend       = tss_window,
    mean_mode    = "w0",
    value_column = "score",
    background   = 0,
    smooth       = TRUE
  )
  order_idx   <- order(rowSums(mat, na.rm = TRUE), decreasing = TRUE)
  mat_ordered <- mat[order_idx, ]

  # Locate the highlighted target's row after sorting
  target_idx  <- which(tss_use$SYMBOL == TARGET$symbol)
  target_rank <- if (length(target_idx) > 0L) {
    which(order_idx == target_idx[1])
  } else NA_integer_

  if (!is.na(target_rank))
    message(TARGET$symbol, " rank: ", target_rank, "/", nrow(mat),
            " (top ", round(target_rank / nrow(mat) * 100, 1), "%)")

  # ── Step 4: color mapping (0 / 50% / 95% three-point ramp) ───
  col_fun <- circlize::colorRamp2(
    c(0,
      stats::quantile(mat_ordered, 0.5,  na.rm = TRUE),
      stats::quantile(mat_ordered, 0.95, na.rm = TRUE)),
    c(palette$heatmap_low, palette$heatmap_mid, palette$heatmap_high)
  )

  # ── Step 5: assemble EnrichedHeatmap ─────────────────────────
  # axis_name_gp = transparent: we draw our own -Nkb / TSS / +Nkb labels
  # via decorate_heatmap_body so that label spacing scales with tss_window
  ht <- EnrichedHeatmap::EnrichedHeatmap(
    mat_ordered,
    name = tf_name, col = col_fun,
    axis_name_gp = grid::gpar(col = "transparent", fontsize = 0.001),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      enriched = EnrichedHeatmap::anno_enriched(
        gp = grid::gpar(col = palette$signal, lwd = 2,
                        fill = palette$signal_fill),
        height = grid::unit(3, "cm")
      ),
      show_annotation_name = FALSE
    ),
    column_title = paste0(tf_name, " Binding at Target Gene TSS (",
                          TARGET$symbol, " highlighted)"),
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold",
                                 col = palette$axis),
    show_row_names = FALSE,
    heatmap_legend_param = list(
      title         = tf_name,
      title_gp      = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp     = grid::gpar(fontsize = 9),
      legend_height = grid::unit(4, "cm")
    ),
    border        = FALSE,
    use_raster    = TRUE,
    raster_quality = 3,
    pos_line      = FALSE
  )

  # Mark the highlighted target with an external link annotation
  if (!is.na(target_rank)) {
    ht <- ht + ComplexHeatmap::rowAnnotation(
      mark = ComplexHeatmap::anno_mark(
        at = target_rank, labels = TARGET$symbol,
        labels_gp = grid::gpar(fontsize = 10, fontface = "bold",
                               col = palette$peak),
        link_gp = grid::gpar(col = palette$peak, lwd = 1.5)
      ),
      show_annotation_name = FALSE
    )
  }

  if (isTRUE(plot)) {
    ComplexHeatmap::draw(ht, padding = grid::unit(c(5, 5, 5, 5), "mm"))

    # Custom x-axis: -Nkb / TSS / +Nkb scaled to tss_window
    win_lab <- .bp_label(tss_window)
    ComplexHeatmap::decorate_heatmap_body(tf_name, {
      labels <- c(paste0("-", win_lab), "TSS", paste0("+", win_lab))
      for (i in seq_along(labels)) {
        grid::grid.text(
          labels[i],
          x = grid::unit(c(0, 0.5, 1)[i], "npc"),
          y = grid::unit(-2, "mm"),
          just = "top",
          gp = grid::gpar(
            fontsize = 9, col = palette$axis,
            fontface = if (labels[i] == "TSS") "bold" else "plain"
          )
        )
      }
    })

    # Top-left rank stamp
    if (!is.na(target_rank)) {
      grid::grid.text(
        paste0(TARGET$symbol, ": #", target_rank, "/", nrow(mat),
               " (top ", round(target_rank / nrow(mat) * 100, 1), "%)"),
        x = grid::unit(0.02, "npc"),
        y = grid::unit(0.98, "npc"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 9, col = palette$peak, fontface = "bold")
      )
    }
  }

  invisible(list(
    matrix    = mat_ordered,
    order     = order_idx,
    rank      = target_rank,
    n_targets = nrow(mat),
    target    = TARGET
  ))
}
