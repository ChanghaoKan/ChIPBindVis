#' Plot a TSS enrichment heatmap of TF ChIP-seq signal (Figure B)
#'
#' Builds an \code{\link[EnrichedHeatmap]{EnrichedHeatmap}} of TF ChIP-seq
#' signal centered on a gene-level TSS for genes with a promoter-associated
#' peak, sorted by total signal. If the user-specified \code{gene_symbol} has
#' a promoter-associated peak, it is highlighted and ranked. Otherwise it is
#' not added to the heatmap, and the result reports that status explicitly.
#' The x-axis is labelled \code{-Nkb / TSS / +Nkb}.
#'
#' Peak-to-gene annotations are obtained with
#' \code{ChIPseeker::annotatePeak}; only rows classified as promoter peaks
#' within \code{tss_region} contribute genes. Each gene contributes a
#' gene-level TSS derived from its TxDb gene range. This is not an isoform-
#' specific or canonical-transcript TSS assignment.
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
#'   query gene, or \code{NA}), \code{n_targets} (total rows),
#'   \code{query_promoter_peak_detected} (logical), \code{query_status}
#'   (character), and \code{target} (resolved gene info).
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
  message("Annotating peaks to extract promoter-associated gene set...")
  peakAnno <- ChIPseeker::annotatePeak(
    peaks,
    tssRegion = tss_region,
    TxDb      = txdb,
    annoDb    = res$orgdb,
    verbose   = FALSE
  )
  anno_df <- as.data.frame(peakAnno)
  target_status <- .promoter_target_status(anno_df, TARGET$symbol)
  target_genes <- target_status$target_genes
  query_promoter_peak_detected <-
    target_status$query_promoter_peak_detected

  if (!query_promoter_peak_detected) {
    message(
      "No promoter-associated peak was assigned to ", TARGET$symbol,
      " within tss_region = c(", paste(tss_region, collapse = ", "),
      "); the query gene will not be added to or ranked in the heatmap."
    )
  }

  # Step 2: derive one gene-level TSS per representable TxDb gene
  # GenomicFeatures::genes() supplies one range per representable gene.
  # Converting those ranges to one-base promoters yields a transparent
  # gene-level TSS proxy rather than arbitrarily selecting one transcript.
  all_genes <- suppressWarnings(GenomicFeatures::genes(txdb))
  # Keep only standard chromosomes for the chosen genome
  keep <- as.character(GenomeInfoDb::seqnames(all_genes)) %in% res$std_chr
  all_genes <- all_genes[keep]
  standard_levels <- intersect(
    res$std_chr, GenomeInfoDb::seqlevels(all_genes)
  )
  all_genes <- GenomeInfoDb::keepSeqlevels(
    all_genes, standard_levels, pruning.mode = "coarse"
  )

  if (length(all_genes) == 0L)
    stop("No gene ranges after chromosome filtering. Check that the TxDb ",
         "matches the genome assembly (", genome, ").", call. = FALSE)

  # Step 3: map Entrez gene identifiers to symbols
  entrez_per_gene <- as.character(all_genes$gene_id)
  if (length(entrez_per_gene) != length(all_genes))
    entrez_per_gene <- as.character(names(all_genes))
  if (length(entrez_per_gene) != length(all_genes) ||
      all(is.na(entrez_per_gene) | !nzchar(entrez_per_gene))) {
    stop("TxDb gene ranges do not expose usable gene identifiers for ",
         "symbol mapping.", call. = FALSE)
  }
  valid_entrez <- unique(
    entrez_per_gene[!is.na(entrez_per_gene) & nzchar(entrez_per_gene)]
  )
  symbol_per_entrez <- if (length(valid_entrez) > 0L) {
    suppressMessages(tryCatch(
      AnnotationDbi::mapIds(
        orgdb,
        keys      = valid_entrez,
        column    = "SYMBOL",
        keytype   = "ENTREZID",
        multiVals = "first"
      ),
      error = function(e) {
        message("Entrez -> symbol mapping failed: ", conditionMessage(e))
        stats::setNames(rep(NA_character_, length(valid_entrez)), valid_entrez)
      }
    ))
  } else {
    stats::setNames(character(0), character(0))
  }
  all_genes$SYMBOL <- unname(symbol_per_entrez[entrez_per_gene])
  all_tss <- GenomicFeatures::promoters(
    all_genes, upstream = 0, downstream = 1
  )
  all_tss$SYMBOL <- all_genes$SYMBOL

  # Filter to promoter-associated genes, one gene-level TSS per symbol.
  tss_use <- all_tss[!is.na(all_tss$SYMBOL) &
                       all_tss$SYMBOL %in% target_genes]
  tss_use <- tss_use[!duplicated(tss_use$SYMBOL)]

  message("Genes with promoter-associated peaks: ", length(target_genes),
          " | gene-level TSS rows used: ", length(tss_use))

  if (length(tss_use) == 0L)
    stop("No genes with promoter-associated peaks could be mapped to TxDb ",
         "gene-level TSSs. Check peaks_file, tss_region, genome, and ",
         "annotation compatibility.", call. = FALSE)

  # ── Step 4: build signal matrix around each TSS ──────────────
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

  # Locate the query gene after sorting. It is never inserted when absent.
  query_idx  <- which(tss_use$SYMBOL == TARGET$symbol)
  query_rank <- if (length(query_idx) > 0L) {
    which(order_idx == query_idx[1])
  } else NA_integer_

  query_status <- if (!query_promoter_peak_detected) {
    "no_promoter_peak_detected"
  } else if (is.na(query_rank)) {
    "promoter_peak_detected_not_ranked"
  } else {
    "promoter_peak_detected_and_ranked"
  }

  if (!is.na(query_rank))
    message(TARGET$symbol, " rank: ", query_rank, "/", nrow(mat),
            " (top ", round(query_rank / nrow(mat) * 100, 1), "%)")
  if (identical(query_status, "promoter_peak_detected_not_ranked"))
    message(TARGET$symbol, " had a promoter-associated peak but could not ",
            "be mapped to a gene-level TxDb TSS, so no rank is reported.")

  # ── Step 5: color mapping (0 / 50% / 95% three-point ramp) ───
  col_fun <- circlize::colorRamp2(
    c(0,
      stats::quantile(mat_ordered, 0.5,  na.rm = TRUE),
      stats::quantile(mat_ordered, 0.95, na.rm = TRUE)),
    c(palette$heatmap_low, palette$heatmap_mid, palette$heatmap_high)
  )

  # ── Step 6: assemble EnrichedHeatmap ─────────────────────────
  # axis_name_gp = transparent: we draw our own -Nkb / TSS / +Nkb labels
  # via decorate_heatmap_body so label spacing scales with tss_window
  query_label <- if (identical(
    query_status, "promoter_peak_detected_and_ranked"
  )) {
    paste0(TARGET$symbol, " highlighted")
  } else if (identical(
    query_status, "promoter_peak_detected_not_ranked"
  )) {
    paste0(TARGET$symbol, ": promoter peak detected; TSS not ranked")
  } else {
    paste0(TARGET$symbol, ": no promoter-associated peak detected")
  }

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
    column_title = paste0(tf_name, " ChIP-seq Signal at Gene-level TSSs (",
                          query_label, ")"),
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold",
                                 col = palette$axis),
    show_row_names = FALSE,
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

  # External link to the highlighted target row
  if (!is.na(query_rank)) {
    ht <- ht + ComplexHeatmap::rowAnnotation(
      mark = ComplexHeatmap::anno_mark(
        at = query_rank, labels = TARGET$symbol,
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
    if (!is.na(query_rank)) {
      grid::grid.text(
        paste0(TARGET$symbol, ": #", query_rank, "/", nrow(mat),
               " (top ", round(query_rank / nrow(mat) * 100, 1), "%)"),
        x = grid::unit(0.02, "npc"),
        y = grid::unit(0.98, "npc"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 9, col = palette$peak, fontface = "bold")
      )
    }
  }

  invisible(list(
    matrix                       = mat_ordered,
    order                        = order_idx,
    rank                         = query_rank,
    n_targets                    = nrow(mat),
    query_promoter_peak_detected = query_promoter_peak_detected,
    query_status                 = query_status,
    target                       = TARGET
  ))
}
