#' Plot a ChIP-seq EnrichedHeatmap across all target gene TSSs (Figure B)
#'
#' Annotates peaks to nearest genes via \pkg{ChIPseeker}, builds a TSS-centered
#' \pkg{EnrichedHeatmap} (\eqn{\pm}2 kb) of the bigWig signal across all target
#' genes, ranks them by total signal, and highlights the query gene.
#'
#' @inheritParams plot_chip_track
#' @param tss_window Half-window (bp) around each TSS for the heatmap matrix
#'   (default \code{2000}).
#' @param tss_region 2-element vector defining the promoter region used by
#'   \pkg{ChIPseeker} when annotating peaks to genes
#'   (default \code{c(-2000, 2000)}).
#' @param plot Logical. If \code{TRUE} (default) the heatmap is drawn on the
#'   active graphics device.
#'
#' @return Invisibly, a list with components:
#' \describe{
#'   \item{matrix}{The signal matrix (rows ordered by total signal, decreasing).}
#'   \item{target_rank}{Integer rank of the query gene (\code{NA} if absent).}
#'   \item{n_targets}{Number of TSSs in the heatmap.}
#'   \item{tss_gr}{The \code{GRanges} of TSSs used.}
#' }
#' @export
plot_chip_heatmap <- function(bigwig_file,
                              peaks_file,
                              gene_symbol,
                              tf_name,
                              tss_window = 2000,
                              tss_region = c(-2000, 2000),
                              palette    = morandi_palette(),
                              plot       = TRUE) {

  TARGET <- get_gene_info(gene_symbol)

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  peaks <- .read_peaks(peaks_file)

  message("Annotating peaks to genes (ChIPseeker)...")
  peakAnno <- ChIPseeker::annotatePeak(
    peaks,
    tssRegion = tss_region,
    TxDb      = txdb,
    annoDb    = "org.Hs.eg.db",
    verbose   = FALSE
  )

  target_genes <- unique(stats::na.omit(as.data.frame(peakAnno)$SYMBOL))
  if (!TARGET$symbol %in% target_genes) {
    target_genes <- c(target_genes, TARGET$symbol)
  }

  ## All TSSs (one per transcript) ----------------------------------------
  all_tss <- GenomicFeatures::promoters(txdb, upstream = 0, downstream = 1)
  all_tss <- all_tss[!duplicated(all_tss$tx_name)]

  standard_chr <- paste0("chr", c(1:22, "X", "Y"))
  all_tss <- all_tss[as.character(GenomeInfoDb::seqnames(all_tss)) %in% standard_chr]
  all_tss <- GenomeInfoDb::keepSeqlevels(all_tss, standard_chr,
                                         pruning.mode = "coarse")

  ## TX -> Entrez -> Symbol -----------------------------------------------
  tx2gene <- AnnotationDbi::select(
    txdb,
    keys     = all_tss$tx_name,
    columns  = c("TXNAME", "GENEID"),
    keytype  = "TXNAME"
  )
  gene2symbol <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys     = unique(stats::na.omit(tx2gene$GENEID)),
    columns  = c("ENTREZID", "SYMBOL"),
    keytype  = "ENTREZID"
  )
  tx2gene <- merge(tx2gene, gene2symbol,
                   by.x = "GENEID", by.y = "ENTREZID", all.x = TRUE)
  all_tss$SYMBOL <- tx2gene$SYMBOL[match(all_tss$tx_name, tx2gene$TXNAME)]

  all_tss_filtered <- all_tss[all_tss$SYMBOL %in% target_genes]
  all_tss_filtered <- all_tss_filtered[!duplicated(all_tss_filtered$SYMBOL)]

  message(tf_name, " target genes: ", length(target_genes))
  message("TSSs in heatmap: ", length(all_tss_filtered))

  ## Signal matrix ---------------------------------------------------------
  message("Computing signal matrix...")
  bw <- rtracklayer::import(bigwig_file, format = "BigWig")
  mat <- EnrichedHeatmap::normalizeToMatrix(
    signal       = bw,
    target       = all_tss_filtered,
    extend       = tss_window,
    mean_mode    = "w0",
    value_column = "score",
    background   = 0,
    smooth       = TRUE
  )

  ## Sort rows by total signal -------------------------------------------
  row_signals <- rowSums(mat, na.rm = TRUE)
  order_idx   <- order(row_signals, decreasing = TRUE)
  mat_ordered <- mat[order_idx, ]

  target_idx <- which(all_tss_filtered$SYMBOL == TARGET$symbol)
  if (length(target_idx) > 0L) {
    target_rank <- which(order_idx == target_idx[1])
    message(TARGET$symbol, " rank among ", tf_name, " targets: ",
            target_rank, "/", nrow(mat),
            " (top ", round(target_rank / nrow(mat) * 100, 1), "%)")
  } else {
    target_rank <- NA_integer_
  }

  ## Color ramp -----------------------------------------------------------
  col_fun <- circlize::colorRamp2(
    c(0,
      stats::quantile(mat_ordered, 0.5,  na.rm = TRUE),
      stats::quantile(mat_ordered, 0.95, na.rm = TRUE)),
    c(palette$heatmap_low, palette$heatmap_mid, palette$heatmap_high)
  )

  ## Heatmap --------------------------------------------------------------
  ht <- EnrichedHeatmap::EnrichedHeatmap(
    mat_ordered,
    name            = tf_name,
    col             = col_fun,
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
    border        = FALSE,
    use_raster    = TRUE,
    raster_quality = 3,
    pos_line      = FALSE
  )

  ## Highlight target row -------------------------------------------------
  if (!is.na(target_rank)) {
    row_anno <- ComplexHeatmap::rowAnnotation(
      mark = ComplexHeatmap::anno_mark(
        at        = target_rank,
        labels    = TARGET$symbol,
        labels_gp = grid::gpar(fontsize = 10, fontface = "bold",
                               col = palette$peak),
        link_gp   = grid::gpar(col = palette$peak, lwd = 1.5)
      ),
      show_annotation_name = FALSE
    )
    ht <- ht + row_anno
  }

  if (isTRUE(plot)) {
    ComplexHeatmap::draw(ht, padding = grid::unit(c(5, 5, 5, 5), "mm"))

    ComplexHeatmap::decorate_heatmap_body(tf_name, {
      grid::grid.text("-2kb", x = grid::unit(0,   "npc"),
                      y = grid::unit(-2, "mm"), just = "top",
                      gp = grid::gpar(fontsize = 9, col = palette$axis))
      grid::grid.text("TSS",  x = grid::unit(0.5, "npc"),
                      y = grid::unit(-2, "mm"), just = "top",
                      gp = grid::gpar(fontsize = 9, col = palette$axis,
                                      fontface = "bold"))
      grid::grid.text("+2kb", x = grid::unit(1,   "npc"),
                      y = grid::unit(-2, "mm"), just = "top",
                      gp = grid::gpar(fontsize = 9, col = palette$axis))
    })

    if (!is.na(target_rank)) {
      grid::grid.text(
        paste0(TARGET$symbol, ": #", target_rank, "/", nrow(mat),
               " (top ", round(target_rank / nrow(mat) * 100, 1), "%)"),
        x  = grid::unit(0.02, "npc"),
        y  = grid::unit(0.98, "npc"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 9, col = palette$peak, fontface = "bold")
      )
    }
  }

  invisible(list(
    matrix      = mat_ordered,
    target_rank = target_rank,
    n_targets   = nrow(mat),
    tss_gr      = all_tss_filtered
  ))
}
