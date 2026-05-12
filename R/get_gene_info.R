#' Resolve a gene symbol to genomic coordinates
#'
#' Looks up a gene symbol in the appropriate annotation databases for the
#' chosen genome assembly and returns its chromosome, start, end, strand,
#' and TSS position.
#'
#' @param gene_symbol Character scalar. Gene symbol (e.g. \code{"KIF18A"} for
#'   human, \code{"Kif18a"} for mouse — case must match the database).
#' @param genome Character scalar. Reference genome assembly. One of
#'   \code{"hg38"} (default), \code{"hg19"}, \code{"mm10"}, \code{"mm39"}.
#'
#' @return A named list with elements \code{symbol}, \code{chr},
#'   \code{start}, \code{end}, \code{strand}, \code{tss}, \code{genome}.
#'
#' @details
#' **Strand and TSS**: For positive-strand (\code{+}) genes the TSS equals
#' the gene \code{start}; for negative-strand (\code{-}) genes it equals
#' \code{end}. Downstream functions center their views on this coordinate,
#' and \pkg{EnrichedHeatmap}'s \code{normalizeToMatrix} automatically
#' reverses the signal matrix for negative-strand targets so that the
#' left side of every heatmap row always represents upstream sequence.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Human gene (hg38)
#' get_gene_info("KIF18A")
#'
#' # Mouse gene (mm10)
#' get_gene_info("Kif18a", genome = "mm10")
#' }
get_gene_info <- function(gene_symbol, genome = "hg38") {
  stopifnot(is.character(gene_symbol), length(gene_symbol) == 1L)

  res   <- .genome_resources(genome)
  txdb  <- .load_txdb(res$txdb)
  orgdb <- .load_orgdb(res$orgdb)

  gene_info <- AnnotationDbi::select(
    orgdb,
    keys    = gene_symbol,
    columns = "ENTREZID",
    keytype = "SYMBOL"
  )
  if (nrow(gene_info) == 0L || is.na(gene_info$ENTREZID[1])) {
    stop("Gene not found: '", gene_symbol, "' in ", res$orgdb,
         ".\nFor mouse, use the mouse-style capitalisation (e.g. 'Kif18a').",
         call. = FALSE)
  }

  target_entrez <- gene_info$ENTREZID[1]
  all_genes     <- GenomicFeatures::genes(txdb)
  target_gene   <- all_genes[all_genes$gene_id == target_entrez]

  if (length(target_gene) == 0L)
    stop("No coordinates for '", gene_symbol, "' (Entrez ", target_entrez,
         ") in ", res$txdb, ".", call. = FALSE)

  chr_   <- as.character(GenomeInfoDb::seqnames(target_gene))
  start_ <- GenomicRanges::start(target_gene)
  end_   <- GenomicRanges::end(target_gene)
  str_   <- as.character(GenomicRanges::strand(target_gene))

  # TSS = start for + strand, end for - strand
  tss_   <- if (str_ == "-") end_ else start_

  list(symbol = gene_symbol, chr = chr_, start = start_, end = end_,
       strand = str_, tss = tss_, genome = genome)
}
