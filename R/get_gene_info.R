#' Resolve an HGNC gene symbol to coordinates (hg38)
#'
#' Looks up an HGNC gene symbol in \pkg{org.Hs.eg.db} and
#' \pkg{TxDb.Hsapiens.UCSC.hg38.knownGene}, and returns its chromosome,
#' start, end, strand, and TSS position.
#'
#' @param gene_symbol Character scalar. HGNC gene symbol (e.g. \code{"KIF18A"}).
#'
#' @return A named list with elements \code{symbol}, \code{chr},
#'   \code{start}, \code{end}, \code{strand}, \code{tss}.
#' @export
#'
#' @examples
#' \dontrun{
#' info <- get_gene_info("KIF18A")
#' info$tss
#' }
get_gene_info <- function(gene_symbol) {
  stopifnot(is.character(gene_symbol), length(gene_symbol) == 1L)

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

  gene_info <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys     = gene_symbol,
    columns  = "ENTREZID",
    keytype  = "SYMBOL"
  )
  if (nrow(gene_info) == 0L || is.na(gene_info$ENTREZID[1])) {
    stop("Gene not found: ", gene_symbol, call. = FALSE)
  }

  target_entrez <- gene_info$ENTREZID[1]
  all_genes     <- GenomicFeatures::genes(txdb)
  target_gene   <- all_genes[all_genes$gene_id == target_entrez]

  if (length(target_gene) == 0L) {
    stop("No coordinates for ", gene_symbol,
         " (ENTREZ ID ", target_entrez, ") in TxDb.", call. = FALSE)
  }

  chr_   <- as.character(GenomeInfoDb::seqnames(target_gene))
  start_ <- GenomicRanges::start(target_gene)
  end_   <- GenomicRanges::end(target_gene)
  str_   <- as.character(GenomicRanges::strand(target_gene))
  tss_   <- if (str_ == "-") end_ else start_

  list(
    symbol = gene_symbol,
    chr    = chr_,
    start  = start_,
    end    = end_,
    strand = str_,
    tss    = tss_
  )
}
