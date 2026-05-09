## Internal helpers (not exported)

#' @keywords internal
#' @noRd
.cytoband_cache_path <- function() {
  cache_dir <- tools::R_user_dir("ChIPBindVis", which = "cache")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  file.path(cache_dir, "cytoBand_hg38.txt.gz")
}

#' Get hg38 cytoBand as a data frame, downloading on first use
#'
#' @param cytoband_file Optional path to a user-supplied cytoBand file
#'   (\code{cytoBand.txt} or \code{cytoBand.txt.gz}). If \code{NULL}, the file
#'   is fetched from UCSC and cached to \code{tools::R_user_dir("ChIPBindVis",
#'   "cache")}.
#' @keywords internal
#' @noRd
.get_cytoband_hg38 <- function(cytoband_file = NULL) {
  if (!is.null(cytoband_file)) {
    if (!file.exists(cytoband_file)) {
      stop("cytoband_file not found: ", cytoband_file, call. = FALSE)
    }
    path <- cytoband_file
  } else {
    path <- .cytoband_cache_path()
    if (!file.exists(path)) {
      url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
      message("Downloading hg38 cytoBand from UCSC (one-time, ~7 KB)...")
      tryCatch(
        utils::download.file(url, path, mode = "wb", quiet = TRUE),
        error = function(e) {
          stop(
            "Failed to download cytoBand from UCSC.\n",
            "Please download manually from:\n  ", url, "\n",
            "and pass the path via cytoband_file = '...'\n",
            "Original error: ", conditionMessage(e),
            call. = FALSE
          )
        }
      )
    }
  }

  con <- gzfile(path)
  on.exit(close(con), add = TRUE)
  cyto <- utils::read.table(
    con, sep = "\t",
    col.names = c("chrom", "chromStart", "chromEnd", "name", "gieStain"),
    stringsAsFactors = FALSE
  )
  cyto
}

#' Read peaks and ensure UCSC-style ("chr1") seqlevels
#'
#' @keywords internal
#' @noRd
.read_peaks <- function(peaks_file) {
  peaks <- rtracklayer::import(peaks_file, format = "narrowPeak")
  if (!any(grepl("^chr", GenomeInfoDb::seqlevels(peaks)))) {
    GenomeInfoDb::seqlevelsStyle(peaks) <- "UCSC"
  }
  peaks
}
