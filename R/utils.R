## ============================================================
##  Internal helpers (not exported)
##  Supports: hg38, hg19, mm10, mm39
## ============================================================

# ── Genome resource table ────────────────────────────────────
#' @keywords internal
#' @noRd
.genome_resources <- function(genome) {
  res <- list(
    hg38 = list(
      txdb         = "TxDb.Hsapiens.UCSC.hg38.knownGene",
      orgdb        = "org.Hs.eg.db",
      std_chr      = paste0("chr", c(1:22, "X", "Y")),
      cytoband_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
    ),
    hg19 = list(
      txdb         = "TxDb.Hsapiens.UCSC.hg19.knownGene",
      orgdb        = "org.Hs.eg.db",
      std_chr      = paste0("chr", c(1:22, "X", "Y")),
      cytoband_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz"
    ),
    mm10 = list(
      txdb         = "TxDb.Mmusculus.UCSC.mm10.knownGene",
      orgdb        = "org.Mm.eg.db",
      std_chr      = paste0("chr", c(1:19, "X", "Y")),
      cytoband_url = "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz"
    ),
    mm39 = list(
      txdb         = "TxDb.Mmusculus.UCSC.mm39.knownGene",
      orgdb        = "org.Mm.eg.db",
      std_chr      = paste0("chr", c(1:19, "X", "Y")),
      cytoband_url = "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/cytoBand.txt.gz"
    )
  )
  if (!genome %in% names(res)) {
    stop("Unsupported genome: '", genome, "'. Choose from: ",
         paste(names(res), collapse = ", "), call. = FALSE)
  }
  res[[genome]]
}

# ── Dynamic package loading ──────────────────────────────────
# TxDb/OrgDb packages export an object with the same name as the package

#' @keywords internal
#' @noRd
.load_txdb <- function(pkg_name) {
  if (!requireNamespace(pkg_name, quietly = TRUE))
    stop("Required package not installed: ", pkg_name, "\n",
         "Install with: BiocManager::install('", pkg_name, "')", call. = FALSE)
  getExportedValue(pkg_name, pkg_name)
}

#' @keywords internal
#' @noRd
.load_orgdb <- function(pkg_name) {
  if (!requireNamespace(pkg_name, quietly = TRUE))
    stop("Required package not installed: ", pkg_name, "\n",
         "Install with: BiocManager::install('", pkg_name, "')", call. = FALSE)
  getExportedValue(pkg_name, pkg_name)
}

# ── CytoBand: genome-aware, per-chromosome filtered ─────────
#' @keywords internal
#' @noRd
.cytoband_cache_path <- function(genome) {
  cache_dir <- tools::R_user_dir("ChIPBindVis", which = "cache")
  if (!dir.exists(cache_dir))
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(cache_dir, paste0("cytoBand_", genome, ".txt.gz"))
}

#' @keywords internal
#' @noRd
.get_cytoband <- function(genome, chromosome, cytoband_file = NULL) {
  res  <- .genome_resources(genome)
  path <- if (!is.null(cytoband_file)) {
    if (!file.exists(cytoband_file))
      stop("cytoband_file not found: ", cytoband_file, call. = FALSE)
    cytoband_file
  } else {
    p <- .cytoband_cache_path(genome)
    if (!file.exists(p)) {
      message("Downloading ", genome, " cytoBand from UCSC (one-time, ~7 KB)...")
      tryCatch(
        utils::download.file(res$cytoband_url, p, mode = "wb", quiet = TRUE),
        error = function(e)
          stop("Failed to download cytoBand.\nDownload manually from:\n  ",
               res$cytoband_url, "\nand pass via cytoband_file = '...'\n",
               conditionMessage(e), call. = FALSE)
      )
    }
    p
  }

  # Pass path directly — read.table handles .gz natively.
  # Previously used gzfile() + on.exit(close()), but read.table closes the
  # connection internally, causing "invalid connection" on the second close.
  cyto <- utils::read.table(
    path, sep = "\t",
    col.names        = c("chrom", "chromStart", "chromEnd", "name", "gieStain"),
    stringsAsFactors = FALSE
  )

  # Filter to single chromosome and reset rownames.
  # Without this, Gviz's IdeogramTrack calculates axis breaks across ALL
  # chromosomes in the table and can produce "breaks not unique" errors.
  cyto_chr           <- cyto[cyto$chrom == chromosome, ]
  rownames(cyto_chr) <- NULL
  cyto_chr
}

# ── Peak reading: narrowPeak + broadPeak auto-detect ─────────
#' @keywords internal
#' @noRd
.read_peaks <- function(peaks_file) {
  # Detect format from filename; fall back to narrowPeak
  fmt <- if (grepl("broad", basename(peaks_file), ignore.case = TRUE)) {
    "broadPeak"
  } else {
    "narrowPeak"
  }

  peaks <- tryCatch(
    rtracklayer::import(peaks_file, format = fmt),
    error = function(e) {
      message("Could not read as ", fmt, ", retrying as BED...")
      rtracklayer::import(peaks_file, format = "BED")
    }
  )

  # Normalise to UCSC-style chr* seqlevels
  if (!any(grepl("^chr", GenomeInfoDb::seqlevels(peaks))))
    GenomeInfoDb::seqlevelsStyle(peaks) <- "UCSC"

  peaks
}

# ── Axis label helper ────────────────────────────────────────
# 2000 → "2kb",  500 → "500bp"
#' @keywords internal
#' @noRd
.bp_label <- function(bp) {
  if (bp %% 1000 == 0) paste0(bp / 1000, "kb") else paste0(bp, "bp")
}
