test_that("only promoter annotations enter the heatmap gene set", {
  annotations <- data.frame(
    SYMBOL = c("TP53", "KIF18A", "MYC", "", NA, "TP53"),
    annotation = c(
      "Promoter (<=1kb)",
      "Distal Intergenic",
      "Promoter (1-2kb)",
      "Promoter (<=1kb)",
      "Promoter (<=1kb)",
      "Promoter (2-3kb)"
    ),
    stringsAsFactors = FALSE
  )

  status <- ChIPBindVis:::.promoter_target_status(
    annotations, query_symbol = "KIF18A"
  )

  expect_identical(status$target_genes, c("TP53", "MYC"))
  expect_false(status$query_promoter_peak_detected)
})

test_that("a promoter-associated query is detected without duplication", {
  annotations <- data.frame(
    SYMBOL = c("KIF18A", "KIF18A"),
    annotation = c("Promoter (<=1kb)", "Promoter (1-2kb)"),
    stringsAsFactors = FALSE
  )

  status <- ChIPBindVis:::.promoter_target_status(
    annotations, query_symbol = "KIF18A"
  )

  expect_identical(status$target_genes, "KIF18A")
  expect_true(status$query_promoter_peak_detected)
})

test_that("required ChIPseeker annotation columns are checked", {
  expect_error(
    ChIPBindVis:::.promoter_target_status(
      data.frame(SYMBOL = "KIF18A"), query_symbol = "KIF18A"
    ),
    "missing required column.*annotation"
  )
})
