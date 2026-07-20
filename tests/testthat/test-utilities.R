test_that("axis labels use base pairs or kilobases", {
  expect_identical(ChIPBindVis:::.bp_label(2000), "2kb")
  expect_identical(ChIPBindVis:::.bp_label(500), "500bp")
})

test_that("unsupported genomes fail clearly", {
  expect_error(
    ChIPBindVis:::.genome_resources("hg18"),
    "Unsupported genome"
  )
})

test_that("the default palette exposes all plotting roles", {
  expect_named(
    morandi_palette(),
    c(
      "signal", "signal_fill", "peak", "gene", "tss", "heatmap_low",
      "heatmap_mid", "heatmap_high", "axis", "bg"
    )
  )
})
