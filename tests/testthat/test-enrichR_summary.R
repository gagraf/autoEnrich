test_that("bad threshold_sequence errors", {
  df <- data.frame(gene = c("A", "B"), value = c(1, 2))
  expect_error(
    enrichR_summary(df, "gene", "value", threshold_sequence = c(0, 1),
                     database = "WikiPathways_2024_Mouse"),
    "threshold_sequence"
  )
})

test_that("bad directions errors", {
  df <- data.frame(gene = c("A", "B"), value = c(1, 2))
  expect_error(
    enrichR_summary(df, "gene", "value", directions = "Sideways",
                     database = "WikiPathways_2024_Mouse"),
    "directions"
  )
})

test_that("too-strict min_genes errors before hitting the network", {
  df <- data.frame(gene = c("A", "B"), value = c(1, 2))
  expect_error(
    enrichR_summary(df, "gene", "value", min_genes = 100,
                     database = "WikiPathways_2024_Mouse"),
    "min_genes"
  )
})
