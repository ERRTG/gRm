test_that("extended restscore table rows are built without per-cell overhead", {
  items <- data.frame(
    label = "a",
    name = "ItemA",
    stringsAsFactors = FALSE
  )
  observed <- matrix(seq_len(240), nrow = 3L)
  expected <- observed + 0.5

  timing <- system.time({
    for (i in seq_len(20L)) {
      rows <- restscore_table_rows(items, 1L, 0L, observed, expected)
    }
  })

  expect_lt(unname(timing[["elapsed"]]), 1)
  expect_equal(nrow(rows), length(observed) * 2L)
  expect_equal(rows$restscore_table, rep(c("observed", "expected"), each = length(observed)))
  expect_equal(rows$item_score[seq_len(6L)], c(0L, 0L, 0L, 0L, 0L, 0L))
  expect_equal(rows$restscore[seq_len(6L)], 0:5)
})
