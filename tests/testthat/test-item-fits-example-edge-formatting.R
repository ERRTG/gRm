test_that("extended restscore table rows are built without per-cell overhead", {
  items <- data.frame(
    label_code = "a",
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

test_that("extended zero-variance rows preserve DIGRAM runtime carry-forward", {
  items <- data.frame(
    label_code = "h",
    name = "D3",
    raw_max = 5L,
    stringsAsFactors = FALSE
  )
  conditional <- list(vector("list", 32L))
  conditional[[1L]][[32L]] <- list(
    probabilities = c(0, 0, 0, 0, 1),
    mean = 4,
    variance = 0
  )
  previous <- as.data.frame(as.list(stats::setNames(
    c(-13.173, -9.797, -6.422, -3.047, 0.328),
    paste0("standardized_", 0:4)
  )))

  row <- extended_item_fit_score_row(
    items = items,
    data = data.frame(D3 = 4L),
    item_index = 1L,
    score = 31L,
    n_score = 1L,
    observed_rows = TRUE,
    item_scores = 31L,
    conditional = conditional,
    previous_row = previous
  )

  expect_equal(unlist(row[paste0("standardized_", 0:4)]), unlist(previous))
  expect_true(all(unlist(row[paste0("residual_", 0:4)]) == 0))
  expect_true(is.nan(row$squared_standardized_observed_average))
  expect_true(is.nan(row$infit_ratio))
  expect_identical(row$outfit_standard_error, 0)
})

test_that("compact outfit denominator includes zero-variance score counts", {
  bundle <- list(
    model = list(
      items = data.frame(name = "I1", raw_max = 5L, stringsAsFactors = FALSE),
      backgrounds = data.frame(),
      least_score = 1L,
      largest_score = 3L,
      max_total_score = 4L
    ),
    data = data.frame(I1 = 3L, score = 3L, status = 1L)
  )
  conditional <- list(vector("list", 5L))
  conditional[[1L]][[4L]] <- list(
    probabilities = c(0, 0, 0, 1, 0),
    mean = 3,
    variance = 0
  )

  accumulators <- conditional_item_fit_complete_accumulators(bundle, conditional)

  expect_identical(accumulators$n_used, 1L)
  expect_identical(accumulators$outfit_sum, 0)
  expect_identical(accumulators$infit_weight, 0)
})
