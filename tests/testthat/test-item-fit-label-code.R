item_fit_label_code_data <- function() {
  data.frame(
    ID = seq_len(16L),
    I1 = c(0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L),
    I2 = c(1L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L),
    I3 = c(0L, 0L, 1L, 1L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L),
    site = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 0L, 1L)
  )
}

test_that("base item_fit uses explicit label_code metadata", {
  old <- options(warnPartialMatchDollar = TRUE)
  on.exit(options(old), add = TRUE)

  analysis <- gRm(
    item_fit_label_code_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1),
    exogenous_levels = list(site = 0:1),
    score_cuts = c(1L, 3L)
  )
  fitted <- fit(gllrm(analysis), max_step = 50L, max_delta = 1e-4)

  expect_warning(
    values <- item_fit(fitted, include_extended = TRUE),
    NA
  )
  expect_equal(attr(values, "values", exact = TRUE)$items$item_label, c("a", "b", "c"))
})

test_that("GLLRM item_fit uses explicit label_code metadata", {
  old <- options(warnPartialMatchDollar = TRUE)
  on.exit(options(old), add = TRUE)

  analysis <- gRm(
    item_fit_label_code_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1),
    exogenous_levels = list(site = 0:1),
    score_cuts = c(1L, 3L)
  )
  fitted <- fit(gllrm(analysis, ld = ~ I1:I2, dif = ~ I3:site), max_step = 50L, max_delta = 1e-4)

  expect_warning(
    values <- item_fit(fitted, include_extended = TRUE),
    NA
  )
  expect_equal(attr(values, "values", exact = TRUE)$items$item_label, c("a", "b", "c"))
})
