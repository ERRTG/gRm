test_that("bundle score window is derived before exogenous missingness is applied", {
  data <- data.frame(
    ID = seq_len(4L),
    I1 = c(0L, 1L, 1L, 0L),
    I2 = c(1L, 0L, 1L, 0L),
    site = c(0L, 1L, NA, 0L)
  )
  analysis <- gRm(
    data,
    items = c("I1", "I2"),
    exogenous = "site",
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1),
    exogenous_levels = list(site = 0:1),
    score_cuts = c(1L, 2L)
  )

  bundle <- build_item_parameters_bundle(analysis$project)

  expect_equal(bundle$model$largest_score, 2L)
  expect_equal(bundle$manifest$ncomplete_items, 4L)
  expect_equal(bundle$manifest$ncomplete_backgrounds, 3L)
  expect_equal(bundle$manifest$ncomplete_item_backgrounds, 3L)
  expect_equal(bundle$manifest$nvalid, 2L)
  expect_equal(bundle$manifest$nmissing_backgrounds, 1L)
})
