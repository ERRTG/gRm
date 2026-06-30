test_that("base item_fit uses the supplied fitted Rasch state", {
  data <- data.frame(
    ID = seq_len(18L),
    I1 = c(0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L),
    I2 = c(0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 2L, 0L, 1L, 1L, 2L, 0L, 0L, 2L, 1L),
    I3 = c(2L, 1L, 0L, 2L, 1L, 0L, 2L, 1L, 0L, 0L, 2L, 1L, 1L, 0L, 2L, 2L, 0L, 1L)
  )
  analysis <- gRm(
    data,
    items = c("I1", "I2", "I3"),
    id = "ID",
    item_levels = list(I1 = 0:2, I2 = 0:2, I3 = 0:2),
    score_cuts = c(1L, 6L)
  )
  supplied <- fit(gllrm(analysis), max_step = 1L, max_delta = 1e-12)

  values <- attr(item_fit(supplied, include_extended = FALSE), "values", exact = TRUE)

  expect_equal(supplied$fit$n_step, 1L)
  expect_false(isTRUE(supplied$fit$converged))
  expect_equal(values$fit$n_step, supplied$fit$n_step)
  expect_equal(values$fit$converged, supplied$fit$converged)
  expect_equal(values$fit$item_gamma, supplied$fit$item_gamma)
})
