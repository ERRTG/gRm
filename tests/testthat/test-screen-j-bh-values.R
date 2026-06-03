test_that("SCREEN J asymptotic partial gamma p-values match source RESULTS storage", {
  screen_j_partial_gamma <- get("screen_j_partial_gamma", envir = asNamespace("gRm"))

  saturated <- array(matrix(c(10, 0, 0, 10), nrow = 2L), dim = c(2L, 2L, 1L))
  saturated_stats <- screen_j_partial_gamma(saturated)
  expect_equal(saturated_stats$p_value, 2 * source_tail_norm(4, TRUE))

  strong <- array(matrix(c(8, 2, 1, 9), nrow = 2L), dim = c(2L, 2L, 1L))
  strong_stats <- screen_j_partial_gamma(strong)
  expect_equal(strong_stats$p_value, 0.000002)
})
