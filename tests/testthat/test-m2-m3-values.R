test_that("M2/M3 Pearson statistics use source PFCHI and skip zero expected cells", {
  observed <- matrix(c(10, 5, 4, 11), nrow = 2)
  expected <- matrix(c(8, 7, 6, 9), nrow = 2)

  actual <- m2_m3_pearson_stat(observed, expected, df = 1L)

  expect_equal(actual$chi_square, sum((observed - expected)^2 / expected))
  expect_equal(actual$degrees_of_freedom, 1L)
  expect_equal(actual$p_value, source_pfchi(1L, actual$chi_square))

  zero_expected <- matrix(c(8, 0, 6, 9), nrow = 2)
  actual_zero <- m2_m3_pearson_stat(observed, zero_expected, df = 1L)
  expect_equal(
    actual_zero$chi_square,
    sum((observed[zero_expected > 0] - zero_expected[zero_expected > 0])^2 / zero_expected[zero_expected > 0])
  )
})

test_that("M2/M3 degrees of freedom match source formulas", {
  expect_equal(m2_m3_df_two_way(c(3L, 4L)), 6L)
  expect_equal(m2_m3_df_three_way(c(2L, 3L, 4L)), 17L)
})

test_that("M2/M3 aggregate helpers sum the intended source row families", {
  tests <- data.frame(
    diagnostic = c("m2_margin", "m2_margin", "m3_margin", "m3_margin"),
    margin_type = c("item_item", "item_score_group", "item_item_item", "item_exogenous_score_group"),
    background_label = c(NA, NA, NA, "f"),
    background_name = c(NA, NA, NA, "Exo"),
    chi_square = c(1, 2, 4, 8),
    degrees_of_freedom = c(3L, 5L, 7L, 11L),
    is_m2 = c(TRUE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  tests$p_value <- c(NA_real_, NA_real_, NA_real_, NA_real_)

  out <- m2_m3_aggregate_values(tests)

  expect_equal(out$m2$chi_square, 3)
  expect_equal(out$m2$degrees_of_freedom, 8L)
  expect_equal(out$m3$chi_square, 15)
  expect_equal(out$m3$degrees_of_freedom, 26L)
  expect_equal(out$item_trait$chi_square, 2)
  expect_equal(out$item_trait$degrees_of_freedom, 5L)
  expect_equal(out$invariance, data.frame(
    background_label = character(),
    background_name = character(),
    chi_square = numeric(),
    degrees_of_freedom = integer(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  ))
})
