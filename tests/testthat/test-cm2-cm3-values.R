test_that("CM2/CM3 Pearson statistics use source PFCHI and skip zero expected cells", {
  observed <- matrix(c(10, 5, 4, 11), nrow = 2)
  expected <- matrix(c(8, 7, 6, 9), nrow = 2)

  actual <- cm2_cm3_pearson_stat(observed, expected, df = 1L)

  expect_equal(actual$chi_square, sum((observed - expected)^2 / expected))
  expect_equal(actual$degrees_of_freedom, 1L)
  expect_equal(actual$p_value, source_pfchi(1L, actual$chi_square))

  zero_expected <- matrix(c(8, 0, 6, 9), nrow = 2)
  actual_zero <- cm2_cm3_pearson_stat(observed, zero_expected, df = 1L)
  expect_equal(
    actual_zero$chi_square,
    sum((observed[zero_expected > 0] - zero_expected[zero_expected > 0])^2 / zero_expected[zero_expected > 0])
  )
})

test_that("CM2/CM3 degrees of freedom match source formulas", {
  expect_equal(cm2_cm3_df_two_way(c(3L, 4L)), 6L)
  expect_equal(cm2_cm3_df_three_way(c(2L, 3L, 4L)), 17L)
})

test_that("CM2/CM3 aggregate helpers sum the intended source row families", {
  tests <- data.frame(
    diagnostic = c("cm2_margin", "cm2_margin", "cm3_margin", "cm3_margin"),
    margin_type = c("item_item", "item_score_group", "item_item_item", "item_exogenous_score_group"),
    background_label = c(NA, NA, NA, "f"),
    background_name = c(NA, NA, NA, "Exo"),
    chi_square = c(1, 2, 4, 8),
    degrees_of_freedom = c(3L, 5L, 7L, 11L),
    is_cm2 = c(TRUE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  tests$p_value <- c(NA_real_, NA_real_, NA_real_, NA_real_)

  out <- cm2_cm3_aggregate_values(tests)

  expect_equal(out$cm2$chi_square, 3)
  expect_equal(out$cm2$degrees_of_freedom, 8L)
  expect_equal(out$cm3$chi_square, 15)
  expect_equal(out$cm3$degrees_of_freedom, 26L)
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

test_that("CM2/CM3 invariance retains every fitted exogenous variable at zero df", {
  context <- list(project = list(
    items = data.frame(
      name = c("I1", "I2"),
      label_code = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    backgrounds = data.frame(
      name = c("X1", "X2"),
      label_code = c("X", "Y"),
      stringsAsFactors = FALSE
    )
  ))
  tests <- data.frame(
    diagnostic = "I1:X1",
    margin_type = "item_exogenous",
    background_label = "X",
    background_name = "X1",
    chi_square = 4,
    degrees_of_freedom = 2L,
    p_value = source_pfchi(2L, 4),
    is_cm2 = TRUE,
    stringsAsFactors = FALSE
  )

  invariance <- cm2_cm3_aggregate_values(tests, context)$invariance
  expect_equal(invariance$background_label, c("X", "Y"))
  expect_equal(invariance$background_name, c("X1", "X2"))
  expect_equal(invariance$chi_square, c(4, 0))
  expect_equal(invariance$degrees_of_freedom, c(2L, 0L))
  expect_equal(invariance$p_value, c(source_pfchi(2L, 4), 0))

  empty <- cm2_cm3_aggregate_values(data.frame(), context)$invariance
  expect_equal(empty$chi_square, c(0, 0))
  expect_equal(empty$degrees_of_freedom, c(0L, 0L))
  expect_equal(empty$p_value, c(0, 0))
})
