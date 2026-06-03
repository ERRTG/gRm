output_accessor_data <- function() {
  data.frame(
    ID = seq_len(7L),
    I1 = c(1L, 1L, 2L, 3L, 1L, 0L, 2L),
    I2 = c(1L, 2L, 2L, 3L, 3L, 1L, NA_integer_),
    I3 = c(1L, 3L, 2L, 3L, 1L, 2L, 2L),
    group = c(1L, NA_integer_, 2L, 1L, 2L, 1L, NA_integer_)
  )
}

test_that("summary is the public output accessor for analysis objects", {
  ia <- gRm(
    output_accessor_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "group",
    id = "ID"
  )

  out <- summary(ia)

  expect_s3_class(out, "summary.gRm_analysis")
  expect_true("data" %in% names(out$tables))
  expect_equal(out$data$n_items, 3L)
  expect_false("tidy" %in% getNamespaceExports("gRm"))
  expect_false("glance" %in% getNamespaceExports("gRm"))
  expect_false("details" %in% getNamespaceExports("gRm"))
})

test_that("summary which exposes model, fit, screen, and diagnostic sections", {
  ia <- gRm(
    output_accessor_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "group",
    id = "ID"
  )
  model <- gllrm(ia)
  fit_obj <- fit(model, max_step = 50L)
  screen_obj <- screen(ia, inference = "asymptotic")

  expect_true(is.data.frame(summary(model, which = "model")$model))
  expect_true(is.data.frame(summary(fit_obj, which = "fit")$fit))
  expect_true(is.data.frame(summary(screen_obj, which = "tests")$tests))
  expect_true(is.data.frame(summary(item_parameters(fit_obj), which = "coefficients")$coefficients))
  expect_true(is.data.frame(summary(item_fit(fit_obj), which = "tests")$tests))
})

test_that("internal value details remain available for source-faithful validation helpers", {
  ia <- gRm(
    output_accessor_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "group",
    id = "ID"
  )
  values <- item_fits_values(ia$project, include_extended = TRUE)

  expect_true(all(c(
    "statistics", "compact", "bh_thresholds",
    "score_n", "score_level_fit", "item_fit_summaries",
    "restscore", "local_restscore", "local_gamma",
    "local_gamma_bh_thresholds"
  ) %in% detail_names(values)))
})
