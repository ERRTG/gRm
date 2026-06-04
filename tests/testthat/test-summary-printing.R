summary_print_data <- function() {
  data.frame(
    ID = seq_len(12L),
    I1 = c(1L, 1L, 2L, 2L, 3L, 3L, 1L, 2L, 3L, 1L, 2L, 3L),
    I2 = c(1L, 2L, 1L, 3L, 2L, 3L, 3L, 1L, 2L, 2L, 3L, 1L),
    I3 = c(2L, 1L, 3L, 1L, 2L, 3L, 2L, 1L, 3L, 3L, 2L, 1L),
    site = c(1L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 2L)
  )
}

summary_print_analysis <- function() {
  gRm(
    summary_print_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID",
    score_cuts = "auto"
  )
}

expect_summary_surface <- function(object, expected_class, which = NULL) {
  out <- if (is.null(which)) summary(object) else summary(object, which = which)

  expect_s3_class(out, expected_class)
  expect_true(any(c("tables", "tests", "coefficients") %in% names(out)))
  expect_output(print(out), "DIGRAM")
  invisible(out)
}

test_that("model and fit summaries use statistical labels", {
  analysis <- summary_print_analysis()
  rasch_model <- gllrm(analysis)
  gllrm_model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I3:site)
  rasch_fit <- fit(rasch_model, max_step = 50L)
  gllrm_fit <- fit(gllrm_model, max_step = 50L)

  expect_output(print(rasch_model), "DIGRAM Rasch specification")
  expect_output(print(gllrm_model), "DIGRAM GLLRM specification")
  expect_output(print(rasch_fit), "DIGRAM Rasch fit")
  expect_output(print(gllrm_fit), "DIGRAM GLLRM fit")

  expect_summary_surface(rasch_model, "summary.gRm_model")
  expect_summary_surface(gllrm_model, "summary.gRm_model")
  expect_summary_surface(rasch_fit, "summary.gRm_fit")
  expect_summary_surface(gllrm_fit, "summary.gRm_fit")
})

test_that("summary which selects designed tables and rejects old accessor names", {
  fitted <- fit(gllrm(summary_print_analysis()), max_step = 50L)

  expect_summary_surface(item_parameters(fitted), "summary.gRm_item_parameters", "coefficients")
  expect_summary_surface(item_fit(fitted), "summary.gRm_item_fit", "tests")
  expect_summary_surface(local_dependence(fitted), "summary.gRm_local_dependence", "tests")
  expect_summary_surface(dif(fitted), "summary.gRm_dif", "tests")
  expect_summary_surface(global_homogeneity(fitted), "summary.gRm_global_homogeneity", "tests")

  expect_error(summary(fitted, which = "details"), "which")
  expect_error(summary(fitted, which = "tidy"), "which")
  expect_error(summary(fitted, which = "glance"), "which")
})

test_that("item-fit tests and items summaries expose distinct source-shaped tables", {
  fitted <- fit(gllrm(summary_print_analysis()), max_step = 50L)
  ifit <- item_fit(fitted, include_extended = TRUE)

  tests <- summary(ifit, which = "tests")$tests
  items <- summary(ifit, which = "items")$items

  expect_s3_class(summary(ifit, which = "tests"), "summary.gRm_item_fit")
  expect_s3_class(summary(ifit, which = "items"), "summary.gRm_item_fit")
  expect_true(is.data.frame(tests))
  expect_true(is.data.frame(items))
  expect_equal(nrow(items), nrow(tests))

  expect_true(all(c(
    "outfit", "outfit_sd", "p_outfit",
    "infit", "infit_sd", "p_infit",
    "observed_gamma", "expected_gamma", "gamma_sd", "p_gamma"
  ) %in% names(tests)))

  expect_true(all(c(
    "outfit_total_n", "outfit_total_observed",
    "outfit_total_expected", "outfit_total_value",
    "infit_observed", "infit_expected",
    "infit_variance", "infit_value"
  ) %in% names(items)))

  expect_false("p_gamma" %in% names(items))
  expect_false("outfit_total_value" %in% names(tests))
})

test_that("screen and score-effect summaries print BH selection context", {
  analysis <- summary_print_analysis()
  screen_result <- screen(analysis, inference = "asymptotic", jobs = 1L)
  effects <- score_effects(analysis, jobs = 1L)

  expect_summary_surface(screen_result, "summary.gRm_screen", "tests")
  expect_summary_surface(effects, "summary.gRm_score_effects", "tests")
  expect_output(print(summary(screen_result, which = "tests")), "Benjamini-Hochberg")
  expect_output(print(summary(effects, which = "tests")), "Benjamini-Hochberg")
})
