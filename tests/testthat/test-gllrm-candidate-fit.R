gllrm_candidate_data <- function() {
  rows <- expand.grid(
    I1 = 0:1,
    I2 = 0:1,
    I3 = 0:1,
    I4 = 0:1,
    X1 = 0:1,
    X2 = 0:1,
    KEEP.OUT.ATTRS = FALSE
  )
  rows <- rows[rep(seq_len(nrow(rows)), 2L), , drop = FALSE]
  rows$ID <- seq_len(nrow(rows))
  rows[, c("ID", "I1", "I2", "I3", "I4", "X1", "X2")]
}

gllrm_candidate_fit <- function(max_step = 20L) {
  analysis <- gRm(
    gllrm_candidate_data(),
    items = c("I1", "I2", "I3", "I4"),
    exogenous = c("X1", "X2"),
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1, I4 = 0:1),
    exogenous_levels = list(X1 = 0:1, X2 = 0:1),
    score_cuts = c(1L, 4L)
  )
  fit(gllrm(analysis, ld = ~ I1:I2 + I2:I4, dif = ~ I2:X2), max_step = max_step)
}

test_that("GLLRM candidate spec builders preserve source term order", {
  fitted <- gllrm_candidate_fit(max_step = 10L)

  ld_spec <- gllrm_candidate_ld_spec(fitted, item1 = 1L, item2 = 3L)
  dif_spec <- gllrm_candidate_dif_spec(fitted, item = 1L, background = 1L)

  expect_equal(
    paste(ld_spec$ld$item1, ld_spec$ld$item2, sep = ":"),
    c("I1:I2", "I1:I3", "I2:I4")
  )
  expect_equal(
    paste(dif_spec$dif$item, dif_spec$dif$exogenous, sep = ":"),
    c("I1:X1", "I2:X2")
  )
})

test_that("lightweight included LD candidate fit matches full refit numerics", {
  fitted <- gllrm_candidate_fit(max_step = 20L)
  base_loglike <- fitted$fit$log_likelihood

  full <- fit_gllrm_with_added_ld_full_refit(
    fitted,
    item1 = 1L,
    item2 = 3L,
    max_step = 20L,
    max_delta = 0.0001
  )
  lightweight <- fit_gllrm_candidate_ld(
    fitted,
    item1 = 1L,
    item2 = 3L,
    max_step = 20L,
    max_delta = 0.0001
  )

  expect_s3_class(lightweight, "gRm_gllrm_candidate_fit")
  expect_true("stop_reason" %in% names(lightweight))
  expect_equal(lightweight$log_likelihood, full$fit$log_likelihood, tolerance = 0)
  expect_equal(lightweight$n_step, full$fit$n_step)
  expect_equal(lightweight$report_delta, full$fit$report_delta, tolerance = 0)
  expect_equal(lightweight$converged, full$convergence$converged)
  expect_equal(lightweight$stop_reason, full$fit$stop_reason)
  expect_equal(lightweight$context$observed_ld, full$fit$context$observed_ld)

  ld_index <- gllrm_context_ld_index(lightweight$context, 1L, 3L)
  full_ld_index <- gllrm_context_ld_index(full$fit$context, 1L, 3L)
  df <- source_ij_observed_df(lightweight$context$observed_ld[[ld_index]])
  full_df <- source_ij_observed_df(full$fit$context$observed_ld[[full_ld_index]])
  clr <- 2 * abs(base_loglike - lightweight$log_likelihood)
  full_clr <- 2 * abs(base_loglike - full$fit$log_likelihood)

  expect_equal(df, full_df)
  expect_equal(clr, full_clr, tolerance = 0)
  expect_equal(source_pfchi(df, clr), source_pfchi(full_df, full_clr), tolerance = 0)
})

test_that("lightweight included DIF candidate fit matches full refit numerics", {
  fitted <- gllrm_candidate_fit(max_step = 20L)
  base_loglike <- fitted$fit$log_likelihood

  full <- fit_gllrm_with_added_dif_full_refit(
    fitted,
    item = 1L,
    background = 1L,
    max_step = 20L,
    max_delta = 0.0001
  )
  lightweight <- fit_gllrm_candidate_dif(
    fitted,
    item = 1L,
    background = 1L,
    max_step = 20L,
    max_delta = 0.0001
  )

  expect_s3_class(lightweight, "gRm_gllrm_candidate_fit")
  expect_true("stop_reason" %in% names(lightweight))
  expect_equal(lightweight$log_likelihood, full$fit$log_likelihood, tolerance = 0)
  expect_equal(lightweight$n_step, full$fit$n_step)
  expect_equal(lightweight$report_delta, full$fit$report_delta, tolerance = 0)
  expect_equal(lightweight$converged, full$convergence$converged)
  expect_equal(lightweight$stop_reason, full$fit$stop_reason)
  expect_equal(lightweight$context$observed_dif, full$fit$context$observed_dif)

  dif_index <- gllrm_context_dif_index(lightweight$context, 1L, 1L)
  full_dif_index <- gllrm_context_dif_index(full$fit$context, 1L, 1L)
  df <- source_ix_observed_df(lightweight$context$observed_dif[[dif_index]])
  full_df <- source_ix_observed_df(full$fit$context$observed_dif[[full_dif_index]])
  clr <- 2 * abs(base_loglike - lightweight$log_likelihood)
  full_clr <- 2 * abs(base_loglike - full$fit$log_likelihood)

  expect_equal(df, full_df)
  expect_equal(clr, full_clr, tolerance = 0)
  expect_equal(source_pfchi(df, clr), source_pfchi(full_df, full_clr), tolerance = 0)
})

test_that("native GLLRM expected margins preserve LD candidate fit numerics", {
  fitted <- gllrm_candidate_fit(max_step = 10L)

  reference <- with_gllrm_expected_margin_backend(calculate_gllrm_joint_expected_margins_r, {
    fit_gllrm_candidate_ld(
      fitted,
      item1 = 1L,
      item2 = 3L,
      max_step = 25L,
      max_delta = 0.0001
    )
  })
  native <- fit_gllrm_candidate_ld(
    fitted,
    item1 = 1L,
    item2 = 3L,
    max_step = 25L,
    max_delta = 0.0001
  )

  expect_equal(native$n_step, reference$n_step)
  expect_equal(native$converged, reference$converged)
  expect_equal(native$stop_reason, reference$stop_reason)
  expect_equal(native$log_likelihood, reference$log_likelihood, tolerance = 1e-10)
  expect_equal(native$report_delta, reference$report_delta, tolerance = 1e-10)
  expect_equal(native$state$item_gamma, reference$state$item_gamma, tolerance = 1e-10)
  expect_equal(native$state$ld_parameters, reference$state$ld_parameters, tolerance = 1e-10)
  expect_equal(native$state$dif_parameters, reference$state$dif_parameters, tolerance = 1e-10)
})

test_that("native GLLRM expected margins preserve DIF candidate fit numerics", {
  fitted <- gllrm_candidate_fit(max_step = 10L)

  reference <- with_gllrm_expected_margin_backend(calculate_gllrm_joint_expected_margins_r, {
    fit_gllrm_candidate_dif(
      fitted,
      item = 1L,
      background = 1L,
      max_step = 25L,
      max_delta = 0.0001
    )
  })
  native <- fit_gllrm_candidate_dif(
    fitted,
    item = 1L,
    background = 1L,
    max_step = 25L,
    max_delta = 0.0001
  )

  expect_equal(native$n_step, reference$n_step)
  expect_equal(native$converged, reference$converged)
  expect_equal(native$stop_reason, reference$stop_reason)
  expect_equal(native$log_likelihood, reference$log_likelihood, tolerance = 1e-10)
  expect_equal(native$report_delta, reference$report_delta, tolerance = 1e-10)
  expect_equal(native$state$item_gamma, reference$state$item_gamma, tolerance = 1e-10)
  expect_equal(native$state$ld_parameters, reference$state$ld_parameters, tolerance = 1e-10)
  expect_equal(native$state$dif_parameters, reference$state$dif_parameters, tolerance = 1e-10)
})

test_that("GLLRM public LD and DIF summaries remain available after candidate refit", {
  fitted <- gllrm_candidate_fit(max_step = 20L)

  ld <- local_dependence(fitted, max_step = 20L, jobs = 1L)
  dif_result <- dif(fitted, max_step = 20L, jobs = 1L)

  ld_tests <- summary(ld)$tests
  dif_tests <- summary(dif_result)$tests

  expect_true(is.data.frame(ld_tests))
  expect_true(is.data.frame(dif_tests))
  expect_true(all(c("Chisq", "Df", "Pr(>Chisq)", "Converged", "delta") %in% names(ld_tests)))
  expect_true(all(c("Chisq", "Df", "Pr(>Chisq)", "Converged", "Stable", "delta") %in% names(dif_tests)))
  expect_false(anyNA(ld_tests$Chisq))
  expect_false(anyNA(dif_tests$Chisq))
})
