m2_m3_bootstrap_base_data <- function() {
  configurations <- expand.grid(
    I3 = 0:1,
    I2 = 0:1,
    I1 = 0:1,
    KEEP.OUT.ATTRS = FALSE
  )[, c("I1", "I2", "I3"), drop = FALSE]
  counts <- c(1L, 2L, 1L, 1L, 1L, 2L, 1L, 3L)
  rows <- configurations[rep(seq_len(nrow(configurations)), counts), , drop = FALSE]
  data.frame(ID = seq_len(nrow(rows)), rows, row.names = NULL)
}

m2_m3_bootstrap_base_fit <- function() {
  analysis <- gRm(
    m2_m3_bootstrap_base_data(),
    items = c("I1", "I2", "I3"),
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1),
    score_cuts = c(1L, 3L)
  )
  fit(gllrm(analysis), max_step = 100L, max_delta = 0.0001)
}

m2_m3_bootstrap_ld_dif_fit <- function() {
  rows <- expand.grid(
    I1 = 0:1,
    I2 = 0:1,
    I3 = 0:1,
    X1 = 1:2,
    X2 = 1:2,
    replicate = 1:2,
    KEEP.OUT.ATTRS = FALSE
  )
  data <- data.frame(
    ID = seq_len(nrow(rows)),
    rows[, c("I1", "I2", "I3", "X1", "X2"), drop = FALSE]
  )
  analysis <- gRm(
    data,
    items = c("I1", "I2", "I3"),
    exogenous = c("X1", "X2"),
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1),
    exogenous_levels = list(X1 = 1:2, X2 = 1:2),
    score_cuts = c(1L, 3L)
  )
  fit(
    gllrm(analysis, ld = ~ I1:I2, dif = ~ I3:X1),
    max_step = 100L,
    max_delta = 0.0001
  )
}

test_that("M2/M3 bootstrap controls are public and strict", {
  expected <- c(
    "bootstrap", "nsim", "seed", "reestimate", "bootstrap_max_step",
    "keep_bootstrap_samples", "resample_score_distribution"
  )
  expect_true(all(expected %in% names(formals(gRm::m2))))
  expect_true(all(expected %in% names(formals(gRm::m3))))

  fitted <- m2_m3_bootstrap_base_fit()
  expect_error(m3(fitted, bootstrap = NA), "`bootstrap`")
  expect_error(m3(fitted, nsim = 0L), "`nsim`")
  expect_error(m3(fitted, seed = 0L), "`seed`")
  expect_error(m3(fitted, seed = 2147483647), "`seed`")
  expect_error(m3(fitted, reestimate = 1L), "`reestimate`")
  expect_error(m3(fitted, bootstrap_max_step = 0L), "`max_step`")
  expect_error(m3(fitted, keep_bootstrap_samples = NA), "`keep_bootstrap_samples`")
  expect_error(
    m3(fitted, resample_score_distribution = NA),
    "`resample_score_distribution`"
  )
})

test_that("ordinary M2/M3 diagnostics do not run bootstrap work", {
  fitted <- m2_m3_bootstrap_base_fit()
  m2_result <- m2(fitted)
  m3_result <- m3(fitted)

  for (result in list(m2_result, m3_result)) {
    expect_identical(result$values$bootstrap$enabled, FALSE)
    expect_identical(result$values$bootstrap$source_status, "not_requested")
    expect_identical(result$values$bootstrap$nsim, 0L)
    expect_identical(result$values$bootstrap$nused, 0L)
  }
  expect_true(any(grepl(
    "Parametric bootstrap: not requested",
    capture.output(print(m3_result)),
    fixed = TRUE
  )))
})

test_that("seeded M3 bootstrap is private, reproducible, and score-conditional", {
  fitted <- m2_m3_bootstrap_base_fit()
  set.seed(104L)
  global_state <- .Random.seed
  first <- m3(
    fitted,
    bootstrap = TRUE,
    nsim = 3L,
    seed = 47L,
    reestimate = FALSE,
    keep_bootstrap_samples = TRUE
  )
  expect_identical(.Random.seed, global_state)
  second <- m3(
    fitted,
    bootstrap = TRUE,
    nsim = 3L,
    seed = 47L,
    reestimate = FALSE,
    keep_bootstrap_samples = TRUE
  )
  expect_identical(.Random.seed, global_state)

  first_bootstrap <- first$values$bootstrap
  second_bootstrap <- second$values$bootstrap
  expect_equal(first_bootstrap$replicates, second_bootstrap$replicates)
  expect_equal(first_bootstrap$aggregate_replicates, second_bootstrap$aggregate_replicates)
  expect_equal(first_bootstrap$margin_replicates, second_bootstrap$margin_replicates)
  expect_equal(first_bootstrap$samples, second_bootstrap$samples)
  expect_identical(first_bootstrap$seed, 47L)
  expect_identical(first_bootstrap$nsim, 3L)
  expect_identical(first_bootstrap$nused, 3L)
  expect_identical(first_bootstrap$rng_draws, 24L)
  expect_identical(length(first_bootstrap$samples), 3L)
  expect_identical(
    first_bootstrap$generation_mode,
    "bounded_joint_harness_equivalent"
  )

  observed_scores <- as.integer(table(factor(rowSums(
    m2_m3_bootstrap_base_data()[, c("I1", "I2", "I3")]
  ), levels = 0:3)))
  for (sample in first_bootstrap$samples) {
    expect_identical(
      as.integer(table(factor(sample$score, levels = 0:3))),
      observed_scores
    )
    expect_true(all(rowSums(sample$items) == sample$score))
  }

  printed <- capture.output(print(first))
  summarized <- summary(first)
  expect_true(any(grepl("Bootstrap Pr", printed, fixed = TRUE)))
  expect_true(any(grepl("3/3 samples accepted", printed, fixed = TRUE)))
  expect_true("bootstrap" %in% names(summarized$tables))
  expect_true(all(c("Bootstrap Pr", "Bootstrap n") %in% names(summarized$tables$aggregates)))
})

test_that("the preserved score-distribution resampling control is an explicit no-op", {
  fitted <- m2_m3_bootstrap_base_fit()
  fixed <- m3(
    fitted,
    bootstrap = TRUE,
    nsim = 2L,
    seed = 59L,
    reestimate = FALSE,
    keep_bootstrap_samples = TRUE,
    resample_score_distribution = FALSE
  )$values$bootstrap
  requested <- m3(
    fitted,
    bootstrap = TRUE,
    nsim = 2L,
    seed = 59L,
    reestimate = FALSE,
    keep_bootstrap_samples = TRUE,
    resample_score_distribution = TRUE
  )$values$bootstrap

  expect_equal(requested$samples, fixed$samples)
  expect_equal(requested$replicates, fixed$replicates)
  expect_equal(requested$aggregate_replicates, fixed$aggregate_replicates)
  expect_equal(requested$margin_replicates, fixed$margin_replicates)
  expect_identical(requested$resample_score_distribution, TRUE)
  expect_identical(fixed$resample_score_distribution, FALSE)
  expect_identical(requested$score_distribution_mode, "fixed_observed_source_noop_branch")
})

test_that("bootstrap generation uses active LD and DIF and preserves source extreme rows", {
  fitted <- m2_m3_bootstrap_ld_dif_fit()
  result <- m3(
    fitted,
    bootstrap = TRUE,
    nsim = 2L,
    seed = 83L,
    reestimate = FALSE,
    keep_bootstrap_samples = TRUE
  )
  bootstrap <- result$values$bootstrap

  expect_true(bootstrap$possible)
  expect_identical(bootstrap$capability$observed$largest_component, 2L)
  expect_equal(bootstrap$capability$observed$dif_probability_strata, 2)
  expect_equal(bootstrap$capability$observed$exogenous_strata, 4)
  expect_identical(bootstrap$capability$bounds$maximum_items, 40L)
  expect_identical(bootstrap$capability$bounds$maximum_component_items, 4L)
  expect_identical(bootstrap$capability$bounds$maximum_item_score, 7L)
  expect_identical(bootstrap$capability$bounds$maximum_dif_probability_strata, 64L)
  expect_identical(bootstrap$capability$bounds$maximum_exogenous_strata, 216L)

  for (sample in bootstrap$samples) {
    extreme <- sample$score %in% c(0L, 3L)
    expect_true(all(sample$backgrounds[extreme, "X2"] == 1L))
    expect_identical(sort(unique(sample$backgrounds[extreme, "X1"])), 1:2)
    expect_identical(sort(unique(sample$backgrounds[!extreme, "X2"])), 1:2)
    expect_identical(
      as.integer(table(factor(sample$score, levels = 0:3))),
      c(8L, 24L, 24L, 8L)
    )
  }
})

test_that("bootstrap refits record source acceptance and convergence state", {
  fitted <- m2_m3_bootstrap_base_fit()
  result <- m3(
    fitted,
    bootstrap = TRUE,
    nsim = 1L,
    seed = 67L,
    reestimate = TRUE,
    bootstrap_max_step = 50L
  )
  bootstrap <- result$values$bootstrap
  replicate <- bootstrap$replicates[1L, , drop = FALSE]

  expect_identical(bootstrap$reestimate, TRUE)
  expect_identical(bootstrap$max_step, 50L)
  expect_identical(bootstrap$acceptance_delta, 0.1)
  expect_true(replicate$reestimated)
  expect_false(is.na(replicate$iterations))
  expect_identical(replicate$accepted, replicate$final_delta < 0.1)
  expect_identical(bootstrap$nused, as.integer(replicate$accepted))
  expect_identical(bootstrap$samples, list())
})

test_that("M2 exposes calibrated aggregate and margin tables without inventing M3", {
  fitted <- m2_m3_bootstrap_base_fit()
  result <- m2(
    fitted,
    bootstrap = TRUE,
    nsim = 1L,
    seed = 71L,
    reestimate = FALSE
  )

  expect_s3_class(result, "gRm_m2")
  expect_true("bootstrap_p_value" %in% names(result$values$aggregate))
  expect_true("bootstrap_p_value" %in% names(result$values$item_trait))
  expect_true("bootstrap_p_value" %in% names(result$values$tests))
  expect_false("M3" %in% result$values$bootstrap$aggregate_summary$diagnostic)
})

test_that("bootstrap exceedance summaries use accepted samples and source p direction", {
  observed <- data.frame(
    diagnostic = c("M2", "M3"),
    p_value = c(0.2, 0.5),
    stringsAsFactors = FALSE
  )
  simulated <- data.frame(
    diagnostic = rep(c("M2", "M3"), each = 3L),
    accepted = rep(c(TRUE, TRUE, FALSE), 2L),
    p_value = c(0.2, 0.3, 0, 0.1, 0.5, 0),
    stringsAsFactors = FALSE
  )

  summary <- m2_m3_bootstrap_summarize_rows(
    observed,
    simulated,
    key = "diagnostic",
    nused = 2L
  )

  expect_identical(summary$bootstrap_extreme_count, c(1L, 2L))
  expect_identical(summary$bootstrap_nused, c(2L, 2L))
  expect_equal(summary$bootstrap_p_value, c(0.5, 1))
})
