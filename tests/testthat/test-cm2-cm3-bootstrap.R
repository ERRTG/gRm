cm2_cm3_bootstrap_base_data <- function() {
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

cm2_cm3_bootstrap_base_fit <- function() {
  analysis <- gRm(
    cm2_cm3_bootstrap_base_data(),
    items = c("I1", "I2", "I3"),
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1),
    score_cuts = c(1L, 3L)
  )
  fit(gllrm(analysis), max_step = 100L, max_delta = 0.0001)
}

cm2_cm3_bootstrap_ld_dif_fit <- function() {
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

cm2_cm3_bootstrap_zero_invariance_fit <- function() {
  rows <- expand.grid(
    I1 = 0:1,
    I2 = 0:1,
    I3 = 0:1,
    X1 = 1:2,
    replicate = 1:2,
    KEEP.OUT.ATTRS = FALSE
  )
  data <- data.frame(ID = seq_len(nrow(rows)), rows)
  analysis <- gRm(
    data,
    items = c("I1", "I2", "I3"),
    exogenous = "X1",
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1),
    exogenous_levels = list(X1 = 1:2),
    score_cuts = c(1L, 3L)
  )
  fit(
    gllrm(analysis, dif = ~ I1:X1 + I2:X1),
    max_step = 100L,
    max_delta = 0.0001
  )
}

test_that("CM2/CM3 bootstrap controls are public and strict", {
  expected <- c(
    "bootstrap", "nsim", "seed", "reestimate", "bootstrap_max_step",
    "bootstrap_jobs", "keep_bootstrap_samples", "resample_score_distribution"
  )
  expect_true(all(expected %in% names(formals(gRm::cm2))))
  expect_true(all(expected %in% names(formals(gRm::cm3))))
  expect_identical(formals(gRm::cm2)$bootstrap_max_step, 5000L)
  expect_identical(formals(gRm::cm3)$bootstrap_max_step, 5000L)
  expect_identical(formals(gRm::cm2)$bootstrap_jobs, 1L)
  expect_identical(formals(gRm::cm3)$bootstrap_jobs, 1L)

  fitted <- cm2_cm3_bootstrap_base_fit()
  expect_error(cm3(fitted, bootstrap = NA), "`bootstrap`")
  expect_error(cm3(fitted, nsim = 0L), "`nsim`")
  expect_error(cm3(fitted, seed = -1), "`seed`")
  expect_error(cm3(fitted, seed = 4294967296), "`seed`")
  expect_error(cm3(fitted, seed = 1.5), "`seed`")
  expect_identical(
    normalize_cm2_cm3_bootstrap_control(
      TRUE, 1L, 0, TRUE, 1L, 1L, FALSE, FALSE
    )$seed,
    0
  )
  expect_identical(
    normalize_cm2_cm3_bootstrap_control(
      TRUE, 1L, 4294967295, TRUE, 1L, 1L, FALSE, FALSE
    )$seed,
    4294967295
  )
  expect_error(cm3(fitted, reestimate = 1L), "`reestimate`")
  expect_error(cm3(fitted, bootstrap_max_step = 0L), "`max_step`")
  expect_error(cm3(fitted, bootstrap_jobs = 0L), "`jobs`")
  expect_error(cm3(fitted, bootstrap_jobs = 1.5), "`jobs`")
  expect_error(cm3(fitted, keep_bootstrap_samples = NA), "`keep_bootstrap_samples`")
  expect_error(
    cm3(fitted, resample_score_distribution = NA),
    "`resample_score_distribution`"
  )
})

test_that("the private bootstrap stream is the exact Delphi 4 uint32 trajectory", {
  rng <- new_cm2_cm3_bootstrap_rng(47)
  expected_states <- c(
    2039495916, 110813341, 3921630994, 2184313691, 57464008
  )
  observed_uniforms <- vapply(expected_states, function(unused) rng$uniform(), numeric(1L))
  expect_identical(observed_uniforms, expected_states / 4294967296)
  expect_identical(rng$state(), expected_states[[length(expected_states)]])
  expect_identical(rng$draws(), 5L)

  # This long trajectory was recovered independently from the Delphi 4 LCG.
  # It exercises enough wrapped multiplications to catch any binary64 product
  # implementation that loses uint32 low bits.
  long_rng <- new_cm2_cm3_bootstrap_rng(61821250)
  for (draw in seq_len(2331644L)) {
    long_rng$uniform()
  }
  expect_identical(long_rng$state(), 3180639750)
  expect_identical(long_rng$draws(), 2331644L)

  set.seed(908L)
  global_state <- .Random.seed
  randomize_rng <- new_cm2_cm3_bootstrap_rng()
  expect_identical(.Random.seed, global_state)
  expect_gte(randomize_rng$seed(), 0)
  expect_lt(randomize_rng$seed(), 86400000)
})

test_that("ordinary CM2/CM3 diagnostics do not run bootstrap work", {
  fitted <- cm2_cm3_bootstrap_base_fit()
  cm2_result <- cm2(fitted)
  cm3_result <- cm3(fitted)

  for (result in list(cm2_result, cm3_result)) {
    expect_identical(result$values$bootstrap$enabled, FALSE)
    expect_identical(result$values$bootstrap$source_status, "not_requested")
    expect_identical(result$values$bootstrap$nsim, 0L)
    expect_identical(result$values$bootstrap$nused, 0L)
  }
  expect_true(any(grepl(
    "Parametric bootstrap: not requested",
    capture.output(print(cm3_result)),
    fixed = TRUE
  )))
})

test_that("seeded CM3 bootstrap is private, reproducible, and score-conditional", {
  fitted <- cm2_cm3_bootstrap_base_fit()
  set.seed(104L)
  global_state <- .Random.seed
  first <- cm3(
    fitted,
    bootstrap = TRUE,
    nsim = 3L,
    seed = 47L,
    reestimate = FALSE,
    keep_bootstrap_samples = TRUE
  )
  expect_identical(.Random.seed, global_state)
  second <- cm3(
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
  expect_identical(first_bootstrap$seed, 47)
  expect_identical(first_bootstrap$nsim, 3L)
  expect_identical(first_bootstrap$nused, 3L)
  expect_identical(first_bootstrap$rng_draws, 48L)
  expect_identical(first_bootstrap$final_rng_state, 1214258367)
  expect_identical(
    first_bootstrap$replicates$rng_start_state,
    c(47, 2666053727, 993162383)
  )
  expect_identical(
    first_bootstrap$replicates$rng_final_state,
    c(2666053727, 993162383, 1214258367)
  )
  expect_identical(first_bootstrap$replicates$rng_draws, rep(16L, 3L))
  expect_identical(length(first_bootstrap$samples), 3L)
  expect_identical(
    first_bootstrap$generation_mode,
    "source_component_conditional"
  )
  expect_identical(first_bootstrap$rng, "delphi4_lcg_134775813_uint32")

  expected_items <- list(
    matrix(c(
      0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L,
      0L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 0L, 1L, 0L, 1L,
      0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L
    ), ncol = 3L, byrow = TRUE),
    matrix(c(
      1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L,
      1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L, 1L,
      0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L
    ), ncol = 3L, byrow = TRUE),
    matrix(c(
      1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L,
      1L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 1L, 0L, 1L, 1L,
      0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L
    ), ncol = 3L, byrow = TRUE)
  )
  dimnames(expected_items[[1L]]) <- dimnames(first_bootstrap$samples[[1L]]$items)
  dimnames(expected_items[[2L]]) <- dimnames(first_bootstrap$samples[[2L]]$items)
  dimnames(expected_items[[3L]]) <- dimnames(first_bootstrap$samples[[3L]]$items)
  expect_identical(lapply(first_bootstrap$samples, `[[`, "items"), expected_items)

  refit_bundle <- cm2_cm3_bootstrap_bundle(
    fitted,
    first_bootstrap$samples[[1L]]
  )
  refit_counts <- rasch_counts(refit_bundle)
  expect_identical(refit_bundle$model$least_score, 0L)
  expect_identical(refit_bundle$model$largest_score, 3L)
  expect_identical(refit_bundle$manifest$nvalid, nrow(refit_bundle$data))
  expect_identical(refit_counts$n_valid, nrow(refit_bundle$data))
  expect_identical(
    unname(refit_counts$score_counts[c("0", "3")]),
    c(1L, 3L)
  )
  for (item_index in seq_len(nrow(refit_counts$item_counts))) {
    expect_identical(
      unname(refit_counts$item_counts[item_index, ]),
      as.integer(table(factor(
        first_bootstrap$samples[[1L]]$items[, item_index],
        levels = 0:1
      )))
    )
    expect_identical(sum(refit_counts$item_counts[item_index, ]), 12L)
  }

  observed_scores <- as.integer(table(factor(rowSums(
    cm2_cm3_bootstrap_base_data()[, c("I1", "I2", "I3")]
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

test_that("fixed-seed bootstrap generation and fit state do not depend on focal selection", {
  fitted <- cm2_cm3_bootstrap_base_fit()
  arguments <- list(
    fit = fitted,
    bootstrap = TRUE,
    nsim = 2L,
    seed = 47L,
    reestimate = FALSE,
    keep_bootstrap_samples = TRUE
  )
  all_items <- do.call(cm3, arguments)$values
  subset <- do.call(cm3, c(arguments, list(items = c("I2", "I1"))))$values

  expect_identical(all_items$bootstrap$samples, subset$bootstrap$samples)
  expect_identical(all_items$bootstrap$replicates, subset$bootstrap$replicates)
  expect_identical(all_items$bootstrap$seed, subset$bootstrap$seed)
  expect_identical(all_items$bootstrap$final_rng_state, subset$bootstrap$final_rng_state)
  expect_identical(all_items$bootstrap$rng_draws, subset$bootstrap$rng_draws)
  expect_identical(all_items$bootstrap$nused, subset$bootstrap$nused)
  expect_false(identical(all_items$margin_specs, subset$margin_specs))
  expect_lt(nrow(subset$tests), nrow(all_items$tests))

  # Canonical command 197 re-estimates the complete current GLLRM. One bounded
  # refit replicate proves that focal selection changes only the evaluated
  # margin list, not the generated sample or full-model fit trajectory.
  refit_arguments <- list(
    fit = fitted,
    bootstrap = TRUE,
    nsim = 1L,
    seed = 71L,
    reestimate = TRUE,
    bootstrap_max_step = 50L,
    keep_bootstrap_samples = TRUE
  )
  all_items_refit <- do.call(cm3, refit_arguments)$values
  subset_refit <- do.call(
    cm3,
    c(refit_arguments, list(items = c("I2", "I1")))
  )$values

  expect_identical(all_items_refit$bootstrap$samples, subset_refit$bootstrap$samples)
  expect_identical(all_items_refit$bootstrap$replicates, subset_refit$bootstrap$replicates)
  expect_identical(all_items_refit$bootstrap$seed, subset_refit$bootstrap$seed)
  expect_identical(
    all_items_refit$bootstrap$final_rng_state,
    subset_refit$bootstrap$final_rng_state
  )
  expect_identical(all_items_refit$bootstrap$rng_draws, subset_refit$bootstrap$rng_draws)
  expect_identical(all_items_refit$bootstrap$nused, subset_refit$bootstrap$nused)
  expect_false(identical(all_items_refit$margin_specs, subset_refit$margin_specs))

  # The independently evaluated selection may change membership, aggregate
  # statistics, BH inputs, and margin indices only. Every ordinary and
  # bootstrap statistic for a margin shared with the all-item run must be
  # byte-for-byte the same after source-order position is removed.
  expect_shared_margin_rows <- function(full, selected, key) {
    encode <- function(frame) {
      do.call(paste, c(frame[key], sep = "\r"))
    }
    index <- match(encode(selected), encode(full))
    expect_false(anyNA(index))
    columns <- setdiff(intersect(names(selected), names(full)), "margin_index")
    selected <- selected[, columns, drop = FALSE]
    full <- full[index, columns, drop = FALSE]
    rownames(selected) <- NULL
    rownames(full) <- NULL
    expect_identical(selected, full)
  }
  expect_shared_margin_rows(
    all_items_refit$tests,
    subset_refit$tests,
    "margin"
  )
  expect_shared_margin_rows(
    all_items_refit$bootstrap$margin_replicates,
    subset_refit$bootstrap$margin_replicates,
    c("margin", "replicate")
  )
  expect_shared_margin_rows(
    all_items_refit$bootstrap$margin_summary,
    subset_refit$bootstrap$margin_summary,
    "margin"
  )
})

test_that("CM2 bootstrap is the two-way projection of CM3 for one selection", {
  fitted <- cm2_cm3_bootstrap_base_fit()
  arguments <- list(
    fit = fitted,
    items = c(2L, 1L),
    bootstrap = TRUE,
    nsim = 1L,
    seed = 71L,
    reestimate = FALSE,
    keep_bootstrap_samples = TRUE
  )
  cm2_result <- do.call(cm2, arguments)$values
  cm3_result <- do.call(cm3, arguments)$values

  expect_identical(cm2_result$bootstrap$samples, cm3_result$bootstrap$samples)
  expect_identical(cm2_result$bootstrap$replicates, cm3_result$bootstrap$replicates)
  expect_equal(cm2_result$aggregate, cm3_result$cm2)
  expect_equal(cm2_result$item_trait, cm3_result$item_trait)
  expect_equal(cm2_result$invariance, cm3_result$invariance)
  expect_equal(
    cm2_result$tests,
    cm3_result$tests[cm3_result$tests$is_cm2, , drop = FALSE]
  )
})

test_that("zero-df invariance is preserved through ordinary and bootstrap paths", {
  fitted <- cm2_cm3_bootstrap_zero_invariance_fit()
  result <- cm3(
    fitted,
    items = c("I2", "I1"),
    bootstrap = TRUE,
    nsim = 1L,
    seed = 9L,
    reestimate = FALSE
  )
  invariance <- result$values$invariance

  expect_equal(invariance$background_name, "X1")
  expect_equal(invariance$chi_square, 0)
  expect_equal(invariance$degrees_of_freedom, 0L)
  expect_equal(invariance$p_value, 0)
  expect_equal(invariance$bootstrap_nused, 1L)
  expect_equal(invariance$bootstrap_extreme_count, 1L)
  expect_equal(invariance$bootstrap_p_value, 1)

  replicate_row <- result$values$bootstrap$aggregate_replicates
  replicate_row <- replicate_row[replicate_row$kind == "invariance", , drop = FALSE]
  expect_equal(replicate_row$chi_square, 0)
  expect_equal(replicate_row$degrees_of_freedom, 0L)
  expect_equal(replicate_row$p_value, 0)
})

test_that("the preserved score-distribution resampling control is an explicit no-op", {
  fitted <- cm2_cm3_bootstrap_base_fit()
  fixed <- cm3(
    fitted,
    bootstrap = TRUE,
    nsim = 2L,
    seed = 59L,
    reestimate = FALSE,
    keep_bootstrap_samples = TRUE,
    resample_score_distribution = FALSE
  )$values$bootstrap
  requested <- cm3(
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
  fitted <- cm2_cm3_bootstrap_ld_dif_fit()
  context <- build_gllrm_context(fitted$spec, fitted$bundle)
  result <- cm3(
    fitted,
    bootstrap = TRUE,
    nsim = 2L,
    seed = 83L,
    reestimate = FALSE,
    keep_bootstrap_samples = TRUE
  )
  bootstrap <- result$values$bootstrap

  distribution <- new_cm2_cm3_bootstrap_distribution_cache(
    context,
    fitted$fit
  )(c(1L, 1L))
  # SKbias1.CollectComprecords loops the first component item outermost. This
  # differs from expand.grid()'s first-column-fastest order and determines the
  # seeded cumulative draw trajectory for tied component-score patterns.
  expect_identical(
    unname(distribution$components[[1L]]$configurations),
    matrix(c(0L, 0L, 0L, 1L, 1L, 0L, 1L, 1L), ncol = 2L, byrow = TRUE)
  )

  expect_true(bootstrap$possible)
  expect_identical(bootstrap$capability$observed$largest_component, 2L)
  expect_equal(bootstrap$capability$observed$dif_probability_strata, 2)
  expect_equal(bootstrap$capability$observed$exogenous_strata, 4)
  expect_identical(bootstrap$capability$bounds$maximum_items, 10L)
  expect_identical(bootstrap$capability$bounds$maximum_component_items, 4L)
  expect_identical(bootstrap$capability$bounds$maximum_multi_item_components, 4L)
  expect_identical(
    bootstrap$capability$bounds$maximum_positive_fixed_score_patterns,
    255L
  )
  expect_identical(bootstrap$capability$bounds$maximum_item_score, 7L)
  expect_identical(bootstrap$capability$bounds$maximum_dif_probability_strata, 64L)
  expect_identical(bootstrap$capability$bounds$maximum_exogenous_strata, 216L)
  expect_identical(bootstrap$capability$bounds$maximum_complete_records, 15000L)
  expect_identical(bootstrap$capability$observed$multi_item_components, 1L)
  expect_identical(
    bootstrap$capability$observed$largest_positive_fixed_score_pattern_count,
    2L
  )

  too_many_items <- context
  too_many_items$n_items <- 11L
  item_capability <- cm2_cm3_bootstrap_capability(too_many_items, fitted$fit)
  expect_false(item_capability$possible)
  expect_identical(item_capability$bounds$maximum_items, 10L)
  expect_true("more than 10 items" %in% item_capability$reasons)

  too_many_multicomponents <- context
  too_many_multicomponents$ld_components_items <- rep(
    list(c(1L, 2L)),
    5L
  )
  multicomponent_capability <- cm2_cm3_bootstrap_capability(
    too_many_multicomponents,
    fitted$fit
  )
  expect_false(multicomponent_capability$possible)
  expect_identical(
    multicomponent_capability$observed$multi_item_components,
    5L
  )
  expect_true(
    "more than four multi-item LD components" %in%
      multicomponent_capability$reasons
  )

  too_many_records <- context
  too_many_records$complete_item_exogenous_rows <- seq_len(15001L)
  record_capability <- cm2_cm3_bootstrap_capability(
    too_many_records,
    fitted$fit
  )
  expect_false(record_capability$possible)
  expect_identical(record_capability$observed$complete_records, 15001L)
  expect_true(
    "more than 15000 complete item/exogenous records" %in%
      record_capability$reasons
  )

  pattern_context <- context
  pattern_context$n_items <- 4L
  pattern_context$n_backgrounds <- 0L
  pattern_context$item_raw_max <- rep(8L, 4L)
  pattern_context$background_raw_max <- integer()
  pattern_context$dif_background_indices <- integer()
  pattern_context$max_total_score <- 28L
  pattern_context$ld_components_items <- list(seq_len(4L))
  pattern_key <- gllrm_component_key(seq_len(4L))
  pattern_configurations <- as.matrix(expand.grid(
    rep(list(0:7), 4L),
    KEEP.OUT.ATTRS = FALSE
  ))
  storage.mode(pattern_configurations) <- "integer"
  pattern_context$component_config_matrices <- stats::setNames(
    list(pattern_configurations),
    pattern_key
  )
  pattern_context$component_config_scores <- stats::setNames(
    list(rowSums(pattern_configurations)),
    pattern_key
  )
  pattern_context$component_ld_local_matrices <- stats::setNames(
    list(matrix(integer(), nrow = 0L, ncol = 3L)),
    pattern_key
  )
  pattern_context$component_ld_local_indices <-
    pattern_context$component_ld_local_matrices
  pattern_context$dif_by_item_matrices <- rep(
    list(matrix(integer(), nrow = 0L, ncol = 2L)),
    4L
  )
  pattern_context$dif_by_item <- pattern_context$dif_by_item_matrices
  pattern_state <- fitted$fit
  pattern_state$item_gamma <- matrix(1, nrow = 4L, ncol = 8L)
  pattern_state$ld_parameters <- list()
  pattern_state$dif_parameters <- list()
  pattern_capability <- cm2_cm3_bootstrap_capability(
    pattern_context,
    pattern_state
  )
  expect_false(pattern_capability$possible)
  expect_gt(
    pattern_capability$observed$largest_positive_fixed_score_pattern_count,
    255L
  )
  expect_true(
    "a fixed component score has more than 255 positive response patterns" %in%
      pattern_capability$reasons
  )

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
  fitted <- cm2_cm3_bootstrap_base_fit()
  result <- cm3(
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
  expect_identical(bootstrap$max_delta, 0.0001)
  expect_identical(bootstrap$acceptance_delta, 0.1)
  expect_true(replicate$reestimated)
  expect_false(is.na(replicate$iterations))
  expect_identical(replicate$accepted, replicate$final_delta < 0.1)
  expect_identical(bootstrap$nused, as.integer(replicate$accepted))
  expect_identical(bootstrap$samples, list())
})

test_that("bootstrap acceptance uses source stopping delta, not recomputed delta", {
  fitted <- cm2_cm3_bootstrap_base_fit()
  refitted_result <- cm3(
    fitted,
    bootstrap = TRUE,
    nsim = 1L,
    seed = 67L,
    reestimate = TRUE,
    bootstrap_max_step = 1L,
    keep_bootstrap_samples = TRUE
  )
  refitted_bootstrap <- refitted_result$values$bootstrap
  sample_bundle <- cm2_cm3_bootstrap_bundle(
    fitted,
    refitted_bootstrap$samples[[1L]]
  )
  spec <- if (is.null(fitted$model)) fitted$spec else fitted$model
  independent_refit <- fit_gllrm(
    spec,
    max_step = 1L,
    max_delta = 0.0001,
    bundle = sample_bundle,
    # SKbias8.MakeIRTcopy/skbias12b.UseItemParameters starts every source CM3
    # refit from the saved observed-data parameter structures.
    initial_parameters = gllrm_output_parameter_state(
      cm2_cm3_fit_context_state(fitted)$context,
      cm2_cm3_fit_context_state(fitted)$state
    )
  )

  # One source update leaves a large stopping discrepancy while the later
  # R-only sufficient-margin recomputation differs. The recorded source field
  # must remain the stopping value, and CM3 must reject on that value.
  expect_lt(independent_refit$state$delta, independent_refit$state$report_delta)
  expect_identical(
    refitted_bootstrap$replicates$report_delta,
    independent_refit$state$report_delta
  )
  expect_identical(
    refitted_bootstrap$replicates$final_delta,
    independent_refit$state$report_delta
  )
  expect_false(refitted_bootstrap$replicates$accepted)
  expect_identical(refitted_bootstrap$nused, 0L)

  limited_fit <- fit(
    gllrm(fitted$analysis),
    max_step = 1L,
    max_delta = 0.0001
  )
  expect_identical(limited_fit$fit$report_delta, 1)
  # The CML score-window correction changes this fixture's recomputed margin
  # discrepancy, but not the contract under test: the source stopping field is
  # the larger report_delta and CM3 must reject on that value.
  expect_lt(limited_fit$fit$delta, limited_fit$fit$report_delta)
  fixed_parameter_bootstrap <- cm3(
    limited_fit,
    bootstrap = TRUE,
    nsim = 1L,
    seed = 67L,
    reestimate = FALSE
  )$values$bootstrap
  expect_identical(fixed_parameter_bootstrap$replicates$final_delta, 1)
  expect_false(fixed_parameter_bootstrap$replicates$accepted)
  expect_identical(fixed_parameter_bootstrap$nused, 0L)
})

test_that("ordered parallel bootstrap refits equal serial refits exactly", {
  skip_on_os("windows")
  detected <- parallel::detectCores(logical = TRUE)
  skip_if(is.na(detected) || detected < 2L, "two fork-capable cores are required")

  fitted <- cm2_cm3_bootstrap_base_fit()
  arguments <- list(
    fit = fitted,
    bootstrap = TRUE,
    nsim = 2L,
    seed = 67L,
    reestimate = TRUE,
    bootstrap_max_step = 50L,
    keep_bootstrap_samples = FALSE
  )
  serial <- do.call(cm3, c(arguments, list(bootstrap_jobs = 1L)))$values$bootstrap

  set.seed(1701L)
  global_state <- .Random.seed
  forked <- do.call(cm3, c(arguments, list(bootstrap_jobs = 4L)))$values$bootstrap
  expect_identical(.Random.seed, global_state)

  # Generation precedes the fork and uses the same private Delphi stream. The
  # ordered worker map performs only source refits and diagnostic arithmetic,
  # so every substantive field is bit-for-bit equal to the serial run.
  expect_identical(forked$seed, serial$seed)
  expect_identical(forked$final_rng_state, serial$final_rng_state)
  expect_identical(forked$rng_draws, serial$rng_draws)
  expect_identical(forked$replicates, serial$replicates)
  expect_identical(forked$aggregate_replicates, serial$aggregate_replicates)
  expect_identical(forked$margin_replicates, serial$margin_replicates)
  expect_identical(forked$aggregate_summary, serial$aggregate_summary)
  expect_identical(forked$margin_summary, serial$margin_summary)
  expect_identical(forked$samples, list())
  expect_identical(serial$requested_bootstrap_jobs, 1L)
  expect_identical(serial$bootstrap_jobs, 1L)
  expect_identical(serial$execution_mode, "serial")
  expect_identical(forked$requested_bootstrap_jobs, 4L)
  expect_identical(forked$bootstrap_jobs, 2L)
  expect_identical(forked$execution_mode, "ordered_fork_parallel")

  detected <- parallel::detectCores(logical = TRUE)
  if (is.na(detected) || detected < 1L) {
    detected <- 1L
  }
  expect_identical(
    cm2_cm3_bootstrap_worker_count(128L, 1000L),
    as.integer(min(128L, 1000L, detected))
  )
})

test_that("CM2 exposes calibrated aggregate and margin tables without inventing CM3", {
  fitted <- cm2_cm3_bootstrap_base_fit()
  result <- cm2(
    fitted,
    bootstrap = TRUE,
    nsim = 1L,
    seed = 71L,
    reestimate = FALSE
  )

  expect_s3_class(result, "gRm_cm2")
  expect_true("bootstrap_p_value" %in% names(result$values$aggregate))
  expect_true("bootstrap_p_value" %in% names(result$values$item_trait))
  expect_true("bootstrap_p_value" %in% names(result$values$tests))
  expect_false("CM3" %in% result$values$bootstrap$aggregate_summary$diagnostic)
})

test_that("bootstrap exceedance summaries use accepted samples and source p direction", {
  observed <- data.frame(
    diagnostic = c("CM2", "CM3"),
    p_value = c(0.2, 0.5),
    stringsAsFactors = FALSE
  )
  simulated <- data.frame(
    diagnostic = rep(c("CM2", "CM3"), each = 3L),
    accepted = rep(c(TRUE, TRUE, FALSE), 2L),
    p_value = c(0.2, 0.3, 0, 0.1, 0.5, 0),
    stringsAsFactors = FALSE
  )

  summary <- cm2_cm3_bootstrap_summarize_rows(
    observed,
    simulated,
    key = "diagnostic",
    nused = 2L
  )

  expect_identical(summary$bootstrap_extreme_count, c(1L, 2L))
  expect_identical(summary$bootstrap_nused, c(2L, 2L))
  expect_equal(summary$bootstrap_p_value, c(0.5, 1))
})
