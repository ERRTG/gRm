cm2_cm3_api_data <- function() {
  rows <- expand.grid(
    I1 = 1:3,
    I2 = 1:3,
    I3 = 1:3,
    site = 1:2,
    KEEP.OUT.ATTRS = FALSE
  )
  data.frame(ID = seq_len(nrow(rows)), rows)
}

cm2_cm3_api_analysis <- function(score_cuts = c(2L, 6L)) {
  gRm(
    cm2_cm3_api_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID",
    score_cuts = score_cuts
  )
}

cm2_cm3_api_fit <- function(score_cuts = c(2L, 6L)) {
  fit(gllrm(cm2_cm3_api_analysis(score_cuts = score_cuts)), max_step = 100L)
}

test_that("cm2 and cm3 expose source-named public diagnostics", {
  expect_true(exists("cm2", envir = asNamespace("gRm"), inherits = FALSE))
  expect_true(exists("cm3", envir = asNamespace("gRm"), inherits = FALSE))
  expect_false("which" %in% names(formals(gRm::cm2)))
  expect_false("which" %in% names(formals(gRm::cm3)))
  expect_true("items" %in% names(formals(gRm::cm2)))
  expect_true("items" %in% names(formals(gRm::cm3)))
  expect_true("score_cuts" %in% names(formals(gRm::cm2)))
  expect_true("score_cuts" %in% names(formals(gRm::cm3)))
})

test_that("cm2 and cm3 require fitted model objects", {
  analysis <- cm2_cm3_api_analysis()
  fit0 <- fit(gllrm(analysis), max_step = 100L)

  expect_error(cm2(analysis), "fitted", fixed = TRUE)
  expect_error(cm3(analysis), "fitted", fixed = TRUE)
  expect_s3_class(cm2(fit0), "gRm_cm2")
  expect_s3_class(cm3(fit0), "gRm_cm3")
})

test_that("cm2 and cm3 score cuts are diagnostic-local source state", {
  fit0 <- cm2_cm3_api_fit(score_cuts = c(2L, 6L))
  analysis_cuts <- fit0$analysis$score_groups
  explicit_cuts <- c(3L, 6L)

  cm2_result <- cm2(fit0, score_cuts = explicit_cuts)
  cm3_result <- cm3(fit0, score_cuts = explicit_cuts)

  expect_equal(cm2_result$metadata$score_cuts, explicit_cuts)
  expect_equal(cm3_result$metadata$score_cuts, explicit_cuts)
  expect_equal(cm2_result$values$score_cuts, explicit_cuts)
  expect_equal(cm3_result$values$score_cuts, explicit_cuts)
  expect_equal(fit0$analysis$score_groups, analysis_cuts)
})

test_that("cm2 and cm3 require one-dimensional score-cut vectors", {
  fit0 <- cm2_cm3_api_fit(score_cuts = c(2L, 6L))

  expect_error(
    cm2(fit0, score_cuts = matrix(c(2L, 6L), nrow = 1L)),
    "one-dimensional"
  )
  expect_error(
    cm3(fit0, score_cuts = array(c(2L, 6L), dim = c(2L, 1L, 1L))),
    "one-dimensional"
  )
})

test_that("cm2 and cm3 normalize selected items to source order", {
  fit0 <- cm2_cm3_api_fit()

  cm2_result <- cm2(fit0, items = c("I3", "I1"))
  cm3_result <- cm3(fit0, items = c(3L, 1L))

  expect_equal(cm2_result$metadata$item_indices, c(1L, 3L))
  expect_equal(cm3_result$metadata$item_indices, c(1L, 3L))
  expect_equal(cm2_result$metadata$items, c("I1", "I3"))
  expect_equal(cm3_result$metadata$items, c("I1", "I3"))
  expect_identical(cm2_result$metadata$selection, cm2_result$values$selection)
  expect_identical(cm3_result$metadata$selection, cm3_result$values$selection)
  expect_equal(cm2_result$metadata$selection$schema_version, 1L)
  expect_equal(cm2_result$metadata$selection$mode, "explicit_names")
  expect_equal(cm2_result$metadata$selection$requested_items, c("I3", "I1"))
  expect_equal(cm3_result$metadata$selection$mode, "explicit_indices")
  expect_equal(cm3_result$metadata$selection$requested_items, c(3L, 1L))
  expect_equal(cm3_result$metadata$selection$resolved_items$item_name, c("I1", "I3"))
  expect_equal(cm3_result$metadata$selection$selected_count, 2L)
  expect_equal(cm3_result$metadata$selection$model_item_count, 3L)
  expect_equal(cm3_result$metadata$selection$exogenous_scope, "all_fitted")
  expect_equal(cm3_result$metadata$selection$exogenous$exogenous_name, "site")
  expect_true(cm3_result$metadata$selection$score_group_included)
  expect_equal(cm3_result$metadata$selection$score_total_scope, "all_fitted_items")
  expect_true(all(cm2_result$values$cm2_bh$n_p_values == nrow(cm2_result$values$tests)))
  expect_true(all(
    cm3_result$values$cm3_bh$n_p_values == sum(!cm3_result$values$tests$is_cm2)
  ))
})

test_that("NULL, exact all-name, and exact all-index selections are numerically equivalent", {
  fit0 <- cm2_cm3_api_fit()

  default <- cm3(fit0)
  named <- cm3(fit0, items = c("I3", "I1", "I2"))
  indexed <- cm3(fit0, items = c(3L, 1L, 2L))

  expect_equal(default$values$tests, named$values$tests)
  expect_equal(default$values$tests, indexed$values$tests)
  expect_equal(default$values$cm2, named$values$cm2)
  expect_equal(default$values$cm3, indexed$values$cm3)
  expect_equal(default$metadata$selection$mode, "default_all")
  expect_equal(named$metadata$selection$mode, "explicit_names")
  expect_equal(indexed$metadata$selection$mode, "explicit_indices")
})

test_that("explicit selection is diagnostic-local and CM2 is the two-way CM3 projection", {
  fit0 <- cm2_cm3_api_fit()
  before <- serialize(fit0, NULL)

  cm2_result <- cm2(fit0, items = c("I3", "I1"))
  cm3_result <- cm3(fit0, items = c("I3", "I1"))

  expect_identical(serialize(fit0, NULL), before)
  expect_equal(
    cm2_result$values$tests,
    cm3_result$values$tests[cm3_result$values$tests$is_cm2, , drop = FALSE]
  )
  expect_equal(cm2_result$values$aggregate, cm3_result$values$cm2)
  expect_equal(cm2_result$values$item_trait, cm3_result$values$item_trait)
  expect_equal(cm2_result$values$invariance, cm3_result$values$invariance)
})

test_that("cm2 and cm3 print and summary use planned headers", {
  fit0 <- cm2_cm3_api_fit()
  cm2_result <- cm2(fit0)
  cm3_result <- cm3(fit0)

  expect_equal(cm2_result$values$cm2_bh$fdr, c(0.05, 0.01, 0.001))
  expect_true(all(cm2_result$values$cm2_bh$n_p_values == cm2_result$values$n_two_way_margins))
  expect_equal(cm3_result$values$cm3_bh$fdr, c(0.05, 0.01, 0.001))
  expect_true(all(cm3_result$values$cm3_bh$n_p_values == cm3_result$values$n_three_way_margins))
  expect_true("bh" %in% names(summary(cm3_result)$tables))

  expect_match(capture.output(print(cm2_result))[[1L]], "^gRm: CM2 fit diagnostic$")
  expect_match(capture.output(print(cm3_result))[[1L]], "^gRm: CM3 fit diagnostic$")
  expect_match(capture.output(summary(cm2_result))[[1L]], "^gRm: CM2 fit diagnostic$")
  expect_match(capture.output(summary(cm3_result))[[1L]], "^gRm: CM3 fit diagnostic$")
  expect_true(any(grepl("^  Score cuts: ", capture.output(print(cm2_result)))))
  expect_true(any(grepl("^Score cuts: ", capture.output(summary(cm3_result)))))
  expect_true(any(grepl("^  Selection mode: default_all$", capture.output(print(cm2_result)))))
  expect_true(any(grepl("^Selection: default_all;", capture.output(summary(cm3_result)))))
  expect_true(any(grepl("score total: all fitted items", capture.output(summary(cm3_result)))))
})

test_that("cm2 and cm3 print and summary use diagnostic numeric formatting", {
  old_digits <- getOption("digits")
  options(digits = 7L)
  on.exit(options(digits = old_digits), add = TRUE)

  stat <- data.frame(
    diagnostic = "CM2",
    chi_square = 12.3456789,
    degrees_of_freedom = 4L,
    p_value = 0.000123456,
    stringsAsFactors = FALSE
  )
  values <- list(
    selected_items = data.frame(item_name = "I1", item_index = 1L),
    exogenous_names = "site",
    score_cuts = c(2L, 6L),
    n_two_way_margins = 1L,
    n_three_way_margins = 0L,
    aggregate = stat,
    cm2 = stat,
    cm3 = transform(stat, diagnostic = "CM3"),
    item_trait = transform(stat, diagnostic = "Item-trait"),
    invariance = data.frame(),
    tests = transform(
      stat,
      margin = "I1 x I2",
      margin_type = "item_item",
      variable1 = "I1",
      variable2 = "I2",
      variable3 = NA_character_,
      is_cm2 = TRUE
    )
  )
  result <- structure(list(values = values), class = "gRm_cm2")

  printed <- capture.output(print(result))
  expect_true(any(grepl("^  CM2 Chisq: 12.3$", printed)))
  expect_false(any(grepl("12.3457", printed, fixed = TRUE)))

  summary_printed <- capture.output(print(summary(result)))
  expect_true(any(grepl("12.3", summary_printed, fixed = TRUE)))
  expect_false(any(grepl("12.3457", summary_printed, fixed = TRUE)))
  expect_true(any(grepl("0.000123", summary_printed, fixed = TRUE)))
})

test_that("CM2 and CM3 format bootstrap and audit fields without changing values", {
  old_options <- options(digits = 7L, scipen = 0L, width = 120L)
  on.exit(options(old_options), add = TRUE)

  bootstrap_p <- c(0, 1 / 3, NA_real_)
  critical_p <- c(0.000123456789, 0, NA_real_)
  seed <- c(4294967295, 60437825, NA_real_)
  acceptance_delta <- c(0.123456789, 0.1, NA_real_)
  raw_display <- data.frame(
    `Bootstrap Pr` = bootstrap_p,
    `Critical p` = critical_p,
    Seed = seed,
    `Acceptance delta` = acceptance_delta,
    check.names = FALSE
  )
  formatted <- format_cm2_cm3_summary_table(raw_display)

  expect_identical(
    formatted[["Bootstrap Pr"]],
    summary_p_values(bootstrap_p, empirical = TRUE)
  )
  expect_identical(formatted[["Bootstrap Pr"]], c("0", "0.333", "NA"))
  expect_identical(formatted[["Critical p"]], summary_p_values(critical_p))
  expect_identical(formatted[["Critical p"]], c("0.000123", "< 2e-16", "NA"))
  expect_identical(formatted$Seed, c("4294967295", "60437825", NA_character_))
  expect_identical(
    formatted[["Acceptance delta"]],
    vapply(acceptance_delta, summary_scalar, character(1L))
  )
  expect_type(raw_display[["Bootstrap Pr"]], "double")
  expect_identical(raw_display$Seed, seed)

  stat <- data.frame(
    diagnostic = "CM2",
    chi_square = 12.3456789,
    degrees_of_freedom = 4L,
    p_value = 0.000123456789,
    bootstrap_p_value = 0,
    bootstrap_nused = 999L,
    stringsAsFactors = FALSE
  )
  values <- list(
    selection = list(mode = "default_all"),
    selected_items = data.frame(item_name = "I1", item_index = 1L),
    exogenous_names = "site",
    score_cuts = c(2L, 6L),
    n_two_way_margins = 1L,
    n_three_way_margins = 1L,
    aggregate = stat,
    cm2 = stat,
    cm3 = transform(stat, diagnostic = "CM3", bootstrap_p_value = 1 / 3),
    item_trait = transform(stat, diagnostic = "Item-trait"),
    invariance = transform(stat, background_name = "site"),
    tests = transform(
      stat,
      margin = "I1 x I2",
      margin_type = "item_item",
      variable1 = "I1",
      variable2 = "I2",
      variable3 = NA_character_,
      is_cm2 = TRUE
    ),
    cm2_bh = data.frame(
      diagnostic = "CM2",
      fdr = 0.05,
      critical_p = 0.000123456789,
      n_p_values = 1L,
      source_limit_reached = FALSE
    ),
    cm3_bh = data.frame(
      diagnostic = "CM3",
      fdr = 0.05,
      critical_p = 0.000123456789,
      n_p_values = 1L,
      source_limit_reached = FALSE
    ),
    bootstrap = list(
      enabled = TRUE,
      possible = TRUE,
      source_status = "complete",
      nsim = 1000L,
      nused = 999L,
      seed = 4294967295,
      reestimate = TRUE,
      bootstrap_jobs = 2L,
      acceptance_delta = 0.123456789
    )
  )

  for (class_name in c("gRm_cm2", "gRm_cm3")) {
    result <- structure(list(values = values), class = class_name)
    compact <- capture.output(print(result))
    summarized <- summary(result)
    aggregate_values <- summarized$tables$aggregates[["Bootstrap Pr"]]
    seed_value <- summarized$tables$bootstrap$Seed
    tables_before_print <- serialize(summarized$tables, NULL, version = 3L)
    full <- capture.output(print(summarized))

    expect_true(any(grepl("Bootstrap Pr: 0$", compact)))
    expect_false(any(grepl("Bootstrap Pr: <", compact, fixed = TRUE)))
    expect_true(any(grepl("0.000123", full, fixed = TRUE)))
    expect_false(any(grepl("0.000123456", full, fixed = TRUE)))
    expect_true(any(grepl(
      "^ +CM2 +12[.]3 +4 +0[.]000123 +0 +999$",
      full
    )))
    expect_true(any(grepl("4294967295", full, fixed = TRUE)))
    expect_false(any(grepl("4.29497e+09", full, fixed = TRUE)))
    expect_type(aggregate_values, "double")
    expect_identical(seed_value, 4294967295)
    expect_identical(
      summarized$tables$aggregates[["Bootstrap Pr"]],
      aggregate_values
    )
    expect_identical(
      serialize(summarized$tables, NULL, version = 3L),
      tables_before_print
    )
  }

  cm3_result <- structure(list(values = values), class = "gRm_cm3")
  cm3_printed <- capture.output(print(cm3_result))
  expect_true(any(grepl("CM3 Bootstrap Pr: 0.333$", cm3_printed)))

  options(digits = 3L)
  low_digits <- capture.output(print(summary(cm3_result)))
  expect_true(any(grepl("4294967295", low_digits, fixed = TRUE)))
  expect_false(any(grepl("4.29e+09", low_digits, fixed = TRUE)))
})

test_that("observed CM2 and CM3 BH blocks preserve the source capacity", {
  rows <- data.frame(p_value = c(rep(1, 85L), 0))
  bh <- cm2_cm3_bh_thresholds(rows, "CM3")

  expect_equal(bh$diagnostic, rep("CM3", 3L))
  expect_equal(bh$fdr, c(0.05, 0.01, 0.001))
  expect_equal(bh$n_p_values, rep(85L, 3L))
  expect_true(all(bh$source_limit_reached))
  expect_equal(bh$critical_p, bh$fdr / 85)

  expect_equal(nrow(cm2_cm3_bh_thresholds(data.frame(), "CM2")), 0L)
})
