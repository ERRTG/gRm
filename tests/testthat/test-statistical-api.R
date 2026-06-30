stat_api_data <- function() {
  data.frame(
    ID = seq_len(12L),
    I1 = c(1L, 1L, 2L, 2L, 3L, 3L, 1L, 2L, 3L, 1L, 2L, 3L),
    I2 = c(1L, 2L, 1L, 3L, 2L, 3L, 3L, 1L, 2L, 2L, 3L, 1L),
    I3 = c(2L, 1L, 3L, 1L, 2L, 3L, 2L, 1L, 3L, 3L, 2L, 1L),
    site = c(1L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 2L)
  )
}

stat_api_analysis <- function() {
  gRm(
    stat_api_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID",
    score_cuts = "auto"
  )
}

test_that("gRm constructs the public analysis object", {
  analysis <- stat_api_analysis()

  expect_s3_class(analysis, "gRm_analysis")
  expect_equal(analysis$items, c("I1", "I2", "I3"))
  expect_equal(analysis$exogenous, "site")
  expect_equal(analysis$id, "ID")
  expect_false(inherits(analysis, "gRm_item_analysis"))
})

test_that("score_cuts define analysis score groups and global homogeneity defaults", {
  analysis <- gRm(
    stat_api_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID",
    score_cuts = c(2L, 6L)
  )
  fitted <- fit(gllrm(analysis), max_step = 50L)
  default <- global_homogeneity(fitted, max_step = 50L)
  override <- global_homogeneity(fitted, score_cuts = c(3L, 6L), max_step = 50L)

  expect_equal(analysis$score_groups, c(2L, 6L))
  expect_equal(default$metadata$score_cuts, c(2L, 6L))
  expect_equal(default$values$score_groups$to_score, c(2L, 6L))
  expect_equal(override$metadata$score_cuts, c(3L, 6L))
  expect_equal(override$values$score_groups$to_score, c(3L, 6L))
  expect_false("groups" %in% names(formals(gRm)))
  expect_false("groups" %in% names(formals(read_digram_project)))
  expect_false("groups" %in% names(formals(global_homogeneity)))
})

test_that("global homogeneity rejects invalid explicit score cut overrides", {
  analysis <- gRm(
    stat_api_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID",
    score_cuts = c(2L, 6L)
  )
  fitted <- fit(gllrm(analysis), max_step = 50L)

  expect_error(
    global_homogeneity(fitted, score_cuts = 2L, max_step = 20L),
    "at least two"
  )
  expect_error(
    global_homogeneity(fitted, score_cuts = c(2.5, 6), max_step = 20L),
    "integer-like"
  )
  expect_error(
    global_homogeneity(fitted, score_cuts = c(2L, 2L), max_step = 20L),
    "strictly increasing"
  )
  expect_error(
    global_homogeneity(fitted, score_cuts = c(-1L, 6L), max_step = 20L),
    "possible score range"
  )
  expect_error(
    global_homogeneity(fitted, score_cuts = c(2L, 7L), max_step = 20L),
    "possible score range"
  )
  expect_error(
    global_homogeneity(fitted, score_cuts = c(0L, 6L), max_step = 20L),
    "at least two global-homogeneity score groups"
  )
  expect_error(
    global_homogeneity(fitted, score_cuts = c(2, .Machine$integer.max + 1), max_step = 20L),
    "integer-like score cuts"
  )
})

test_that("gllrm returns one model class for Rasch and GLLRM specifications", {
  analysis <- stat_api_analysis()
  rasch <- gllrm(analysis)
  gllrm_model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I3:site)

  expect_s3_class(rasch, "gRm_model")
  expect_s3_class(gllrm_model, "gRm_model")
  expect_equal(rasch$model_type, "rasch")
  expect_equal(gllrm_model$model_type, "gllrm")
  expect_equal(gllrm_model$terms$ld$item1, "I1")
  expect_equal(gllrm_model$terms$ld$item2, "I2")
  expect_equal(gllrm_model$terms$dif$item, "I3")
  expect_equal(gllrm_model$terms$dif$exogenous, "site")
})

test_that("gllrm formulas support backticked non-syntactic variable names", {
  data <- data.frame(
    ID = seq_len(8L),
    check.names = FALSE,
    "item one" = c(0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L),
    "item two" = c(1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L),
    "item three" = c(0L, 0L, 1L, 1L, 1L, 0L, 1L, 0L),
    "clinic group" = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  )
  analysis <- gRm(
    data,
    items = c("item one", "item two", "item three"),
    exogenous = "clinic group",
    id = "ID",
    score_cuts = c(1L, 3L)
  )

  model <- gllrm(
    analysis,
    ld = ~ `item one`:`item two`,
    dif = ~ `item three`:`clinic group`
  )

  expect_equal(model$terms$ld$item1, "item one")
  expect_equal(model$terms$ld$item2, "item two")
  expect_equal(model$terms$dif$item, "item three")
  expect_equal(model$terms$dif$exogenous, "clinic group")
})

test_that("gllrm formulas reject expanded model syntax", {
  analysis <- stat_api_analysis()

  expect_equal(gllrm(analysis, ld = ~ 0 + I1:I2)$terms$ld$item1, "I1")
  expect_equal(gllrm(analysis, ld = ~ -1 + I1:I2)$terms$ld$item2, "I2")
  expect_error(
    gllrm(analysis, ld = ~ I1 * I2),
    "allow only `:` interaction terms joined by"
  )
  expect_error(
    gllrm(analysis, dif = ~ I3:(site + I1)),
    "allow only `:` interaction terms joined by"
  )
})

test_that("fit returns one public fit class with model-type metadata", {
  analysis <- stat_api_analysis()
  rasch_fit <- fit(gllrm(analysis), max_step = 50L)
  gllrm_fit <- fit(gllrm(analysis, ld = ~ I1:I2, dif = ~ I3:site), max_step = 50L)

  expect_s3_class(rasch_fit, "gRm_fit")
  expect_s3_class(gllrm_fit, "gRm_fit")
  expect_equal(rasch_fit$model_type, "rasch")
  expect_equal(gllrm_fit$model_type, "gllrm")
  expect_false(inherits(rasch_fit, "gRm_gllrm_fit"))
})

test_that("fit controls are validated before fitting paths run", {
  analysis <- stat_api_analysis()
  model <- gllrm(analysis)
  gllrm_model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I3:site)

  bad_steps <- list(0L, -1L, NA_integer_, 1.5, Inf, c(1L, 2L))
  for (bad_step in bad_steps) {
    expect_error(
      fit(model, max_step = bad_step),
      "`max_step` must be a single positive integer-like value"
    )
  }

  bad_deltas <- list(0, -0.1, NA_real_, Inf, c(0.1, 0.2))
  for (bad_delta in bad_deltas) {
    expect_error(
      fit(model, max_step = 1L, max_delta = bad_delta),
      "`max_delta` must be a single positive finite number"
    )
  }

  expect_no_error(fit(model, max_step = 1, max_delta = 0.1))
  expect_error(
    fit(gllrm_model, max_step = 0L),
    "`max_step` must be a single positive integer-like value"
  )
  expect_error(
    fit(gllrm_model, max_step = 1L, max_delta = 0),
    "`max_delta` must be a single positive finite number"
  )
})

test_that("diagnostic refit controls use the same validation rules", {
  fitted <- fit(gllrm(stat_api_analysis()), max_step = 50L)

  expect_error(
    local_dependence(fitted, max_step = 0L),
    "`max_step` must be a single positive integer-like value"
  )
  expect_error(
    dif(fitted, max_step = 0L),
    "`max_step` must be a single positive integer-like value"
  )
  expect_error(
    global_homogeneity(fitted, max_step = 0L),
    "`max_step` must be a single positive integer-like value"
  )

  expect_error(
    local_dependence(fitted, max_delta = 0),
    "`max_delta` must be a single positive finite number"
  )
  expect_error(
    dif(fitted, max_delta = 0),
    "`max_delta` must be a single positive finite number"
  )
  expect_error(
    global_homogeneity(fitted, max_delta = 0),
    "`max_delta` must be a single positive finite number"
  )
})

test_that("screen keeps model discovery separate from score effects", {
  analysis <- stat_api_analysis()
  result <- screen(analysis, inference = "asymptotic")
  effects <- score_effects(analysis)

  expect_s3_class(result, "gRm_screen")
  expect_s3_class(effects, "gRm_score_effects")
  expect_s3_class(effects, "gRm_direct_table")
  expect_false("score_effects" %in% names(result))
  result_summary <- summary(result)
  expect_true(is.data.frame(result_summary$local_dependence))
  expect_true(is.data.frame(result_summary$dif))
  expect_true(is.data.frame(result_summary$score_effects))
  expect_true(is.data.frame(result_summary$selected))
  expect_true(is.data.frame(attr(result_summary, "bh", exact = TRUE)))
  expect_true(is.data.frame(effects))
  expect_true(is.data.frame(attr(effects, "selected", exact = TRUE)))
  expect_true(is.data.frame(attr(effects, "bh", exact = TRUE)))
})

test_that("screen maps exact inference modes to DIGRAM source command state", {
  analysis <- stat_api_analysis()

  exact <- screen(analysis, inference = "exact", nsim = 25L, seed = 9L)
  repeated <- screen(analysis, inference = "repeated", nsim = 25L, seed = 9L)

  expect_equal(exact$exact_state$command_no, 2L)
  expect_equal(exact$exact_state$seq_limit, 25L)
  expect_equal(repeated$exact_state$command_no, 74L)
  expect_equal(repeated$exact_state$seq_limit, 20L)
  expect_equal(repeated$exact_state$seq_alpha, 0.001)
  expect_error(
    screen(analysis, inference = "sequential", nsim = 25L, seed = 9L),
    "should be one of"
  )
})

test_that("score_effects maps exact inference modes to DIGRAM source command state", {
  analysis <- stat_api_analysis()

  exact <- score_effects(analysis, inference = "exact", nsim = 25L, seed = 9L)
  repeated <- score_effects(analysis, inference = "repeated", nsim = 25L, seed = 9L)

  expect_true(all(c(
    "Exact Pr(>Chisq)", "Exact Pr(Gamma+)", "Exact Pr(|Gamma|)", "Simulations"
  ) %in% names(exact)))
  expect_true(all(c(
    "Exact Pr(>Chisq)", "Exact Pr(Gamma+)", "Exact Pr(|Gamma|)", "Simulations"
  ) %in% names(repeated)))
  exact_metadata <- attr(exact, "metadata", exact = TRUE)
  repeated_metadata <- attr(repeated, "metadata", exact = TRUE)
  expect_equal(exact_metadata$exact_state$command_no, 2L)
  expect_equal(exact_metadata$exact_state$seq_limit, 25L)
  expect_equal(repeated_metadata$exact_state$command_no, 74L)
  expect_equal(repeated_metadata$exact_state$seq_limit, 20L)
  expect_equal(repeated_metadata$exact_state$seq_alpha, 0.001)
  expect_error(
    score_effects(analysis, inference = "sequential", nsim = 25L, seed = 9L),
    "should be one of"
  )
})

test_that("public exact inference arguments reject malformed scalars", {
  analysis <- stat_api_analysis()

  expect_error(
    screen(analysis, inference = "exact", nsim = 1.9),
    "`nsim` must be a single non-negative integer-like value"
  )
  expect_error(
    screen(analysis, inference = "exact", nsim = c(25L, 999L)),
    "`nsim` must be a single non-negative integer-like value"
  )
  expect_error(
    screen(analysis, inference = "exact", seed = 9.9),
    "`seed` must be a single integer-like value"
  )
  expect_error(
    screen(analysis, inference = "repeated", critlevel = 25.5),
    "`critlevel` must be a single non-negative integer-like value"
  )
  expect_error(
    screen(analysis, inference = "repeated", risk = c(1L, 2L)),
    "`risk` must be a single non-negative integer-like value"
  )
  expect_error(
    score_effects(analysis, inference = "exact", nsim = 1.9),
    "`nsim` must be a single non-negative integer-like value"
  )
  expect_error(
    score_effects(analysis, inference = "repeated", risk = c(1L, 2L)),
    "`risk` must be a single non-negative integer-like value"
  )
})

test_that("non-parallel APIs do not expose jobs", {
  analysis <- stat_api_analysis()

  expect_false("jobs" %in% names(formals(score_effects)))
  expect_false("jobs" %in% names(formals(screen.gRm_analysis)))
  expect_false("jobs" %in% names(formals(screen.gRm_project)))
  expect_false("jobs" %in% names(formals(global_homogeneity)))
  expect_error(
    score_effects(analysis, jobs = 1L),
    "unused"
  )
  expect_error(
    screen(analysis, inference = "asymptotic", jobs = 1L),
    "reserved"
  )
  expect_error(
    global_homogeneity(NULL, jobs = 1L),
    "reserved"
  )
})

test_that("public integer-like validators reject unsafe coercions", {
  huge <- .Machine$integer.max + 1

  expect_error(
    gRm(
      stat_api_data(),
      items = c("I1", "I2", "I3"),
      exogenous = "site",
      id = "ID",
      score_cuts = c(2, huge)
    ),
    "integer-like values"
  )
})

test_that("score_effects does not expose the source score cap as a public argument", {
  analysis <- stat_api_analysis()

  expect_false("score_cap" %in% names(formals(score_effects)))
  expect_error(
    score_effects(analysis, score_cap = 2L),
    "unused"
  )

  effects <- score_effects(analysis)
  expect_equal(attr(effects, "metadata", exact = TRUE)$score_cap, 56L)
  expect_equal(attr(effects, "values", exact = TRUE)$score_cap, 56L)
})

test_that("post-fit diagnostic functions return public diagnostic classes", {
  fitted <- fit(gllrm(stat_api_analysis()), max_step = 50L)

  item_tests <- item_fit(fitted)
  ld_tests <- local_dependence(fitted)
  dif_tests <- dif(fitted)
  gh_tests <- global_homogeneity(fitted)
  ari_table <- ari(fitted)

  expect_false("item_parameters" %in% getNamespaceExports("gRm"))
  expect_s3_class(item_tests, "gRm_item_fit")
  expect_s3_class(ld_tests, "gRm_local_dependence")
  expect_s3_class(dif_tests, "gRm_dif")
  expect_s3_class(gh_tests, "gRm_global_homogeneity")
  expect_s3_class(ari_table, "gRm_ari")
  expect_true(is.data.frame(summary(fitted, which = "parameters")$parameters))
  expect_true(is.data.frame(summary(fitted, which = "thresholds")$thresholds))
  expect_true(is.data.frame(item_tests))
  expect_s3_class(attr(item_tests, "values", exact = TRUE), "gRm_item_fits_values")
  expect_true(is.data.frame(summary(ld_tests)$tests))
  expect_true(is.data.frame(summary(dif_tests)$tests))
  expect_true(is.data.frame(summary(gh_tests, which = "test")$test))
  expect_true(is.data.frame(ari_table))
})

test_that("removed helpers and accessors are not available as public or namespace functions", {
  removed <- c(
    "sum_score",
    "score_groups_auto",
    "score_groups_cut",
    "check",
    "model_terms",
    "tidy",
    "glance",
    "details",
    "detail_names",
    "status",
    "source_trace",
    "validation",
    "unmodeled",
    "gRm_warnings"
  )

  expect_false(any(removed %in% getNamespaceExports("gRm")))

  ns <- asNamespace("gRm")
  namespace_absent <- c(
    "sum_score",
    "score_groups_auto",
    "score_groups_cut",
    "validate_score_spec",
    "validate_score_group_spec",
    "resolve_gRm_score_groups",
    "tidy",
    "glance",
    "details",
    "detail_names"
  )
  for (name in namespace_absent) {
    expect_false(exists(name, envir = ns, inherits = FALSE), info = name)
  }
})
