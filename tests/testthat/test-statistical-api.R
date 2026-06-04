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

test_that("screen keeps model discovery separate from score effects", {
  analysis <- stat_api_analysis()
  result <- screen(analysis, inference = "asymptotic", jobs = 1L)
  effects <- score_effects(analysis, jobs = 1L)

  expect_s3_class(result, "gRm_screen")
  expect_s3_class(effects, "gRm_score_effects")
  expect_false("score_effects" %in% names(result))
  expect_true(is.data.frame(summary(result)$tests))
  expect_true(is.data.frame(summary(effects)$tests))
})

test_that("screen maps exact inference modes to DIGRAM source command state", {
  analysis <- stat_api_analysis()

  exact <- screen(analysis, inference = "exact", nsim = 25L, seed = 9L, jobs = 1L)
  repeated <- screen(analysis, inference = "repeated", nsim = 25L, seed = 9L, jobs = 1L)

  expect_equal(exact$exact_state$command_no, 2L)
  expect_equal(exact$exact_state$seq_limit, 25L)
  expect_equal(repeated$exact_state$command_no, 74L)
  expect_equal(repeated$exact_state$seq_limit, 20L)
  expect_equal(repeated$exact_state$seq_alpha, 0.001)
  expect_error(
    screen(analysis, inference = "sequential", nsim = 25L, seed = 9L, jobs = 1L),
    "should be one of"
  )
})

test_that("score_effects maps exact inference modes to DIGRAM source command state", {
  analysis <- stat_api_analysis()

  exact <- score_effects(analysis, inference = "exact", nsim = 25L, seed = 9L, jobs = 1L)
  repeated <- score_effects(analysis, inference = "repeated", nsim = 25L, seed = 9L, jobs = 1L)

  expect_equal(exact$metadata$exact_state$command_no, 2L)
  expect_equal(exact$metadata$exact_state$seq_limit, 25L)
  expect_equal(repeated$metadata$exact_state$command_no, 74L)
  expect_equal(repeated$metadata$exact_state$seq_limit, 20L)
  expect_equal(repeated$metadata$exact_state$seq_alpha, 0.001)
  expect_error(
    score_effects(analysis, inference = "sequential", nsim = 25L, seed = 9L, jobs = 1L),
    "should be one of"
  )
})

test_that("post-fit diagnostic functions return public diagnostic classes", {
  fitted <- fit(gllrm(stat_api_analysis()), max_step = 50L)

  parameters <- item_parameters(fitted)
  item_tests <- item_fit(fitted)
  ld_tests <- local_dependence(fitted)
  dif_tests <- dif(fitted)
  gh_tests <- global_homogeneity(fitted)

  expect_s3_class(parameters, "gRm_item_parameters")
  expect_s3_class(item_tests, "gRm_item_fit")
  expect_s3_class(ld_tests, "gRm_local_dependence")
  expect_s3_class(dif_tests, "gRm_dif")
  expect_s3_class(gh_tests, "gRm_global_homogeneity")
  expect_true(all(vapply(
    list(parameters, item_tests, ld_tests, dif_tests, gh_tests),
    function(x) is.data.frame(summary(x)$tests),
    logical(1L)
  )))
})

test_that("removed accessors are not available as public functions", {
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
})
