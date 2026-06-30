incomplete_item_fit_data <- function(include_incomplete = TRUE) {
  complete <- data.frame(
    ID = seq_len(8L),
    I1 = c(0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L),
    I2 = c(1L, 0L, 1L, 0L, 1L, 1L, 0L, 1L),
    I3 = c(0L, 0L, 1L, 1L, 1L, 0L, 1L, 0L),
    site = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  )
  if (!isTRUE(include_incomplete)) {
    return(complete)
  }
  rbind(
    complete,
    data.frame(ID = 9L, I1 = 1L, I2 = 0L, I3 = NA_integer_, site = 1L)
  )
}

incomplete_item_fit_analysis <- function(data) {
  gRm(
    data,
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1),
    exogenous_levels = list(site = 0:1)
  )
}

empty_incomplete_item_fit_records <- function(bundle) {
  records <- as.data.frame(
    stats::setNames(
      replicate(nrow(bundle$model$items), integer(), simplify = FALSE),
      bundle$model$items$name
    ),
    stringsAsFactors = FALSE
  )
  list(
    records = records,
    count = integer(),
    score = integer(),
    max_score = integer()
  )
}

test_that("item fits do not synthesize incomplete response records", {
  incomplete <- incomplete_item_fit_analysis(incomplete_item_fit_data(TRUE))

  bundle <- build_item_parameters_bundle(incomplete$project)
  incomplete_records <- collect_source_incomplete_records(bundle)
  expect_equal(nrow(incomplete_records$records), 1L)
  expect_equal(incomplete_records$count, 1L)
  expect_equal(incomplete_records$score, 1L)
  expect_equal(incomplete_records$max_score, 2L)

  fit <- fit_rasch_base(bundle, max_step = 50L)
  conditional <- item_conditional_moments(bundle, fit$item_gamma, include_probabilities = TRUE)
  empty_records <- empty_incomplete_item_fit_records(bundle)
  expected_fit <- calculate_conditional_item_fit_values(
    bundle,
    fit,
    conditional = conditional,
    incomplete = empty_records
  )
  expected_gamma <- calculate_item_restscore_gamma_values(
    bundle,
    fit,
    conditional = conditional,
    incomplete = empty_records
  )

  values <- item_fits_values(incomplete$project, max_step = 50L, include_extended = FALSE)

  expect_false(values$incomplete_records_used)
  expect_equal(values$incomplete_records_status, "not_source_backed")
  expect_equal(values$items$outfit, expected_fit$outfit)
  expect_equal(values$items$infit, expected_fit$infit)
  expect_equal(values$items$observed_gamma, expected_gamma$observed_gamma)
  expect_equal(values$items$expected_gamma, expected_gamma$expected_gamma)
  expect_equal(values$items$gamma_sd, expected_gamma$gamma_sd)
})

test_that("item fit incomplete-record restriction is documented", {
  item_fit_docs <- readLines(repo_path("gRm", "R", "api-results.R"), warn = FALSE)
  package_docs <- readLines(repo_path("gRm", "R", "gRm-package.R"), warn = FALSE)

  expect_true(any(grepl("Incomplete item fit records", item_fit_docs, fixed = TRUE)))
  expect_true(any(grepl("NincompleteRecs", item_fit_docs, fixed = TRUE)))
  expect_true(any(grepl("incomplete_records_used = FALSE", item_fit_docs, fixed = TRUE)))
  expect_true(any(grepl("incomplete_records_status = \"not_source_backed\"", item_fit_docs, fixed = TRUE)))
  expect_true(any(grepl("Incomplete item fit records", package_docs, fixed = TRUE)))
  expect_true(any(grepl("not synthesized", package_docs, fixed = TRUE)))
})
