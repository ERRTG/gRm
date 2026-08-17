m2_m3_api_data <- function() {
  rows <- expand.grid(
    I1 = 1:3,
    I2 = 1:3,
    I3 = 1:3,
    site = 1:2,
    KEEP.OUT.ATTRS = FALSE
  )
  data.frame(ID = seq_len(nrow(rows)), rows)
}

m2_m3_api_analysis <- function(score_cuts = c(2L, 6L)) {
  gRm(
    m2_m3_api_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID",
    score_cuts = score_cuts
  )
}

m2_m3_api_fit <- function(score_cuts = c(2L, 6L)) {
  fit(gllrm(m2_m3_api_analysis(score_cuts = score_cuts)), max_step = 100L)
}

test_that("m2 and m3 expose source-named public diagnostics", {
  expect_true(exists("m2", envir = asNamespace("gRm"), inherits = FALSE))
  expect_true(exists("m3", envir = asNamespace("gRm"), inherits = FALSE))
  expect_false("which" %in% names(formals(gRm::m2)))
  expect_false("which" %in% names(formals(gRm::m3)))
  expect_true("items" %in% names(formals(gRm::m2)))
  expect_true("items" %in% names(formals(gRm::m3)))
  expect_true("score_cuts" %in% names(formals(gRm::m2)))
  expect_true("score_cuts" %in% names(formals(gRm::m3)))
})

test_that("m2 and m3 require fitted model objects", {
  analysis <- m2_m3_api_analysis()
  fit0 <- fit(gllrm(analysis), max_step = 100L)

  expect_error(m2(analysis), "fitted", fixed = TRUE)
  expect_error(m3(analysis), "fitted", fixed = TRUE)
  expect_s3_class(m2(fit0), "gRm_m2")
  expect_s3_class(m3(fit0), "gRm_m3")
})

test_that("m2 and m3 score cuts are diagnostic-local source state", {
  fit0 <- m2_m3_api_fit(score_cuts = c(2L, 6L))
  analysis_cuts <- fit0$analysis$score_groups
  explicit_cuts <- c(3L, 6L)

  m2_result <- m2(fit0, score_cuts = explicit_cuts)
  m3_result <- m3(fit0, score_cuts = explicit_cuts)

  expect_equal(m2_result$metadata$score_cuts, explicit_cuts)
  expect_equal(m3_result$metadata$score_cuts, explicit_cuts)
  expect_equal(m2_result$values$score_cuts, explicit_cuts)
  expect_equal(m3_result$values$score_cuts, explicit_cuts)
  expect_equal(fit0$analysis$score_groups, analysis_cuts)
})

test_that("m2 and m3 normalize selected items to source order", {
  fit0 <- m2_m3_api_fit()

  m2_result <- m2(fit0, items = c("I3", "I1"))
  m3_result <- m3(fit0, items = c(3L, 1L))

  expect_equal(m2_result$metadata$item_indices, c(1L, 3L))
  expect_equal(m3_result$metadata$item_indices, c(1L, 3L))
  expect_equal(m2_result$metadata$items, c("I1", "I3"))
  expect_equal(m3_result$metadata$items, c("I1", "I3"))
})

test_that("m2 and m3 print and summary use planned headers", {
  fit0 <- m2_m3_api_fit()
  m2_result <- m2(fit0)
  m3_result <- m3(fit0)

  expect_match(capture.output(print(m2_result))[[1L]], "^gRm: M2 fit diagnostic$")
  expect_match(capture.output(print(m3_result))[[1L]], "^gRm: M3 fit diagnostic$")
  expect_match(capture.output(summary(m2_result))[[1L]], "^gRm: M2 fit diagnostic$")
  expect_match(capture.output(summary(m3_result))[[1L]], "^gRm: M3 fit diagnostic$")
  expect_true(any(grepl("^  Score cuts: ", capture.output(print(m2_result)))))
  expect_true(any(grepl("^Score cuts: ", capture.output(summary(m3_result)))))
})

test_that("m2 and m3 print and summary use diagnostic numeric formatting", {
  old_digits <- getOption("digits")
  options(digits = 7L)
  on.exit(options(digits = old_digits), add = TRUE)

  stat <- data.frame(
    diagnostic = "M2",
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
    m2 = stat,
    m3 = transform(stat, diagnostic = "M3"),
    item_trait = transform(stat, diagnostic = "Item-trait"),
    invariance = data.frame(),
    tests = transform(
      stat,
      margin = "I1 x I2",
      margin_type = "item_item",
      variable1 = "I1",
      variable2 = "I2",
      variable3 = NA_character_,
      is_m2 = TRUE
    )
  )
  result <- structure(list(values = values), class = "gRm_m2")

  printed <- capture.output(print(result))
  expect_true(any(grepl("^  M2 Chisq: 12.3$", printed)))
  expect_false(any(grepl("12.3457", printed, fixed = TRUE)))

  summary_printed <- capture.output(print(summary(result)))
  expect_true(any(grepl("12.3", summary_printed, fixed = TRUE)))
  expect_false(any(grepl("12.3457", summary_printed, fixed = TRUE)))
  expect_true(any(grepl("0.000123", summary_printed, fixed = TRUE)))
})
