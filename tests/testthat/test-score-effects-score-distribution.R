test_that("score-effect score distribution handles no known item scores", {
  exo_select_score_distribution <- get(
    "exo_select_score_distribution",
    envir = asNamespace("gRm")
  )

  score_distribution <- exo_select_score_distribution(
    scores = c(-1L, -1L, -1L),
    known = c(FALSE, FALSE, FALSE),
    min_score = 0L,
    max_score = 2L,
    largest_possible_score = 2L
  )

  expect_equal(nrow(score_distribution$distribution), 0L)
  expect_equal(score_distribution$summary$n, 0L)
  expect_true(is.na(score_distribution$summary$mean))
  expect_true(is.na(score_distribution$summary$variance))
  expect_true(is.na(score_distribution$summary$sd))
  expect_true(is.na(score_distribution$summary$skewness))
  expect_equal(score_distribution$summary$below, 0L)
  expect_equal(score_distribution$summary$above, 0L)
  expect_equal(score_distribution$summary$missing, 3L)

  default_empty <- expect_silent(exo_select_score_distribution(integer()))
  expect_equal(default_empty$summary$n, 0L)
  expect_equal(default_empty$summary$missing, 0L)
})

test_that("score_effects score distribution requires known exogenous values", {
  data <- data.frame(
    ID = seq_len(6L),
    I1 = c(0L, 1L, 1L, NA, 0L, 1L),
    I2 = c(0L, 0L, 1L, 1L, 1L, 0L),
    site = c(0L, 1L, NA, 1L, 0L, 1L)
  )
  analysis <- gRm(
    data,
    items = c("I1", "I2"),
    exogenous = "site",
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1),
    exogenous_levels = list(site = 0:1),
    score_cuts = c(1L, 2L)
  )

  effects <- score_effects(analysis)
  values <- attr(effects, "values", exact = TRUE)

  expect_equal(values$score_summary$n, 4L)
  expect_equal(values$score_summary$missing, 0L)
  expect_equal(values$score_distribution$score, 0:1)
  expect_equal(values$score_distribution$count, c(1L, 3L))
})
