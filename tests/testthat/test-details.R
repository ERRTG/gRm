result_output_data <- function() {
  data.frame(
    ID = seq_len(12L),
    I1 = c(1L, 1L, 2L, 2L, 3L, 3L, 1L, 2L, 3L, 1L, 2L, 3L),
    I2 = c(1L, 2L, 1L, 3L, 2L, 3L, 3L, 1L, 2L, 2L, 3L, 1L),
    I3 = c(2L, 1L, 3L, 1L, 2L, 3L, 2L, 1L, 3L, 3L, 2L, 1L),
    site = c(1L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 2L)
  )
}

test_that("analysis summaries expose input metadata without public details", {
  ia <- gRm(
    result_output_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID"
  )

  out <- summary(ia)

  expect_s3_class(out, "summary.gRm_analysis")
  expect_true("data" %in% names(out$tables))
  expect_equal(out$data$n_items, 3L)
  expect_equal(out$data$n_exogenous, 1L)
  expect_error(summary(ia, which = "not_a_table"), "which")
})

test_that("screen summaries expose selected model terms and BH metadata", {
  ia <- gRm(result_output_data(), items = c("I1", "I2"), exogenous = "site", id = "ID")
  scr <- screen(ia, inference = "asymptotic")

  all_terms <- summary(scr, which = "all")$all
  bh <- summary(scr, which = "bh")$bh
  tests <- summary(scr, which = "tests")$tests

  expect_true(is.data.frame(all_terms))
  expect_true(is.data.frame(bh))
  expect_true(is.data.frame(tests))
  expect_true("Benjamini-Hochberg" %in% tests$selection)
})

test_that("diagnostic summaries expose result tables through summary which", {
  ia <- gRm(result_output_data(), items = c("I1", "I2"), exogenous = "site", id = "ID")
  fitted <- fit(gllrm(ia), max_step = 50L)

  ld <- local_dependence(fitted, jobs = 1L)
  dif_result <- dif(fitted, jobs = 1L)

  expect_true(is.data.frame(summary(ld, which = "tests")$tests))
  expect_true(is.data.frame(summary(dif_result, which = "tests")$tests))
  expect_error(summary(ld, which = "details"), "which")
})
