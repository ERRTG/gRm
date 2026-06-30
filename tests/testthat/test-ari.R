ari_boundary_fit <- function() {
  data <- data.frame(
    ID = seq_len(12L),
    I1 = c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L),
    I2 = c(0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 1L),
    I3 = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 0L)
  )
  fit(gllrm(gRm(data, items = c("I1", "I2", "I3"), id = "ID", score_cuts = c(1L, 3L))), max_step = 50L)
}

ari_mixed_fit <- function() {
  data <- data.frame(
    ID = seq_len(12L),
    B = c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L),
    T = c(0L, 0L, 1L, 1L, 2L, 2L, 0L, 1L, 2L, 0L, 1L, 2L)
  )
  fit(gllrm(gRm(
    data,
    items = c("B", "T"),
    id = "ID",
    item_levels = list(B = 0:1, T = 0:2),
    score_cuts = c(1L, 3L)
  )), max_step = 50L)
}

test_that("ari rejects non-fit objects", {
  expect_error(ari(gRm(data.frame(I1 = c(0L, 1L), I2 = c(1L, 0L)), items = c("I1", "I2"), score_cuts = c(1L, 2L))), "fitted DIGRAM model")
  expect_error(ari(list()), "fitted DIGRAM model")
})

test_that("ari returns the source-shaped ARI data frame", {
  out <- ari(ari_boundary_fit())

  expect_s3_class(out, "gRm_ari")
  expect_s3_class(out, "data.frame")
  expect_equal(
    names(out),
    c(
      "ItemNo", "Item", "Score", "n",
      "Obs0", "Obs1",
      "ObsMean", "ObsVar",
      "Exp0", "Exp1",
      "ExpMean", "ExpVar", "z"
    )
  )
  expect_gt(nrow(out), 0L)
  expect_true(is.numeric(out$ExpMean))
  expect_false(any(format(out$ExpMean, digits = 16L) == format(round(out$ExpMean, 4L), digits = 16L)))
})

test_that("ari excludes source score boundary rows", {
  out <- ari(ari_boundary_fit())

  expect_false(0L %in% out$Score)
  expect_false(3L %in% out$Score)
  expect_equal(sort(unique(out$Score)), c(1L, 2L))
})

test_that("ari uses row counts when the source-shaped bundle has counts", {
  fitted <- ari_boundary_fit()
  fitted$bundle$data$count <- ifelse(fitted$bundle$data$score == 1L, 2L, 1L)

  out <- ari(fitted)

  expect_equal(unique(out$n[out$Score == 1L]), 8)
  expect_equal(unique(out$n[out$Score == 2L]), 4)
})

test_that("ari keeps global category columns for mixed item maxima", {
  out <- ari(ari_mixed_fit())

  expect_true(all(c("Obs0", "Obs1", "Obs2", "Exp0", "Exp1", "Exp2") %in% names(out)))
  binary_rows <- out$Item == "B"
  expect_true(all(out$Obs2[binary_rows] == 0))
  expect_true(all(out$Exp2[binary_rows] == 0))
  expect_true(any(out$Obs2[out$Item == "T"] > 0))
})

test_that("print.gRm_ari is compact and invisible", {
  out <- ari(ari_boundary_fit())

  expect_output(returned <- print(out), "gRm: ARI item score curves")
  expect_identical(returned, out)
})
