test_that("public fitting rejects analyses with no source-valid rows", {
  cases <- list(
    all_missing = data.frame(
      ID = seq_len(3L),
      I1 = c(NA_integer_, NA_integer_, NA_integer_),
      I2 = c(NA_integer_, NA_integer_, NA_integer_)
    ),
    no_complete_item_rows = data.frame(
      ID = seq_len(3L),
      I1 = c(0L, NA_integer_, 1L),
      I2 = c(NA_integer_, 1L, NA_integer_)
    ),
    lower_boundary_only = data.frame(
      ID = seq_len(3L),
      I1 = c(0L, 0L, 0L),
      I2 = c(0L, 0L, 0L)
    )
  )

  for (case_data in cases) {
    analysis <- gRm(
      case_data,
      items = c("I1", "I2"),
      id = "ID",
      item_levels = list(I1 = 0:1, I2 = 0:1),
      score_cuts = c(1L, 2L)
    )

    bundle <- build_item_parameters_bundle(analysis$project)
    expect_equal(bundle$manifest$nvalid, 0L)
    expect_error(
      fit(gllrm(analysis), max_step = 20L),
      "at least one source-valid complete response pattern"
    )
    expect_error(
      fit(gllrm(analysis, ld = ~ I1:I2), max_step = 20L),
      "at least one source-valid complete response pattern"
    )
  }
})

test_that("Rasch fitting still accepts source-valid complete rows", {
  analysis <- gRm(
    data.frame(
      ID = seq_len(4L),
      I1 = c(0L, 0L, 1L, 1L),
      I2 = c(0L, 1L, 0L, 1L)
    ),
    items = c("I1", "I2"),
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1),
    score_cuts = c(1L, 2L)
  )

  bundle <- build_item_parameters_bundle(analysis$project)
  expect_gt(bundle$manifest$nvalid, 0L)
  expect_s3_class(fit(gllrm(analysis), max_step = 20L), "gRm_fit")
})
