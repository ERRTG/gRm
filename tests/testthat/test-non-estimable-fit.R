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

test_that("public fitting rejects upper-boundary and one-category support", {
  cases <- list(
    upper_boundary_only = list(
      data = data.frame(
        ID = seq_len(3L),
        I1 = rep(1L, 3L),
        I2 = rep(1L, 3L)
      ),
      nvalid = 0L,
      error = "at least one source-valid complete response pattern"
    ),
    one_category_per_item = list(
      data = data.frame(
        ID = seq_len(4L),
        I1 = rep(0L, 4L),
        I2 = rep(1L, 4L)
      ),
      nvalid = 4L,
      error = "supported response categories"
    ),
    mixed_valid_invalid = list(
      data = data.frame(
        ID = seq_len(6L),
        I1 = rep(c(1L, 0L), 3L),
        I2 = rep(c(0L, 1L), 3L),
        I3 = rep(0L, 6L)
      ),
      nvalid = 6L,
      error = "supported response categories"
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    case_data <- case$data
    item_names <- setdiff(names(case_data), "ID")
    levels <- stats::setNames(rep(list(0:1), length(item_names)), item_names)
    analysis <- gRm(
      case_data,
      items = item_names,
      id = "ID",
      item_levels = levels,
      score_cuts = c(1L, length(item_names))
    )
    bundle <- build_item_parameters_bundle(analysis$project)

    expect_identical(bundle$manifest$nvalid, case$nvalid, info = case_name)
    expect_error(
      fit(gllrm(analysis), max_step = 20L),
      case$error,
      info = case_name
    )
    expect_error(
      fit(gllrm(analysis, ld = stats::as.formula(paste0("~ ", item_names[[1L]], ":", item_names[[2L]]))), max_step = 20L),
      case$error,
      info = paste(case_name, "GLLRM")
    )
  }
})

test_that("public fitting rejects a one-item conditional model", {
  analysis <- gRm(
    data.frame(ID = seq_len(4L), I1 = c(0L, 1L, 0L, 1L)),
    items = "I1",
    id = "ID",
    item_levels = list(I1 = 0:1),
    score_cuts = c(0L, 1L)
  )

  expect_error(
    fit(gllrm(analysis), max_step = 20L),
    "at least one source-valid complete response pattern"
  )
})

test_that("sparse observed patterns do not create a false rank rejection", {
  data <- data.frame(
    ID = seq_len(8L),
    I1 = rep(c(1L, 0L), each = 4L),
    I2 = rep(c(1L, 0L), each = 4L),
    I3 = rep(c(0L, 1L), each = 4L),
    I4 = rep(c(0L, 1L), each = 4L)
  )
  analysis <- gRm(
    data,
    items = c("I1", "I2", "I3", "I4"),
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1, I4 = 0:1),
    score_cuts = c(2L, 4L)
  )

  fitted <- fit(gllrm(analysis), max_step = 20L)
  expect_s3_class(fitted, "gRm_fit")
  expect_true(fitted$convergence$converged)
})

test_that("interaction terms require positive source reference-free support", {
  data <- data.frame(
    ID = seq_len(9L),
    I1 = rep(c(0L, 0L, 1L), 3L),
    I2 = rep(c(0L, 1L, 0L), 3L),
    I3 = rep(c(1L, 0L, 0L), 3L)
  )
  analysis <- gRm(
    data,
    items = c("I1", "I2", "I3"),
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1),
    score_cuts = c(1L, 3L)
  )

  expect_error(
    fit(gllrm(analysis, ld = ~ I1:I2), max_step = 20L),
    "without estimable observed support"
  )
})

test_that("constant total score remains estimable with supported margins", {
  data <- data.frame(
    ID = seq_len(9L),
    I1 = rep(c(1L, 0L, 0L), 3L),
    I2 = rep(c(0L, 1L, 0L), 3L),
    I3 = rep(c(0L, 0L, 1L), 3L)
  )
  analysis <- gRm(
    data,
    items = c("I1", "I2", "I3"),
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1),
    score_cuts = c(1L, 3L)
  )

  fitted <- fit(gllrm(analysis), max_step = 50L)
  expect_s3_class(fitted, "gRm_fit")
  expect_equal(fitted$values$n_parameters, 2L)
  expect_false(any(fitted$values$thresholds <= -999999, na.rm = TRUE))
})

test_that("public fitting rejects holes in declared item-category support", {
  data <- data.frame(
    ID = seq_len(8L),
    I1 = rep(c(0L, 2L), 4L),
    I2 = rep(c(1L, 0L), 4L)
  )
  analysis <- gRm(
    data,
    items = c("I1", "I2"),
    id = "ID",
    item_levels = list(I1 = 0:2, I2 = 0:1),
    score_cuts = c(1L, 3L)
  )

  expect_error(
    fit(gllrm(analysis), max_step = 20L),
    "missing category support"
  )
})
