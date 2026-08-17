public_argument_objects <- function() {
  data <- expand.grid(
    I1 = 0:1,
    I2 = 0:1,
    I3 = 0:1,
    X = 0:1,
    KEEP.OUT.ATTRS = FALSE
  )
  data$ID <- seq_len(nrow(data))
  analysis <- gRm(
    data[c("ID", "I1", "I2", "I3", "X")],
    items = c("I1", "I2", "I3"),
    exogenous = "X",
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1),
    exogenous_levels = list(X = 0:1),
    score_cuts = c(1L, 3L)
  )
  model <- gllrm(analysis)
  fitted <- fit(model, max_step = 50L)
  list(data = data, analysis = analysis, model = model, fitted = fitted)
}

test_that("public workflows reject named unsupported dots and report their names", {
  objects <- public_argument_objects()
  fitted <- objects$fitted

  expect_error(fit(objects$model, max_stpe = 1L), "`max_stpe`")
  expect_error(screen(objects$analysis, n_sim = 1L), "`n_sim`")
  expect_error(ari(fitted, confdence = 0.9), "`confdence`")
  expect_error(item_fit(fitted, includ_extended = FALSE), "`includ_extended`")
  expect_error(local_dependence(fitted, max_stpe = 1L), "`max_stpe`")
  expect_error(dif(fitted, max_stpe = 1L), "`max_stpe`")
  expect_error(global_homogeneity(fitted, max_stpe = 1L), "`max_stpe`")
  expect_error(m2(fitted, itmes = "I1"), "`itmes`")
  expect_error(m3(fitted, itmes = "I1"), "`itmes`")
  expect_error(update(objects$model, ldd = ~ I1:I2), "`ldd`")
  expect_error(model_graph(objects$model, layuot = "fr"), "`layuot`")
  expect_error(logLik(fitted, na_rm = TRUE), "`na_rm`")
  expect_error(summary(fitted, detials = TRUE), "`detials`")
  expect_error(print(fitted, detials = TRUE), "`detials`")
})

test_that("exported functions without dots let R reject misspellings", {
  objects <- public_argument_objects()
  expect_error(
    gRm(objects$data, items = c("I1", "I2", "I3"), itemz = "I1"),
    "unused argument.*itemz"
  )
  expect_error(gllrm(objects$analysis, ldd = ~ I1:I2), "unused argument.*ldd")
  expect_error(score_effects(objects$analysis, n_sim = 1L), "unused argument.*n_sim")
})

test_that("item_fit include_extended is one strict logical control", {
  fitted <- public_argument_objects()$fitted
  for (bad in list(NA, c(TRUE, FALSE), 1, "TRUE", NULL)) {
    expect_error(
      item_fit(fitted, include_extended = bad),
      "`include_extended` must be a single non-missing logical value"
    )
  }

  included <- item_fit(fitted, include_extended = TRUE)
  omitted <- item_fit(fitted, include_extended = FALSE)
  expect_true(attr(included, "metadata", exact = TRUE)$include_extended)
  expect_false(attr(omitted, "metadata", exact = TRUE)$include_extended)
})

test_that("public plotting logical controls use the same scalar contract", {
  objects <- public_argument_objects()
  ari_table <- ari(objects$fitted)

  for (bad in list(NA, c(TRUE, FALSE), 1, "TRUE", NULL)) {
    expect_error(
      plot(ari_table, show_expected = bad),
      "`show_expected` must be a single non-missing logical value"
    )
    expect_error(
      plot(objects$model, rescale = bad),
      "`rescale` must be a single non-missing logical value"
    )
  }
  expect_error(plot(ari_table, class_sze = 4L), "`class_sze`")
})
