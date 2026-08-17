output_accessor_data <- function() {
  complete <- expand.grid(
    I1 = 0:2,
    I2 = 0:2,
    I3 = 0:2,
    site = 0:1,
    KEEP.OUT.ATTRS = FALSE
  )
  incomplete <- data.frame(
    I1 = c(0L, 1L),
    I2 = c(NA_integer_, 2L),
    I3 = c(1L, 2L),
    site = c(1L, NA_integer_)
  )
  data <- rbind(complete, incomplete)
  data$ID <- seq_len(nrow(data))
  data[c("ID", "I1", "I2", "I3", "site")]
}

test_that("summary is the public output accessor for analysis objects", {
  ia <- gRm(
    output_accessor_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID"
  )

  out <- summary(ia)

  expect_s3_class(out, "summary.gRm_analysis")
  expect_true("data" %in% names(out$tables))
  expect_equal(out$data$n_items, 3L)
  expect_false("tidy" %in% getNamespaceExports("gRm"))
  expect_false("glance" %in% getNamespaceExports("gRm"))
  expect_false("details" %in% getNamespaceExports("gRm"))
  expect_false("detail_names" %in% getNamespaceExports("gRm"))
})

test_that("summary table helpers do not expose obsolete details naming", {
  namespace_names <- ls(asNamespace("gRm"), all.names = TRUE)
  expect_false(any(namespace_names %in% c(
    "empty_detail",
    "bind_detail_rows",
    "detail_matrix_long",
    "detail_named_vector",
    "item_gamma_detail",
    "item_threshold_detail",
    "item_score_matrix_detail",
    "item_effect_detail",
    "item_parameter_detail_tables",
    "item_parameters",
    "summary.gRm_item_parameters",
    "print.gRm_item_parameters"
  )))
})

test_that("summary which exposes model, fit, screen, and diagnostic sections", {
  ia <- gRm(
    output_accessor_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID"
  )
  model <- gllrm(ia)
  fit_obj <- fit(model, max_step = 50L)
  screen_obj <- screen(ia, inference = "asymptotic")

  expect_true(is.data.frame(summary(model)$model))
  expect_true(is.data.frame(summary(fit_obj)$parameters))
  expect_true(is.data.frame(summary(fit_obj, which = "thresholds")$thresholds))
  expect_true(is.data.frame(summary(screen_obj)$local_dependence))
  expect_true(is.data.frame(summary(fit_obj, which = "parameters")$parameters))
  expect_true(is.data.frame(item_fit(fit_obj)))
})

test_that("item fit items summary is empty when extended item summaries are not computed", {
  ia <- gRm(
    output_accessor_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID"
  )
  fit_obj <- fit(gllrm(ia), max_step = 50L)

  tests <- item_fit(fit_obj, include_extended = FALSE)
  items <- item_fit(fit_obj, which = "items", include_extended = FALSE)

  expect_s3_class(tests, "gRm_item_fit")
  expect_s3_class(items, "gRm_item_fit")
  expect_true(is.data.frame(tests))
  expect_gt(nrow(tests), 0L)
  expect_equal(names(tests), c(
    "Item",
    "Outfit",
    "Outfit SE",
    "Pr(>Outfit)",
    "Outfit FDR",
    "Infit",
    "Infit SE",
    "Pr(>Infit)",
    "Infit FDR",
    "Observed gamma",
    "Expected gamma",
    "Gamma SE",
    "Pr(>Gamma)",
    "Gamma FDR",
    "Gamma direction"
  ))
  expect_true(is.data.frame(items))
  expect_equal(nrow(items), 0L)
  expect_equal(names(items), c(
    "Item",
    "Outfit N",
    "Outfit observed",
    "Outfit expected",
    "Outfit total",
    "Infit observed",
    "Infit expected",
    "Infit variance",
    "Infit ratio"
  ))
  expect_false(identical(names(items), names(tests)))
})

test_that("item fit result object retains computed values needed by validators", {
  ia <- gRm(
    output_accessor_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID"
  )
  values <- item_fits_values(ia$project, include_extended = TRUE)
  public <- item_fit(fit(gllrm(ia), max_step = 50L), include_extended = TRUE)

  expect_s3_class(values, "gRm_item_fits_values")
  expect_true(is.data.frame(values$items))
  expect_true(is.data.frame(values$side_file))
  expect_true(is.numeric(values$bh_limits))
  expect_true(is.list(values$extended))
  expect_true(is.data.frame(values$extended$score_n))
  expect_true(is.data.frame(values$extended$scores))
  expect_true(is.data.frame(values$extended$summaries))
  expect_true(is.data.frame(values$extended$restscore_tables))
  expect_true(is.data.frame(values$extended$local_restscore_tables))
  expect_true(is.data.frame(values$extended$local_gamma$rows))
  expect_true(is.data.frame(values$extended$local_gamma$limits))
  expect_s3_class(public, "gRm_item_fit")
  expect_s3_class(attr(public, "values", exact = TRUE), "gRm_item_fits_values")
  expect_true(is.data.frame(attr(public, "bh", exact = TRUE)))
  expect_true(all(c("fdr_5", "fdr_1") %in% attr(public, "bh", exact = TRUE)$threshold))
  expect_false("values" %in% names(public))
  expect_null(attr(public, "analysis", exact = TRUE))
  expect_null(attr(public, "fit", exact = TRUE))
  expect_true(all(c("item_label", "outfit_fdr", "infit_fdr", "gamma_fdr") %in% names(attr(public, "values", exact = TRUE)$items)))
  expect_true(all(c(
    "item_label", "item_name", "outfit_total_value", "infit_value"
  ) %in% names(attr(public, "values", exact = TRUE)$extended$summaries)))
})
