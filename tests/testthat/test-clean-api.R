test_that("namespace exports only the statistical modeling API", {
  exports <- getNamespaceExports("gRm")
  expected <- c(
    "gRm",
    "read_digram_project",
    "gllrm",
    "fit",
    "screen",
    "score_effects",
    "item_parameters",
    "item_fit",
    "local_dependence",
    "dif",
    "global_homogeneity"
  )
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

  expect_setequal(exports, expected)
  expect_false(any(removed %in% exports))
})

test_that("internal null-coalescing helper is package-owned and not exported", {
  expect_true(exists("%||%", envir = asNamespace("gRm"), inherits = FALSE))
  expect_false("%||%" %in% getNamespaceExports("gRm"))
})
