test_that("gRm package scaffold is present", {
  root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  if (!dir.exists(file.path(root, "R"))) {
    root <- system.file(package = "gRm")
  }

  expect_true(file.exists(file.path(root, "DESCRIPTION")))
  expect_true(file.exists(file.path(root, "NAMESPACE")))
  expect_true(dir.exists(file.path(root, "R")))
  expect_true(file.exists(file.path(root, "R", "gRm-package.R")))
  expect_true(file.exists(file.path(root, "R", "api-constructors.R")))
  expect_true(file.exists(file.path(root, "R", "api-model-spec.R")))
  expect_true(file.exists(file.path(root, "R", "api-summary.R")))
  expect_true(file.exists(file.path(root, "tests", "testthat.R")))
})
