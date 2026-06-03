test_that("R DIGRAM package scaffold is present", {
  root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  if (!dir.exists(file.path(root, "R"))) {
    root <- system.file(package = "gRm")
  }

  expect_true(file.exists(file.path(root, "DESCRIPTION")))
  expect_true(dir.exists(file.path(root, "R")))
  expect_true(
    dir.exists(file.path(root, "inst", "extdata")) ||
      dir.exists(file.path(root, "extdata"))
  )
})
