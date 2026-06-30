test_that("gRm package scaffold is present", {
  root <- package_root()

  expect_true(file.exists(file.path(root, "DESCRIPTION")))
  expect_true(file.exists(file.path(root, "NAMESPACE")))
  expect_true(dir.exists(file.path(root, "R")))
  expect_true(file.exists(file.path(root, "R", "gRm-package.R")))
  expect_true(file.exists(file.path(root, "R", "api-constructors.R")))
  expect_true(file.exists(file.path(root, "R", "api-model-spec.R")))
  expect_true(file.exists(file.path(root, "R", "api-summary.R")))
  expect_true(file.exists(file.path(root, "tests", "testthat.R")))
})

test_that("compiled build products are cleaned from the package source", {
  root <- package_root()

  cleanup <- file.path(root, "cleanup")
  install_hook <- file.path(root, "src", "install.libs.R")
  expect_true(file.exists(cleanup))
  expect_true(file.exists(install_hook))
  if (file.exists(install_hook)) {
    hook_text <- paste(readLines(install_hook, warn = FALSE), collapse = "\n")
    expect_match(hook_text, "R_PACKAGE_DIR", fixed = TRUE)
    expect_match(hook_text, "unlink", fixed = TRUE)
  }
  if (file.exists(cleanup)) {
    cleanup_text <- paste(readLines(cleanup, warn = FALSE), collapse = "\n")
    expect_match(cleanup_text, "src/[*][.]o", fixed = FALSE)
    expect_match(cleanup_text, "src/[*][.]so", fixed = FALSE)
    expect_match(cleanup_text, "src/symbols.rds", fixed = TRUE)
    status <- system2(cleanup, stdout = TRUE, stderr = TRUE)
    status_code <- attr(status, "status")
    if (is.null(status_code)) {
      status_code <- 0L
    }
    expect_equal(status_code, 0L)
  }

  build_products <- list.files(
    file.path(root, "src"),
    pattern = "[.](o|so|dll|dylib)$|^symbols[.]rds$",
    full.names = FALSE
  )
  expect_equal(build_products, character())
})
