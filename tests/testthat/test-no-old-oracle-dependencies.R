test_that("old oracle parity tests have been retired", {
  retired_tests <- c(
    "test-bfi-raw-csv-project.R",
    "test-dif-tests-values-vs-pascal.R",
    "test-exo-select-values-vs-pascal.R",
    "test-gamma-values-vs-gRm.R",
    "test-gamma-values-vs-pascal.R",
    "test-global-homogeneity-values-vs-gRm.R",
    "test-global-homogeneity-vs-pascal-oracle.R",
    "test-item-fits-extended-values-vs-gRm.R",
    "test-item-fits-extended-values-vs-pascal.R",
    "test-item-fits-values-vs-pascal.R",
    "test-item-parameter-values-vs-pascal.R",
    "test-items-select-values-vs-gRm.R",
    "test-local-independence-values-vs-pascal.R",
    "test-oracle-gap-manifest.R",
    "test-project-interface.R",
    "test-rasch-base-fit-vs-pascal.R",
    "test-read-gRm-project.R",
    "test-render-dif-tests-vs-gRm.R",
    "test-render-exo-select-vs-gRm.R",
    "test-render-gamma-vs-gRm.R",
    "test-render-global-homogeneity-vs-gRm.R",
    "test-render-global-invariance-vs-gRm.R",
    "test-render-item-fits-extended-vs-gRm.R",
    "test-render-item-fits-vs-gRm.R",
    "test-render-item-parameters-vs-gRm.R",
    "test-render-items-select-vs-gRm.R",
    "test-render-local-independence-vs-gRm.R",
    "test-render-local-invariance-vs-bfi.R",
    "test-render-screen-j-vs-legacy-fixture.R",
    "test-render-screen-vs-gRm.R",
    "test-screen-values-vs-pascal.R"
  )
  retired_paths <- repo_path("gRm", "tests", "testthat", retired_tests)

  expect_false(any(file.exists(retired_paths)))
})

test_that("report reconstruction tests no longer validate against old data", {
  test_files <- list.files(
    repo_path("gRm", "tests", "testthat"),
    pattern = "[.]R$",
    full.names = TRUE
  )
  test_files <- test_files[basename(test_files) != "test-no-old-oracle-dependencies.R"]
  reconstruction_files <- grepl(
    "-vs-|oracle|numeric-reconstruction|report-reconstruction|^test-render-|values",
    basename(test_files)
  )
  test_files <- test_files[reconstruction_files]

  old_oracle_patterns <- c(
    "data/legacy_analysis",
    "data/bfi_analysis",
    "repo_path\\(\"data\"",
    "file\\.path\\([^\\n]*(repo_root\\(\\)|root)[^\\n]*,\\s*\"data\"",
    "read_digram_project\\(",
    "read_raw_csv_project\\(",
    "screen_J\\.txt",
    "screen_J_exact\\.txt",
    "ItemParameters\\.txt",
    "ItemFits\\.txt",
    "DIF-tests\\.txt"
  )
  old_oracle_regex <- paste(old_oracle_patterns, collapse = "|")

  old_oracle_files <- vapply(test_files, function(path) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    grepl(old_oracle_regex, text)
  }, logical(1))

  offender_names <- basename(test_files[old_oracle_files])
  expect_equal(
    offender_names,
    character(),
    info = paste(offender_names, collapse = "\n")
  )
})

test_that("production R code does not read examples output reports", {
  r_files <- list.files(repo_path("gRm", "R"), pattern = "[.]R$", full.names = TRUE)
  r_files <- r_files[basename(r_files) != "examples_validation.R"]
  text <- paste(unlist(lapply(r_files, readLines, warn = FALSE)), collapse = "\n")

  expect_false(grepl("examples/.*/output", text))
  expect_false(grepl("readLines\\([^\\n]*(check-|exo|item-|items|screen).*\\.txt", text))
})

test_that("production R executable paths do not depend on old data or output files", {
  r_dir <- repo_path("gRm", "R")
  validation_only <- c("examples_manifest.R", "examples_validation.R", "pascal_reference.R")
  r_files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  r_files <- r_files[!basename(r_files) %in% validation_only]

  old_oracle_regex <- paste(c(
    "data/legacy_analysis",
    "data/bfi_analysis",
    "read_digram_project\\(",
    "read_raw_csv_project\\(",
    "screen_J\\.txt",
    "screen_J_exact\\.txt",
    "ItemParameters\\.txt",
    "ItemFits\\.txt",
    "DIF-tests\\.txt"
  ), collapse = "|")

  offenders <- character()
  for (path in r_files) {
    lines <- readLines(path, warn = FALSE)
    executable_lines <- lines[!grepl("^\\s*#", lines)]
    hits <- grep(old_oracle_regex, executable_lines, value = TRUE)
    if (length(hits) > 0L) {
      offenders <- c(offenders, paste0(basename(path), ": ", trimws(hits)))
    }
  }

  expect_equal(offenders, character(), info = paste(offenders, collapse = "\n"))
})
