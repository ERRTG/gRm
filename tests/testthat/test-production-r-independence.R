test_that("production R files do not call Pascal reference helpers", {
  r_dir <- repo_path("gRm", "R")
  production_files <- setdiff(
    list.files(r_dir, pattern = "[.]R$", full.names = TRUE),
    file.path(r_dir, "pascal_reference.R")
  )

  offending_lines <- character()
  for (path in production_files) {
    lines <- readLines(path, warn = FALSE)
    non_comment_lines <- lines[!grepl("^\\s*#", lines)]
    hits <- grep(
      "run_pascal_[A-Za-z0-9_]*_reference\\s*\\(|system2\\s*\\(|SOURCE_GLLRM|LEGACY_[A-Z_]*REPORT",
      non_comment_lines,
      value = TRUE
    )
    if (length(hits) > 0L) {
      offending_lines <- c(
        offending_lines,
        paste0(basename(path), ": ", trimws(hits))
      )
    }
  }

  expect_equal(offending_lines, character())
})

test_that("production R executable paths do not branch on validation fixture names", {
  r_dir <- repo_path("gRm", "R")
  production_files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)

  fixture_regex <- paste(c(
    "legacy_disinh",
    "legacy_emo",
    "legacy_impuls",
    "legacy_movtiva",
    "legacy_soccog"
  ), collapse = "|")

  offending_lines <- character()
  for (path in production_files) {
    lines <- readLines(path, warn = FALSE)
    non_comment_lines <- lines[!grepl("^\\s*#", lines)]
    hits <- grep(fixture_regex, non_comment_lines, value = TRUE)
    if (length(hits) > 0L) {
      offending_lines <- c(
        offending_lines,
        paste0(basename(path), ": ", trimws(hits))
      )
    }
  }

  expect_equal(offending_lines, character(), info = paste(offending_lines, collapse = "\n"))
})
