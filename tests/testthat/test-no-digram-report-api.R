test_that("package does not export DIGRAM report or external validation APIs", {
  namespace_lines <- readLines(repo_path("gRm", "NAMESPACE"), warn = FALSE)
  old_report <- function(prefix = "", suffix = "") paste0(prefix, "gRm", "_", "report", suffix)
  forbidden_exports <- c(
    old_report("available_"),
    old_report(),
    old_report(suffix = "s"),
    old_report(suffix = "_registry"),
    old_report("write_"),
    old_report("write_", "s"),
    "validate",
    paste0("validate", "_examples")
  )
  forbidden_exports <- paste0("export(", forbidden_exports, ")")

  expect_equal(intersect(namespace_lines, forbidden_exports), character())
})

test_that("package R code contains no report generation layer", {
  r_dir <- test_path("..", "..", "R")
  r_files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  text_by_file <- lapply(r_files, readLines, warn = FALSE)
  names(text_by_file) <- basename(r_files)

  forbidden_files <- c(
    "api-reports.R",
    "report_functions.R",
    "report_registry.R",
    "report_spec.R",
    "run_reports.R",
    "render_report_set.R",
    "write_reports.R",
    "api-validation.R",
    "examples_validation.R",
    "pascal_reference.R"
  )
  expect_false(any(file.exists(file.path(r_dir, forbidden_files))))

  forbidden_patterns <- paste(c(
    paste0("new_", "gRm", "_", "report"),
    paste0("gRm", "_", "report", "_set"),
    paste0("gRm", "_", "report", "\\s*<-"),
    paste0("gRm", "_", "report", "s\\s*<-"),
    "validate\\s*<-",
    paste0("write_", "gRm", "_", "report"),
    paste0("write_", "gRm", "_", "report", "s"),
    paste0("validate", "_examples"),
    paste0("gRm", "_examples", "_validation"),
    paste0("gRm", "_examples", "_validation", "_matrix"),
    paste0("gRm", "_validate", "_example", "_numeric", "_tokens"),
    paste0("oracle", " report"),
    paste0("validation", "_corpus", "\\s*=")
  ), collapse = "|")

  offenders <- character()
  for (file in names(text_by_file)) {
    lines <- text_by_file[[file]]
    hits <- grep(forbidden_patterns, lines, value = TRUE)
    hits <- hits[!grepl("^\\s*#", hits)]
    if (length(hits)) {
      offenders <- c(offenders, paste0(file, ": ", trimws(hits)))
    }
  }

  expect_equal(offenders, character(), info = paste(offenders, collapse = "\n"))
})
