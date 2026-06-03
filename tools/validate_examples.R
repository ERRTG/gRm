#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0L) {
  sub("^--file=", "", script_arg[[1L]])
} else {
  "gRm/tools/validate_examples.R"
}
repo <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
setwd(repo)

pkgload::load_all("gRm", quiet = TRUE)

files <- character()
jobs <- 128L
if (length(args) > 0L) {
  for (arg in args) {
    if (grepl("^--jobs=", arg)) {
      jobs <- as.integer(sub("^--jobs=", "", arg))
    } else if (nzchar(arg)) {
      files <- c(files, strsplit(arg, ",", fixed = TRUE)[[1L]])
    }
  }
}

manifest <- gRm_examples_manifest("examples")
matrix <- gRm_examples_validation_matrix(manifest)
if (length(files) > 0L) {
  matrix <- matrix[matrix$file %in% files | matrix$report_id %in% files, , drop = FALSE]
}

summaries <- list()
failures <- list()

cat("project_id\treport_id\tfile\tgenerated_tokens\texpected_tokens\tcompared_tokens\tmismatches\tfirst_bad_index\tfirst_bad_generated\tfirst_bad_expected\tfirst_bad_generated_line\tfirst_bad_expected_line\n")
flush.console()

for (project_id in unique(matrix$project_id)) {
  project_row <- manifest[manifest$project_id == project_id, , drop = FALSE]
  project <- read_gRm_example_project(project_row)
  project_rows <- matrix[matrix$project_id == project_id, , drop = FALSE]

  for (row_index in seq_len(nrow(project_rows))) {
    row <- project_rows[row_index, , drop = FALSE]
    generated <- gRm_text_lines(gRm_run_examples_report(
      project = project,
      report_id = row$report_id,
      inference = row$inference,
      extended = row$extended,
      nsim = row$nsim,
      seed = row$seed,
      jobs = jobs
    )$text)
    comparison <- gRm_validate_example_numeric_tokens(generated, row$expected_path[[1L]])
    summary <- gRm_summarize_example_numeric_tokens(
      comparison,
      project_id = row$project_id[[1L]],
      report_id = row$report_id[[1L]],
      file = row$file[[1L]]
    )
    summaries[[length(summaries) + 1L]] <- summary
    bad <- comparison[!comparison$match, , drop = FALSE]
    if (nrow(bad) > 0L) {
      failures[[paste(row$project_id[[1L]], row$file[[1L]], sep = "/")]] <- bad
    }
    cat(paste(summary[1, c(
      "project_id",
      "report_id",
      "file",
      "generated_tokens",
      "expected_tokens",
      "compared_tokens",
      "mismatches",
      "first_bad_index",
      "first_bad_generated",
      "first_bad_expected",
      "first_bad_generated_line",
      "first_bad_expected_line"
    )], collapse = "\t"), "\n", sep = "")
    flush.console()
  }
}

result <- list(summary = do.call(rbind, summaries), failures = failures)
bad_summary <- result$summary[result$summary$mismatches > 0L, , drop = FALSE]

cat("\n")
cat("reports=", nrow(result$summary), "\n", sep = "")
cat("bad_reports=", nrow(bad_summary), "\n", sep = "")
cat("bad_tokens=", sum(bad_summary$mismatches), "\n", sep = "")

if (nrow(bad_summary) > 0L) {
  cat("\n")
  cat(format_examples_numeric_validation_failures(result, max_reports = nrow(bad_summary), max_tokens = 8L), sep = "\n")
  quit(status = 1L)
}
