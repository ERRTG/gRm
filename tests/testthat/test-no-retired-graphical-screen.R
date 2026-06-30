test_that("retired graphical SCREEN implementation is absent", {
  r_dir <- repo_path("gRm", "R")
  man_dir <- repo_path("gRm", "man")
  graphical_screen_r <- file.path(r_dir, "screen_values.R")
  graphical_screen_docs <- file.path(
    man_dir,
    c(
      "screen_values.Rd",
      "screen_two_way_matrix.Rd",
      "screen_hidden_matrix.Rd",
      "screen_final_model.Rd",
      "screen_partial_gamma_matrix.Rd"
    )
  )

  expect_false(file.exists(graphical_screen_r))
  expect_false(any(file.exists(graphical_screen_docs)))

  r_files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  executable_lines <- unlist(lapply(r_files, function(path) {
    lines <- readLines(path, warn = FALSE)
    lines[!grepl("^\\s*#", lines)]
  }), use.names = FALSE)

  expect_false(any(grepl("Rcpp::sourceCpp\\s*\\(", executable_lines)))
  expect_false(any(grepl("requireNamespace\\(\"Rcpp\"", executable_lines, fixed = FALSE)))

  description <- readLines(repo_path("gRm", "DESCRIPTION"), warn = FALSE)
  expect_false(any(grepl("^\\s*Rcpp\\s*,?\\s*$", description)))
})
