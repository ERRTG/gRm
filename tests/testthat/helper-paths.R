repo_root <- function() {
  normalizePath(file.path(test_path(), "..", "..", ".."), mustWork = TRUE)
}

repo_path <- function(...) {
  file.path(repo_root(), ...)
}

source_gRm <- function() {
  r_dir <- repo_path("gRm", "R")
  files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  for (file in sort(files)) {
    sys.source(file, envir = globalenv())
  }
}

source_gRm()

normalize_report_text <- function(text) {
  text <- paste(text, collapse = "\n")
  text <- gsub("\r\n?", "\n", text)
  text <- gsub("-0\\.000", "0.000", text)
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  lines <- gsub("[[:space:]]+", " ", trimws(lines))
  lines[lines != ""]
}
