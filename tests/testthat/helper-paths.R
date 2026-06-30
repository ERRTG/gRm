repo_root <- function() {
  dirname(package_root())
}

repo_path <- function(...) {
  file.path(repo_root(), ...)
}

is_package_source_root <- function(path) {
  file.exists(file.path(path, "DESCRIPTION")) &&
    file.exists(file.path(path, "R", "api-summary.R"))
}

package_root <- local({
  root <- NULL

  function() {
    if (!is.null(root)) {
      return(root)
    }

    current_test_file <- getOption("testthat_path", NA_character_)
    if (is.na(current_test_file)) {
      stop("Could not locate the gRm package source root.", call. = FALSE)
    }

    source_root <- normalizePath(file.path(dirname(current_test_file), "..", ".."), mustWork = FALSE)
    if (!is_package_source_root(source_root)) {
      if (grepl("[.]Rcheck$", basename(source_root))) {
        r_check_roots <- normalizePath(
          c(
            file.path(source_root, "00_pkg_src", "gRm"),
            file.path(source_root, "..", "gRm")
          ),
          mustWork = FALSE
        )
        source_matches <- vapply(r_check_roots, is_package_source_root, logical(1))
        if (any(source_matches)) {
          source_root <- r_check_roots[which(source_matches)[[1L]]]
        }
      }
    }

    if (is_package_source_root(source_root)) {
      root <<- normalizePath(source_root, mustWork = TRUE)
      return(root)
    }

    stop("Could not locate the gRm package source root.", call. = FALSE)
  }
})

package_path <- function(...) {
  file.path(package_root(), ...)
}

source_gRm <- function() {
  r_dir <- package_path("R")
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
