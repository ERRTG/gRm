#!/usr/bin/env Rscript

# Regenerate the package function inventory and relocate source-trace rows after
# responsibility-preserving file splits. The curated numerical manifest remains
# authoritative for deciding which helpers implement numerical source logic.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) {
  stop("Run this file with Rscript.", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
package_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

top_level_functions <- function(path) {
  parsed <- parse(path, keep.source = TRUE)
  references <- attr(parsed, "srcref")
  rows <- list()
  for (index in seq_along(parsed)) {
    expression <- parsed[[index]]
    if (
      is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) && is.call(expression[[3L]]) &&
      identical(expression[[3L]][[1L]], as.name("function"))
    ) {
      formal_names <- names(as.list(expression[[3L]][[2L]]))
      if (is.null(formal_names)) {
        formal_names <- character()
      }
      rows[[length(rows) + 1L]] <- list(
        name = as.character(expression[[2L]]),
        line = as.integer(references[[index]][1L]),
        formals = formal_names
      )
    }
  }
  rows
}

roxygen_block <- function(lines, function_line) {
  if (function_line <= 1L || !grepl("^#'", lines[[function_line - 1L]])) {
    return(character())
  }
  first <- function_line - 1L
  while (first > 1L && grepl("^#'", lines[[first - 1L]])) {
    first <- first - 1L
  }
  lines[first:(function_line - 1L)]
}

namespace <- readLines(file.path(package_root, "NAMESPACE"), warn = FALSE)
exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", namespace, value = TRUE))
methods <- sub(
  "^S3method\\(([^,]+),([^\\)]+)\\)$",
  "\\1.\\2",
  grep("^S3method\\(", namespace, value = TRUE)
)

manifest_path <- file.path(package_root, "tools", "numerical-source-trace.csv")
manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
if (anyDuplicated(manifest$r_function)) {
  stop("The numerical source-trace manifest has duplicate function names.", call. = FALSE)
}

r_files <- sort(list.files(file.path(package_root, "R"), pattern = "[.]R$", full.names = TRUE))
observed <- list()
for (path in r_files) {
  lines <- readLines(path, warn = FALSE)
  entries <- top_level_functions(path)
  for (entry in entries) {
    block <- roxygen_block(lines, entry$line)
    observed[[length(observed) + 1L]] <- data.frame(
      r_function = entry$name,
      r_file = file.path("R", basename(path)),
      line = entry$line,
      exported = entry$name %in% exports,
      s3_method = entry$name %in% methods,
      numerical = entry$name %in% manifest$r_function,
      roxygen_attached = length(block) > 0L,
      keywords_internal = any(grepl("@keywords\\s+internal(?:\\s|$)", block, perl = TRUE)),
      formals = paste(entry$formals, collapse = ";"),
      source_trace_manifest = entry$name %in% manifest$r_function,
      stringsAsFactors = FALSE
    )
  }
}
inventory <- do.call(rbind, observed)
if (anyDuplicated(inventory$r_function)) {
  duplicates <- unique(inventory$r_function[duplicated(inventory$r_function)])
  stop("Duplicate top-level function names: ", paste(duplicates, collapse = ", "), call. = FALSE)
}
missing_manifest <- setdiff(manifest$r_function, inventory$r_function)
if (length(missing_manifest)) {
  stop(
    "Manifest functions missing from R/: ", paste(missing_manifest, collapse = ", "),
    call. = FALSE
  )
}

inventory <- inventory[order(inventory$r_file, inventory$line), , drop = FALSE]
rownames(inventory) <- NULL
utils::write.csv(
  inventory,
  file.path(package_root, "tools", "r-function-inventory.csv"),
  row.names = FALSE,
  na = ""
)

location <- setNames(file.path("gRm", inventory$r_file), inventory$r_function)
manifest$r_file <- unname(location[manifest$r_function])
manifest <- manifest[order(manifest$r_file, manifest$r_function), , drop = FALSE]
rownames(manifest) <- NULL
utils::write.csv(manifest, manifest_path, row.names = FALSE, na = "")

cat(
  sprintf(
    "Wrote %d function rows and relocated %d numerical source-trace rows.\n",
    nrow(inventory), nrow(manifest)
  )
)
