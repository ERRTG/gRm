#' DIGRAM import-project input and output helpers
#'
#' These helpers read or write legacy DIGRAM CSV/IMP/IMV import bundles. Public
#' users normally enter through gRm() for in-memory data or read_digram_project()
#' for an existing DIGRAM import project.
#'
#' Read a CSV file as DIGRAM data
#'
#' @param csv_path Path to the source CSV file.
#' @param items Character vector of item column names.
#' @param exo Character vector of exogeneous/person-factor column names.
#' @param idvar Identifier column name. Pass `NULL` to use the first CSV
#'   column.
#' @param output_dir Directory where `DIGRAM.csv`, `DIGRAM.imp`, and
#'   `DIGRAM.imv` should be written. Pass `NULL` when `save_digram_files` is
#'   `FALSE`.
#' @param save_digram_files Whether to write the three DIGRAM import files.
#'   If `TRUE`, `output_dir` must be supplied.
#' @param name DIGRAM project name used for saved files.
#' @param digram_folder Folder path written to `DIGRAM.imp`.
#' @param na.strings Strings treated as missing values by [utils::read.csv()].
#' @return A `gRm_data`/`gRm_project` object.
#' @keywords internal
#' @noRd
read_digram_csv <- function(csv_path,
                            items,
                            exo,
                            idvar,
                            output_dir,
                            save_digram_files,
                            name = "DIGRAM",
                            digram_folder = ".",
                            na.strings = "NA") {
  csv_path <- normalizePath(csv_path, mustWork = TRUE)
  data <- utils::read.csv(
    csv_path,
    sep = ",",
    header = TRUE,
    na.strings = na.strings,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (is.null(idvar)) {
    idvar <- names(data)[[1L]]
  }
  if (!idvar %in% names(data)) {
    stop("idvar not found in CSV data: ", idvar, call. = FALSE)
  }

  project <- build_gRm_internal_project(
    data = data,
    items = items,
    exo = exo,
    item_labels = NULL,
    exo_labels = NULL,
    item_max = NULL,
    exo_max = NULL,
    item_levels = "observed",
    exo_levels = "observed",
    paths = list(input_dir = dirname(csv_path), csv = csv_path)
  )
  project$source_data <- data
  project$import <- list(
    loader = "read_digram_csv",
    csv_path = csv_path,
    item_columns = project$items$name,
    exo_columns = project$backgrounds$name,
    idvar = idvar
  )
  project$source_trace <- c(
    data = "read_digram_csv",
    structure = "build_gRm_internal_project",
    import_bundle = "write_digram_import_files"
  )

  if (isTRUE(save_digram_files)) {
    if (is.null(output_dir)) {
      stop("output_dir must be supplied when save_digram_files is TRUE.", call. = FALSE)
    }
    written <- write_digram_import_files(
      data = data,
      project = project,
      idvar = idvar,
      output_dir = output_dir,
      name = name,
      digram_folder = digram_folder,
      line_ending = "\r\n"
    )
    project$import$output_dir <- normalizePath(output_dir, mustWork = TRUE)
    project$import$written_files <- written
  }

  project
}

#' Read a DIGRAM import-file directory
#'
#' This is the file-directory data entry point for the package. It reads a
#' directory containing `DIGRAM.csv`, `DIGRAM.imp`, and `DIGRAM.imv`, then
#' constructs the same internal DIGRAM project representation as the CSV
#' data-entry helper.
#'
#' `DIGRAM.csv`, `DIGRAM.imp`, and `DIGRAM.imv` do not encode the item versus
#' exogeneous split explicitly. The current reader therefore uses the explicit
#' `items` and `exo` declarations to assign roles while using `DIGRAM.imv` for
#' DIGRAM aliases, category maxima, and names.
#'
#' @param input_dir Directory containing `DIGRAM.csv`, `DIGRAM.imp`, and
#'   `DIGRAM.imv`.
#' @param items Character vector of item variable names.
#' @param exo Character vector of exogeneous/person-factor variable names.
#' @param idvar Identifier column name in `DIGRAM.csv`. Pass `NULL` to use the
#'   first CSV column.
#' @param name DIGRAM project name prefix.
#' @param na.strings Strings treated as missing values by [utils::read.csv()].
#' @return A `gRm_data`/`gRm_project` object.
#' @keywords internal
#' @noRd
read_digram_files <- function(input_dir,
                              items,
                              exo,
                              idvar,
                              name = "DIGRAM",
                              na.strings = "NA") {
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  imp_path <- find_digram_import_file(input_dir, paste0(name, ".imp"))
  if (!file.exists(imp_path)) {
    stop("Missing DIGRAM import file(s): ", imp_path, call. = FALSE)
  }
  imp <- read_digram_imp(imp_path)
  csv_path <- find_digram_import_file(imp$path, paste0(imp$project_name, ".csv"))
  imv_path <- find_digram_import_file(imp$path, paste0(imp$project_name, ".imv"))
  missing <- c(csv_path, imv_path)[!file.exists(c(csv_path, imv_path))]
  if (length(missing) > 0L) {
    stop("Missing DIGRAM import file(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  metadata <- read_digram_imv(imv_path)
  imv_level_markers <- attr(metadata, "level_markers", exact = TRUE) %||% character()
  declared <- c(items, exo)
  absent <- setdiff(declared, metadata$name)
  if (length(absent) > 0L) {
    stop("Declared variables not found in DIGRAM.imv: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  duplicate_names <- metadata$name[duplicated(metadata$name)]
  if (length(duplicate_names) > 0L) {
    stop("Duplicate variable names in DIGRAM.imv: ", paste(unique(duplicate_names), collapse = ", "), call. = FALSE)
  }
  metadata <- metadata[match(declared, metadata$name), , drop = FALSE]

  data <- utils::read.csv(
    csv_path,
    header = TRUE,
    na.strings = na.strings,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (is.null(idvar)) {
    idvar <- names(data)[[1L]]
  }
  if (!idvar %in% names(data)) {
    stop("idvar not found in DIGRAM.csv: ", idvar, call. = FALSE)
  }

  project <- build_gRm_internal_project(
    data = data,
    items = items,
    exo = exo,
    item_labels = metadata$label_code[seq_along(items)],
    exo_labels = metadata$label_code[length(items) + seq_along(exo)],
    item_max = metadata$raw_max[seq_along(items)],
    exo_max = metadata$raw_max[length(items) + seq_along(exo)],
    item_levels = NULL,
    exo_levels = NULL,
    paths = list(input_dir = imp$path, csv = csv_path, imp = imp_path, imv = imv_path)
  )
  project$source_data <- data
  project$import <- list(
    loader = "read_digram_files",
    digram_dir = input_dir,
    source_dir = imp$path,
    project_name = imp$project_name,
    csv_path = csv_path,
    imp_path = imp_path,
    imv_path = imv_path,
    imv_level_markers = imv_level_markers,
    item_columns = project$items$name,
    exo_columns = project$backgrounds$name,
    idvar = idvar
  )
  project$source_trace <- c(
    data = "read_digram_files",
    structure = "build_gRm_internal_project",
    import_bundle = "DIGRAM.csv/DIGRAM.imp/DIGRAM.imv"
  )
  project
}

#' Internal find digram import file helper
#'
#' Supports the digram project io implementation while preserving its internal contract.
#' @param directory Internal `directory` value used by this helper.
#' @param filename Internal `filename` value used by this helper.
#' @return The internal `find_digram_import_file()` computation result.
#' @keywords internal
#' @noRd
find_digram_import_file <- function(directory, filename) {
  candidate <- file.path(directory, filename)
  if (file.exists(candidate)) {
    return(normalizePath(candidate, mustWork = TRUE))
  }
  files <- list.files(directory, all.files = TRUE, no.. = TRUE, full.names = FALSE)
  matched <- which(tolower(files) == tolower(filename))
  if (length(matched) > 0L) {
    return(normalizePath(file.path(directory, files[[matched[[1L]]]]), mustWork = TRUE))
  }
  candidate
}

#' Internal read digram imp helper
#'
#' Supports the digram project io implementation while preserving its internal contract.
#' @param imp_path Internal `imp_path` value used by this helper.
#' @return The internal `read_digram_imp()` computation result.
#' @keywords internal
#' @noRd
read_digram_imp <- function(imp_path) {
  lines <- readLines(imp_path, warn = FALSE)
  if (length(lines) < 2L) {
    stop("Malformed DIGRAM.imp file: expected path and project name.", call. = FALSE)
  }
  project_name <- trimws(lines[[2L]])
  if (identical(project_name, "") || identical(project_name, "-")) {
    stop("Malformed DIGRAM.imp file: missing project name.", call. = FALSE)
  }
  list(
    path = resolve_digram_imp_path(lines[[1L]], dirname(imp_path)),
    project_name = project_name
  )
}

#' Internal resolve digram imp path helper
#'
#' Supports the digram project io implementation while preserving its internal contract.
#' @param path Filesystem path used by the helper.
#' @param base_dir Internal `base_dir` value used by this helper.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
resolve_digram_imp_path <- function(path, base_dir) {
  path <- trimws(path)
  if (identical(path, "") || identical(path, "-")) {
    path <- "."
  }
  path <- gsub("\\", "/", path, fixed = TRUE)
  is_windows_absolute <- grepl("^[A-Za-z]:/", path)
  is_absolute <- grepl("^/", path) || is_windows_absolute
  candidate <- if (is_absolute) path else file.path(base_dir, path)
  if (dir.exists(candidate)) {
    return(normalizePath(candidate, mustWork = TRUE))
  }

  local_candidate <- file.path(base_dir, basename(path))
  if (!identical(local_candidate, candidate) && dir.exists(local_candidate)) {
    return(normalizePath(local_candidate, mustWork = TRUE))
  }

  normalizePath(candidate, mustWork = FALSE)
}

#' Internal write digram import files helper
#'
#' Supports the digram project io implementation while preserving its internal contract.
#' @param data Input data for the computation.
#' @param project Encoded gRm project.
#' @param idvar Internal `idvar` value used by this helper.
#' @param output_dir Internal `output_dir` value used by this helper.
#' @param name Internal name or label.
#' @param digram_folder Internal `digram_folder` value used by this helper.
#' @param line_ending Internal `line_ending` value used by this helper.
#' @return The internal `write_digram_import_files()` computation result.
#' @keywords internal
#' @noRd
write_digram_import_files <- function(data,
                                      project,
                                      idvar,
                                      output_dir,
                                      name,
                                      digram_folder,
                                      line_ending) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  out <- digram_import_data(data, project, idvar)
  csv_path <- file.path(output_dir, paste0(name, ".csv"))
  utils::write.table(
    out,
    csv_path,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "",
    eol = line_ending
  )

  imp_path <- file.path(output_dir, paste0(name, ".imp"))
  writeLines(
    c(ensure_digram_folder(digram_folder), name, "-", "-"),
    imp_path,
    sep = line_ending,
    useBytes = TRUE
  )

  imv_path <- file.path(output_dir, paste0(name, ".imv"))
  writeLines(
    digram_imv_lines(project),
    imv_path,
    sep = line_ending,
    useBytes = TRUE
  )

  c(csv = csv_path, imp = imp_path, imv = imv_path)
}

#' Internal digram import data helper
#'
#' Supports the digram project io implementation while preserving its internal contract.
#' @param data Input data for the computation.
#' @param project Encoded gRm project.
#' @param idvar Internal `idvar` value used by this helper.
#' @return The internal `digram_import_data()` computation result.
#' @keywords internal
#' @noRd
digram_import_data <- function(data, project, idvar) {
  encoded_positions <- c(project$items$position, project$backgrounds$position)
  encoded_names <- c(project$items$name, project$backgrounds$name)
  encoded <- as.data.frame(
    project$raw_data[, encoded_positions, drop = FALSE],
    stringsAsFactors = FALSE
  )
  names(encoded) <- encoded_names

  # DIGRAM import CSV files represent internal missing values as blank cells.
  encoded[encoded == -999L] <- NA_integer_

  cbind(
    data.frame(id = data[[idvar]], stringsAsFactors = FALSE),
    encoded,
    stringsAsFactors = FALSE
  )
}

#' Internal ensure digram folder helper
#'
#' Supports the digram project io implementation while preserving its internal contract.
#' @param path Filesystem path used by the helper.
#' @return The internal `ensure_digram_folder()` computation result.
#' @keywords internal
#' @noRd
ensure_digram_folder <- function(path) {
  path <- gsub("/", "\\\\", path, fixed = TRUE)
  if (grepl("\\\\$", path)) path else paste0(path, "\\")
}

#' Internal digram imv lines helper
#'
#' Supports the digram project io implementation while preserving its internal contract.
#' @param project Encoded gRm project.
#' @return The internal `digram_imv_lines()` computation result.
#' @keywords internal
#' @noRd
digram_imv_lines <- function(project) {
  variables <- rbind(project$items, project$backgrounds)
  labels <- digram_category_labels()
  vapply(seq_len(nrow(variables)), function(index) {
    variable <- variables[index, , drop = FALSE]
    values <- seq_len(variable$raw_max[[1L]])
    parts <- unlist(lapply(values, function(value) {
      label <- labels[[as.character(value)]]
      if (is.null(label)) {
        stop("No DIGRAM label defined for category ", value, ".", call. = FALSE)
      }
      c(as.character(value), label)
    }), use.names = FALSE)
    paste(c(variable$label_code, variable$name, parts), collapse = ",")
  }, character(1L))
}

#' Internal digram category labels helper
#'
#' Supports the digram project io implementation while preserving its internal contract.
#' @return The internal `digram_category_labels()` computation result.
#' @keywords internal
#' @noRd
digram_category_labels <- function() {
  labels <- c(
    "one", "two", "three", "four", "five", "six", "seven", "eight",
    "nine", "ten", "eleven", "twelve", "thirteen", "fourteen",
    "fifteen", "sixteen", "seventeen", "eighteen"
  )
  stats::setNames(labels, as.character(seq_along(labels)))
}

#' Internal read digram imv helper
#'
#' Supports the digram project io implementation while preserving its internal contract.
#' @param imv_path Internal `imv_path` value used by this helper.
#' @return The internal `read_digram_imv()` computation result.
#' @keywords internal
#' @noRd
read_digram_imv <- function(imv_path) {
  lines <- readLines(imv_path, warn = FALSE)
  nonblank <- lines[nzchar(trimws(lines))]
  marker <- startsWith(trimws(nonblank), "<")
  rows <- lapply(nonblank[!marker], function(line) {
    fields <- strsplit(line, ",", fixed = TRUE)[[1L]]
    if (length(fields) < 4L || length(fields) %% 2L != 0L) {
      stop("Malformed DIGRAM.imv row: ", line, call. = FALSE)
    }
    value_fields <- trimws(fields[seq(3L, length(fields), by = 2L)])
    valid_integer_text <- grepl("^[+]?[0-9]+$", value_fields)
    if (!all(valid_integer_text)) {
      stop(
        "DIGRAM.imv category declarations must be contiguous one-based positive integer codes: ",
        line,
        call. = FALSE
      )
    }
    values <- suppressWarnings(as.integer(value_fields))
    expected_values <- seq_len(max(values, na.rm = TRUE))
    if (anyNA(values) ||
      any(values < 1L) ||
      anyDuplicated(values) ||
      !identical(sort(values), expected_values)) {
      stop(
        "DIGRAM.imv category declarations must be contiguous one-based positive integer codes: ",
        line,
        call. = FALSE
      )
    }
    data.frame(
      label_code = fields[[1L]],
      name = fields[[2L]],
      raw_max = max(values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- if (length(rows) == 0L) {
    data.frame(
      label_code = character(),
      name = character(),
      raw_max = integer(),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, rows)
  }
  attr(out, "level_markers") <- nonblank[marker]
  out
}
