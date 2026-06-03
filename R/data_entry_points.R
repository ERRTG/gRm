#' Read a CSV file as DIGRAM data
#'
#' This is the CSV data entry point for the package. It reads a raw CSV file,
#' uses explicit item and exogeneous column declarations, and returns the
#' internal DIGRAM project representation used by the R implementation.
#' Optionally, it also writes the three DIGRAM import-bundle files
#' `DIGRAM.csv`, `DIGRAM.imp`, and `DIGRAM.imv`.
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
read_digram_files <- function(input_dir,
                              items,
                              exo,
                              idvar,
                              name = "DIGRAM",
                              na.strings = "NA") {
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  csv_path <- file.path(input_dir, paste0(name, ".csv"))
  imp_path <- file.path(input_dir, paste0(name, ".imp"))
  imv_path <- file.path(input_dir, paste0(name, ".imv"))
  missing <- c(csv_path, imp_path, imv_path)[!file.exists(c(csv_path, imp_path, imv_path))]
  if (length(missing) > 0L) {
    stop("Missing DIGRAM import file(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  metadata <- read_digram_imv(imv_path)
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
    paths = list(input_dir = input_dir, csv = csv_path, imp = imp_path, imv = imv_path)
  )
  project$source_data <- data
  project$import <- list(
    loader = "read_digram_files",
    digram_dir = input_dir,
    csv_path = csv_path,
    imp_path = imp_path,
    imv_path = imv_path,
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

#' @keywords internal
build_gRm_internal_project <- function(data,
                                          items,
                                          exo,
                                          item_labels,
                                          exo_labels,
                                          item_max,
                                          exo_max,
                                          item_levels = "observed",
                                          exo_levels = "observed",
                                          paths) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  if (length(items) == 0L) {
    stop("At least one item column is required.", call. = FALSE)
  }
  if (anyDuplicated(c(items, exo))) {
    stop("Item and exogeneous column declarations must be unique.", call. = FALSE)
  }
  missing <- setdiff(c(items, exo), names(data))
  if (length(missing) > 0L) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(items) + length(exo) > length(gRm_entry_aliases())) {
    stop("DIGRAM aliases are only defined for 50 variables.", call. = FALSE)
  }

  aliases <- gRm_entry_aliases()
  if (is.null(item_labels)) {
    item_labels <- aliases[seq_along(items)]
  }
  if (is.null(exo_labels)) {
    exo_labels <- aliases[length(items) + seq_along(exo)]
  }
  check_gRm_entry_labels(item_labels, length(items), "item_labels")
  check_gRm_entry_labels(exo_labels, length(exo), "exo_labels")

  item_levels <- resolve_gRm_entry_levels(data[items], item_levels, item_max, "item_levels")
  exo_levels <- resolve_gRm_entry_levels(data[exo], exo_levels, exo_max, "exogenous_levels")
  item_max <- lengths(item_levels)
  exo_max <- lengths(exo_levels)

  selected <- c(items, exo)
  entry_levels <- c(item_levels, exo_levels)
  raw_data <- encode_gRm_entry_data(data[selected], entry_levels, "levels")
  colnames(raw_data) <- NULL

  specs <- data.frame(
    label_code = c(item_labels, exo_labels),
    position = seq_along(selected),
    raw_max = unname(c(item_max, exo_max)),
    vtype = rep.int(3L, length(selected)),
    name = selected,
    is_item = c(rep(TRUE, length(items)), rep(FALSE, length(exo))),
    stringsAsFactors = FALSE
  )

  project <- list(
    paths = paths,
    variables = specs,
    items = specs[specs$is_item, , drop = FALSE],
    backgrounds = specs[!specs$is_item, , drop = FALSE],
    raw_data = raw_data
  )
  class(project) <- c("gRm_data", "gRm_project", "list")
  project
}

#' @keywords internal
gRm_entry_aliases <- function() {
  c(letters, LETTERS[1:24])
}

#' @keywords internal
check_gRm_entry_labels <- function(labels, expected, name) {
  if (length(labels) != expected) {
    stop(name, " must have length ", expected, ".", call. = FALSE)
  }
  if (anyDuplicated(labels)) {
    stop(name, " must not contain duplicates.", call. = FALSE)
  }
  invisible(NULL)
}

#' @keywords internal
resolve_gRm_entry_levels <- function(data, explicit, maxima, name) {
  n <- ncol(data)
  if (n == 0L) {
    return(stats::setNames(list(), names(data)))
  }
  variable_names <- names(data)
  if (!is.null(maxima)) {
    maxima <- recycle_gRm_entry_integer(maxima, n, paste0(name, " maxima"), variable_names = variable_names)
    if (any(is.na(maxima) | maxima < 1L)) {
      stop(name, " maxima must be positive integers.", call. = FALSE)
    }
    levels <- lapply(maxima, seq_len)
    names(levels) <- variable_names
    validate_gRm_entry_levels_observed(data, levels, name)
    return(levels)
  }

  recycled <- recycle_gRm_entry_levels(explicit %||% "observed", n, name, variable_names)
  levels <- vector("list", n)
  names(levels) <- variable_names
  for (index in seq_len(n)) {
    spec <- recycled[[index]]
    levels[[index]] <- if (is.character(spec) && length(spec) == 1L && identical(spec, "observed")) {
      observed_gRm_entry_levels(data[[index]], name, variable_names[[index]])
    } else {
      validate_gRm_entry_level_vector(spec, name, variable_names[[index]])
    }
  }
  validate_gRm_entry_levels_observed(data, levels, name)
  levels
}

#' @keywords internal
recycle_gRm_entry_integer <- function(x, n, name, variable_names = NULL) {
  if (n == 0L) {
    return(integer())
  }
  if (!is.null(names(x)) && any(nzchar(names(x)))) {
    x_names <- names(x)
    if (any(!nzchar(x_names))) {
      stop(name, " must be fully named or fully unnamed.", call. = FALSE)
    }
    if (is.null(variable_names) || !setequal(x_names, variable_names)) {
      stop(name, " names must match declared variables.", call. = FALSE)
    }
    x <- x[variable_names]
  } else if (length(x) == 1L) {
    x <- rep(x, n)
  } else if (length(x) != n) {
    stop(name, " must have length 1 or ", n, ".", call. = FALSE)
  }
  as.integer(x)
}

#' @keywords internal
recycle_gRm_entry_levels <- function(x, n, name, variable_names) {
  if (n == 0L) {
    return(list())
  }
  if (is.list(x)) {
    if (!is.null(names(x)) && any(nzchar(names(x)))) {
      x_names <- names(x)
      if (any(!nzchar(x_names))) {
        stop(name, " must be fully named or fully unnamed.", call. = FALSE)
      }
      if (!setequal(x_names, variable_names)) {
        stop(name, " names must match declared variables.", call. = FALSE)
      }
      return(x[variable_names])
    }
    if (length(x) == 1L) {
      return(rep(x, n))
    }
    if (length(x) != n) {
      stop(name, " must have length 1 or ", n, ".", call. = FALSE)
    }
    return(x)
  }
  rep(list(x), n)
}

#' @keywords internal
observed_gRm_entry_levels <- function(column, name, variable_name) {
  values <- column[!is.na(column)]
  if (length(values) == 0L) {
    stop("Cannot infer ", name, " for ", variable_name, " from an all-missing column.", call. = FALSE)
  }
  if (is.factor(column)) {
    levels <- levels(column)[levels(column) %in% as.character(values)]
  } else {
    levels <- unique(values)
    if (is.numeric(values) || is.integer(values) || is.logical(values)) {
      levels <- sort(levels)
    } else {
      levels <- sort(as.character(levels))
    }
  }
  validate_gRm_entry_level_vector(levels, name, variable_name)
}

#' @keywords internal
validate_gRm_entry_level_vector <- function(levels, name, variable_name) {
  if (length(levels) == 0L) {
    stop(name, " for ", variable_name, " must contain at least one level.", call. = FALSE)
  }
  if (anyNA(levels)) {
    stop(name, " for ", variable_name, " must not contain missing levels.", call. = FALSE)
  }
  keys <- gRm_entry_level_key(levels)
  if (anyDuplicated(keys)) {
    stop(name, " for ", variable_name, " must not contain duplicate levels.", call. = FALSE)
  }
  levels
}

#' @keywords internal
gRm_entry_level_key <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  as.character(x)
}

#' @keywords internal
validate_gRm_entry_levels_observed <- function(data, levels, name) {
  for (index in seq_along(data)) {
    values <- data[[index]][!is.na(data[[index]])]
    if (length(values) == 0L) {
      next
    }
    value_keys <- gRm_entry_level_key(values)
    level_keys <- gRm_entry_level_key(levels[[index]])
    missing <- setdiff(unique(value_keys), level_keys)
    if (length(missing) > 0L) {
      stop(
        name,
        " for ",
        names(data)[[index]],
        " do not cover observed values: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

#' @keywords internal
encode_gRm_entry_data <- function(data, levels, name) {
  raw_data <- matrix(-999L, nrow = nrow(data), ncol = ncol(data))
  for (index in seq_along(data)) {
    column <- data[[index]]
    present <- !is.na(column)
    if (!any(present)) {
      next
    }
    encoded <- match(gRm_entry_level_key(column[present]), gRm_entry_level_key(levels[[index]]))
    if (anyNA(encoded)) {
      stop(name, " for ", names(data)[[index]], " do not cover observed values.", call. = FALSE)
    }
    raw_data[present, index] <- as.integer(encoded)
  }
  raw_data
}

#' @keywords internal
validate_gRm_variable_names <- function(variables, available, name, required) {
  if (is.null(variables)) {
    variables <- character()
  }
  if (!is.character(variables)) {
    stop(name, " must be a character vector.", call. = FALSE)
  }
  if (required && length(variables) == 0L) {
    stop(name, " must contain at least one variable.", call. = FALSE)
  }
  if (anyNA(variables) || any(variables == "")) {
    stop(name, " must contain non-missing column names.", call. = FALSE)
  }
  missing <- setdiff(variables, available)
  if (length(missing) > 0L) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  variables
}

#' @keywords internal
write_digram_import_files <- function(data,
                                      project,
                                      idvar,
                                      output_dir,
                                      name,
                                      digram_folder,
                                      line_ending) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  out <- data[c(idvar, project$items$name, project$backgrounds$name)]
  names(out)[[1L]] <- "id"
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

#' @keywords internal
ensure_digram_folder <- function(path) {
  path <- gsub("/", "\\\\", path, fixed = TRUE)
  if (grepl("\\\\$", path)) path else paste0(path, "\\")
}

#' @keywords internal
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

#' @keywords internal
digram_category_labels <- function() {
  labels <- c(
    "zero", "one", "two", "three", "four", "five", "six", "seven",
    "eight", "nine", "ten", "eleven", "twelve", "thirteen",
    "fourteen", "fifteen", "sixteen", "seventeen", "eighteen"
  )
  stats::setNames(labels, as.character(seq_along(labels) - 1L))
}

#' @keywords internal
read_digram_imv <- function(imv_path) {
  lines <- readLines(imv_path, warn = FALSE)
  rows <- lapply(lines[nzchar(trimws(lines))], function(line) {
    fields <- strsplit(line, ",", fixed = TRUE)[[1L]]
    if (length(fields) < 4L || length(fields) %% 2L != 0L) {
      stop("Malformed DIGRAM.imv row: ", line, call. = FALSE)
    }
    values <- suppressWarnings(as.integer(fields[seq(3L, length(fields), by = 2L)]))
    data.frame(
      label_code = fields[[1L]],
      name = fields[[2L]],
      raw_max = max(values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
