# Internal project construction and category encoding
#
# These helpers convert in-memory R data and declared item/exogenous roles into
# the one-based internal project representation used by the source-faithful
# estimation code. DIGRAM file reading and writing lives in digram_project_io.R.

#' Internal build gRm internal project helper
#'
#' Supports the project input implementation while preserving its internal contract.
#' @param data Input data for the computation.
#' @param items Item selection or item metadata.
#' @param exo Internal `exo` value used by this helper.
#' @param item_labels Internal `item_labels` value used by this helper.
#' @param exo_labels Internal `exo_labels` value used by this helper.
#' @param item_max Internal `item_max` value used by this helper.
#' @param exo_max Internal `exo_max` value used by this helper.
#' @param item_levels Internal `item_levels` value used by this helper.
#' @param exo_levels Internal `exo_levels` value used by this helper.
#' @param paths Internal `paths` value used by this helper.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
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
    category_levels = list(
      items = item_levels,
      backgrounds = exo_levels
    ),
    raw_data = raw_data
  )
  class(project) <- c("gRm_data", "gRm_project", "list")
  project
}

#' Internal gRm entry aliases helper
#'
#' Supports the project input implementation while preserving its internal contract.
#' @return The internal `gRm_entry_aliases()` computation result.
#' @keywords internal
#' @noRd
gRm_entry_aliases <- function() {
  c(letters, LETTERS[1:24])
}

#' Internal check gRm entry labels helper
#'
#' Supports the project input implementation while preserving its internal contract.
#' @param labels Source-facing labels.
#' @param expected Internal `expected` value used by this helper.
#' @param name Internal name or label.
#' @return The internal `check_gRm_entry_labels()` computation result.
#' @keywords internal
#' @noRd
check_gRm_entry_labels <- function(labels, expected, name) {
  if (length(labels) != expected) {
    stop(name, " must have length ", expected, ".", call. = FALSE)
  }
  if (anyDuplicated(labels)) {
    stop(name, " must not contain duplicates.", call. = FALSE)
  }
  invisible(NULL)
}

#' Internal resolve gRm entry levels helper
#'
#' Supports the project input implementation while preserving its internal contract.
#' @param data Input data for the computation.
#' @param explicit Internal `explicit` value used by this helper.
#' @param maxima Internal `maxima` value used by this helper.
#' @param name Internal name or label.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
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

#' Internal recycle gRm entry integer helper
#'
#' Supports the project input implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @param n Internal `n` value used by this helper.
#' @param name Internal name or label.
#' @param variable_names Internal `variable_names` value used by this helper.
#' @return The internal `recycle_gRm_entry_integer()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal recycle gRm entry levels helper
#'
#' Supports the project input implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @param n Internal `n` value used by this helper.
#' @param name Internal name or label.
#' @param variable_names Internal `variable_names` value used by this helper.
#' @return The internal `recycle_gRm_entry_levels()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal observed gRm entry levels helper
#'
#' Supports the project input implementation while preserving its internal contract.
#' @param column Internal `column` value used by this helper.
#' @param name Internal name or label.
#' @param variable_name Internal `variable_name` value used by this helper.
#' @return The internal `observed_gRm_entry_levels()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal validate gRm entry level vector helper
#'
#' Supports the project input implementation while preserving its internal contract.
#' @param levels Internal `levels` value used by this helper.
#' @param name Internal name or label.
#' @param variable_name Internal `variable_name` value used by this helper.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
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

#' Internal gRm entry level key helper
#'
#' Supports the project input implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @return The internal `gRm_entry_level_key()` computation result.
#' @keywords internal
#' @noRd
gRm_entry_level_key <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  as.character(x)
}

#' Internal validate gRm entry levels observed helper
#'
#' Supports the project input implementation while preserving its internal contract.
#' @param data Input data for the computation.
#' @param levels Internal `levels` value used by this helper.
#' @param name Internal name or label.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
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

#' Internal encode gRm entry data helper
#'
#' Supports the project input implementation while preserving its internal contract.
#' @param data Input data for the computation.
#' @param levels Internal `levels` value used by this helper.
#' @param name Internal name or label.
#' @return The internal `encode_gRm_entry_data()` computation result.
#' @keywords internal
#' @noRd
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
