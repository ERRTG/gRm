#' Create a gRm item-analysis project
#'
#' @param data A data frame containing item and optional exogeneous variables.
#' @param items Character vector of item column names.
#' @param exogenous Character vector of exogeneous/person-factor column names.
#' @param id Optional identifier column name.
#' @param item_levels Item response levels. Use `"observed"` to infer sorted
#'   observed levels per item, an atomic vector to apply the same levels to all
#'   items, or a named list with one entry per item. Levels are mapped to
#'   DIGRAM's one-based raw categories before fitting.
#' @param exogenous_levels Exogeneous/person-factor levels, with the same
#'   conventions as `item_levels`.
#' @param score_cuts `"auto"` for source-faithful automatic score cuts, or an
#'   integer-like vector of explicit upper total-score cut values. The cuts are
#'   stored with the analysis and later define score groups as consecutive
#'   intervals: the first group runs from the source-valid lowest score through
#'   the first cut, the next group starts at the following score and runs through
#'   the next cut, and so on.
#' @param name Optional project name.
#' @return A `gRm_analysis` object.
#' @export
#' @examples
#' data <- data.frame(
#'   ID = 1:6,
#'   I1 = c(0, 1, 0, 1, 0, 1),
#'   I2 = c(1, 0, 1, 0, 1, 0),
#'   site = c(0, 0, 1, 1, 0, 1)
#' )
#' analysis <- gRm(data, items = c("I1", "I2"), exogenous = "site", id = "ID")
#' summary(analysis)
gRm <- function(data,
                    items,
                    exogenous = character(),
                    id = NULL,
                    item_levels = "observed",
                    exogenous_levels = "observed",
                    score_cuts = "auto",
                    name = NULL) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (missing(items) || length(items) == 0L) {
    stop("`items` must name at least one item column.", call. = FALSE)
  }
  if (!is.null(id) && !id %in% names(data)) {
    stop("`id` was not found in `data`: ", id, call. = FALSE)
  }

  project <- build_gRm_internal_project(
    data = data,
    items = as.character(items),
    exo = as.character(exogenous %||% character()),
    item_labels = NULL,
    exo_labels = NULL,
    item_max = NULL,
    exo_max = NULL,
    item_levels = item_levels,
    exo_levels = exogenous_levels,
    paths = list(input_dir = NA_character_, csv = NA_character_)
  )
  project$source_data <- data
  project$import <- list(
    loader = "gRm",
    item_columns = project$items$name,
    exo_columns = project$backgrounds$name,
    idvar = id
  )
  project$source_trace <- c(
    data = "gRm",
    structure = "build_gRm_internal_project"
  )

  new_gRm_analysis(
    project = project,
    data = data,
    id = id,
    score_cuts = score_cuts,
    name = name %||% "gRm",
    call = match.call()
  )
}

#' Read an existing DIGRAM import project
#'
#' @param path Directory containing `DIGRAM.csv`, `DIGRAM.imp`, and `DIGRAM.imv`.
#' @param items Character vector of item variable names.
#' @param exogenous Character vector of exogeneous/person-factor variable names.
#' @param id Identifier column name. Defaults to the first CSV column.
#' @param score_cuts `"auto"` for source-faithful automatic score cuts, or an
#'   integer-like vector of explicit upper total-score cut values. The cuts are
#'   stored with the analysis and later define score groups as consecutive
#'   intervals: the first group runs from the source-valid lowest score through
#'   the first cut, the next group starts at the following score and runs through
#'   the next cut, and so on.
#' @param name DIGRAM file prefix.
#' @return A `gRm_analysis` object.
#' @export
#' @examples
#' \dontrun{
#' analysis <- read_digram_project(
#'   "path/to/DIGRAM",
#'   items = c("I1", "I2"),
#'   exogenous = "site",
#'   id = "id"
#' )
#' }
read_digram_project <- function(path,
                                items,
                                exogenous = character(),
                                id = NULL,
                                score_cuts = "auto",
                                name = "DIGRAM") {
  path <- normalizePath(path, mustWork = TRUE)
  if (missing(items) || is.null(items) || length(items) == 0L) {
    stop("`items` must be supplied.", call. = FALSE)
  }
  roles <- resolve_gRm_project_roles(path, items, exogenous, id)
  project <- read_digram_files(
    input_dir = path,
    items = roles$items,
    exo = roles$exogenous,
    idvar = roles$id,
    name = name
  )
  new_gRm_analysis(
    project = project,
    data = project$source_data,
    id = project$import$idvar,
    score_cuts = score_cuts,
    name = basename(path),
    call = match.call()
  )
}

#' Sum-score specification
#'
#' @return A score specification object.
#' @noRd
sum_score <- function() {
  structure(list(type = "sum_score"), class = "gRm_score_spec")
}

#' Source-faithful automatic score groups
#'
#' @return A score-group specification object.
#' @noRd
score_groups_auto <- function() {
  structure(list(type = "auto"), class = "gRm_score_group_spec")
}

#' Explicit score-group cut values
#'
#' @param cut_values Integer-like score cut values.
#' @return A score-group specification object.
#' @noRd
score_groups_cut <- function(cut_values) {
  if (missing(cut_values) || length(cut_values) == 0L) {
    stop("`cut_values` must contain at least one score cut.", call. = FALSE)
  }
  if (anyNA(cut_values) || any(cut_values != as.integer(cut_values))) {
    stop("`cut_values` must be integer-like and non-missing.", call. = FALSE)
  }
  cuts <- as.integer(cut_values)
  if (is.unsorted(cuts, strictly = TRUE)) {
    stop("`cut_values` must be strictly increasing.", call. = FALSE)
  }
  structure(list(type = "cut", cuts = cuts), class = "gRm_score_group_spec")
}

new_gRm_item_analysis <- function(project, data, id, score = sum_score(), score_cuts = score_groups_auto(), name, call) {
  new_gRm_analysis(project, data, id, score_cuts, name, call)
}

new_gRm_analysis <- function(project, data, id, score_cuts, name, call) {
  out <- list(
    data = data,
    project = project,
    name = name,
    items = project$items$name,
    exogenous = project$backgrounds$name,
    id = id,
    score_groups = normalize_gRm_score_cuts(score_cuts, project),
    source_trace = c(project$source_trace %||% character(), api = "gRm_analysis"),
    validation = list(status = "not_validated", corpus = NA_character_),
    unmodeled = character(),
    warnings = character(),
    call = call
  )
  class(out) <- c("gRm_analysis", "list")
  out
}

validate_score_spec <- function(score) {
  if (!inherits(score, "gRm_score_spec") || !identical(score$type, "sum_score")) {
    stop("Only `sum_score()` is currently supported.", call. = FALSE)
  }
  score
}

validate_score_group_spec <- function(score_cuts) {
  if (!inherits(score_cuts, "gRm_score_group_spec")) {
    stop("`score_cuts` must be created by `score_groups_auto()` or `score_groups_cut()`.", call. = FALSE)
  }
  score_cuts
}

resolve_gRm_score_groups <- function(project, score_cuts) {
  normalize_gRm_score_cuts(score_cuts, project)
}

normalize_gRm_score_cuts <- function(score_cuts, project) {
  if (identical(score_cuts, "auto")) {
    return(resolve_auto_score_groups(project))
  }
  if (is.numeric(score_cuts) || is.integer(score_cuts)) {
    return(resolve_explicit_score_groups(project, score_cuts))
  }
  if (inherits(score_cuts, "gRm_score_group_spec")) {
    if (identical(score_cuts$type, "auto")) {
      return(resolve_auto_score_groups(project))
    }
    if (identical(score_cuts$type, "cut")) {
      return(resolve_explicit_score_groups(project, score_cuts$cuts))
    }
  }
  stop("`score_cuts` must be \"auto\" or an integer-like vector of score cuts.", call. = FALSE)
}

resolve_auto_score_groups <- function(project) {
  tryCatch(
    as.integer(items_select_values(project)$score_groups$to_score),
    error = function(e) integer()
  )
}

resolve_explicit_score_groups <- function(project, score_cuts) {
  if (length(score_cuts) == 0L || anyNA(score_cuts) || any(score_cuts != as.integer(score_cuts))) {
    stop("`score_cuts` cut values must be non-missing integer-like values.", call. = FALSE)
  }
  cuts <- as.integer(score_cuts)
  if (is.unsorted(cuts, strictly = TRUE)) {
    stop("`score_cuts` cut values must be strictly increasing.", call. = FALSE)
  }
  max_score <- sum(project$items$raw_max - 1L)
  if (any(cuts < 0L | cuts > max_score)) {
    stop("`score_cuts` cut values must lie within the possible score range 0..", max_score, ".", call. = FALSE)
  }
  cuts
}

resolve_gRm_project_roles <- function(path, items, exogenous, id) {
  list(items = as.character(items), exogenous = as.character(exogenous %||% character()), id = id)
}

#' @export
print.gRm_analysis <- function(x, ...) {
  cat("<gRm_analysis>\n")
  cat("  rows: ", nrow(x$project$raw_data), "\n", sep = "")
  cat("  items: ", paste(x$items, collapse = ", "), "\n", sep = "")
  cat("  exogenous: ", paste(x$exogenous, collapse = ", "), "\n", sep = "")
  cat("  score range: 0..", sum(x$project$items$raw_max - 1L), "\n", sep = "")
  invisible(x)
}

#' @export
print.gRm_item_analysis <- function(x, ...) {
  print.gRm_analysis(x, ...)
}
