#' Create a gRm item-analysis project
#'
#' @param data A data frame containing item and optional exogeneous variables.
#'   The current package version supports ordinal variables only. Nominal and
#'   mixed DIGRAM variable types are not implemented.
#' @param items Character vector of item column names.
#' @param exogenous Character vector of exogeneous/person-factor column names.
#'   Exogeneous variables are interpreted as ordinal source variables in the
#'   current package version. Multi-category nominal exogeneous variables are
#'   outside the currently source-faithful scope.
#' @param id Optional identifier column name.
#' @param item_levels Item response levels. Use `"observed"` to infer sorted
#'   observed levels per item, an atomic vector to apply the same levels to all
#'   items, or a named list with one entry per item. Levels are mapped to
#'   DIGRAM's one-based raw categories before fitting.
#' @param exogenous_levels Exogeneous/person-factor levels, with the same
#'   conventions as `item_levels`. These levels define category order for
#'   ordinal analysis; they do not declare nominal variable types.
#' @param score_cuts `"auto"` for source-faithful automatic score cuts, or an
#'   integer-like vector of explicit upper total-score cut values. The cuts are
#'   stored with the analysis and later define score groups as consecutive
#'   intervals: the first group runs from the source-valid lowest score through
#'   the first cut, the next group starts at the following score and runs through
#'   the next cut, and so on. Automatic cuts must define at least two usable
#'   source score groups; otherwise construction fails with a non-estimable
#'   score-group error.
#' @param name Optional project name.
#' @details
#' `item_levels` and `exogenous_levels` are construction-time encoding
#' controls. They define how user-facing data values are mapped to DIGRAM's
#' internal one-based category codes. For items, this order also defines the
#' score scale: category 1 is scored as 0, category 2 as 1, and so on. After the
#' `gRm_analysis` object is built, later analyses use the encoded
#' `project$raw_data` and variable category counts, not the original level
#' objects.
#' @section Workflow:
#' A `gRm_analysis` stores the encoded DIGRAM project and score-group setup.
#' From there, either specify a model manually with [gllrm()] or run [screen()]
#' to discover candidate LD and DIF terms. [screen()] returns a screening
#' result; pass it to [gllrm()] to create the selected model. Model diagnostics
#' require a fitted model from [fit()]. [score_effects()] is analysis-level and
#' can be called directly on the `gRm_analysis` object.
#'
#' ```
#' manual_model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I1:site)
#' screened <- screen(analysis)
#' screened_model <- gllrm(screened)
#'
#' fitted <- fit(screened_model)
#'
#' summary(fitted)
#' item_fit(fitted)
#' local_dependence(fitted)
#' dif(fitted)
#' global_homogeneity(fitted)
#' ```
#' @return A `gRm_analysis` object.
#' @export
#' @examples
#' data <- data.frame(
#'   ID = 1:6,
#'   I1 = c(0, 1, 0, 1, 0, 1),
#'   I2 = c(1, 0, 1, 0, 1, 0),
#'   site = c(0, 0, 1, 1, 0, 1)
#' )
#' analysis <- gRm(
#'   data,
#'   items = c("I1", "I2"),
#'   exogenous = "site",
#'   id = "ID",
#'   score_cuts = c(1L, 2L)
#' )
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
    data_name = gRm_data_argument_label(substitute(data)),
    score_cuts = score_cuts,
    name = name %||% "gRm",
    call = match.call()
  )
}

#' Read an existing DIGRAM import project
#'
#' @param path Directory containing the initial `DIGRAM.imp` import file. The
#'   reader follows the source convention that this file names the project
#'   directory and project prefix used to locate the corresponding `.csv` and
#'   `.imv` files.
#'   The current package reader supports the ordinal-analysis subset of DIGRAM
#'   projects. Nominal and mixed DIGRAM variable-type behavior is not
#'   implemented. Category codes in `DIGRAM.imv` must be contiguous one-based
#'   codes matching the values in `DIGRAM.csv`; zero-based, non-contiguous, or
#'   separately recoded historical category mappings are outside the currently
#'   source-faithful import subset. Recursive-level marker rows in `.imv` files
#'   are accepted as structural metadata rather than variables.
#' @param items Character vector of item variable names.
#' @param exogenous Character vector of exogeneous/person-factor variable names.
#'   Exogeneous variables are interpreted as ordinal source variables in the
#'   current package version. Multi-category nominal exogeneous variables are
#'   outside the currently source-faithful scope.
#' @param id Optional identifier column name. Defaults to the first CSV column
#'   when omitted.
#' @param score_cuts `"auto"` for source-faithful automatic score cuts, or an
#'   integer-like vector of explicit upper total-score cut values. The cuts are
#'   stored with the analysis and later define score groups as consecutive
#'   intervals: the first group runs from the source-valid lowest score through
#'   the first cut, the next group starts at the following score and runs through
#'   the next cut, and so on. Automatic cuts must define at least two usable
#'   source score groups; otherwise construction fails with a non-estimable
#'   score-group error.
#' @param name Prefix of the initial `.imp` file to open. Defaults to
#'   `"DIGRAM"`. The analysis name and the `.csv`/`.imv` file prefix are then
#'   read from that `.imp` file.
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
  items <- as.character(items)
  exogenous <- as.character(exogenous %||% character())
  project <- read_digram_files(
    input_dir = path,
    items = items,
    exo = exogenous,
    idvar = id,
    name = name
  )
  new_gRm_analysis(
    project = project,
    data = project$source_data,
    id = project$import$idvar,
    data_name = project$import$project_name %||% basename(path),
    score_cuts = score_cuts,
    name = project$import$project_name %||% basename(path),
    call = match.call()
  )
}

#' Internal normalize public integer like helper
#'
#' Supports the api constructors implementation while preserving its internal contract.
#' @param value Value to validate or transform.
#' @param message Internal `message` value used by this helper.
#' @param min_length Internal `min_length` value used by this helper.
#' @param scalar Internal `scalar` value used by this helper.
#' @param lower Internal `lower` value used by this helper.
#' @param upper Internal `upper` value used by this helper.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
normalize_public_integer_like <- function(value,
                                          message,
                                          min_length = 1L,
                                          scalar = FALSE,
                                          lower = NULL,
                                          upper = .Machine$integer.max) {
  fail <- function() stop(message, call. = FALSE)

  if (!is.numeric(value) && !is.integer(value)) {
    fail()
  }
  if (isTRUE(scalar)) {
    if (length(value) != 1L) {
      fail()
    }
  } else if (length(value) < min_length) {
    fail()
  }
  if (
    anyNA(value) ||
      any(!is.finite(value)) ||
      any(value != floor(value))
  ) {
    fail()
  }
  if (!is.null(lower) && any(value < lower)) {
    fail()
  }
  if (!is.null(upper) && any(value > upper)) {
    fail()
  }
  if (any(value < -.Machine$integer.max)) {
    fail()
  }

  as.integer(value)
}

#' Internal gRm data argument label helper
#'
#' Supports the api constructors implementation while preserving its internal contract.
#' @param expr Internal `expr` value used by this helper.
#' @return The internal `gRm_data_argument_label()` computation result.
#' @keywords internal
#' @noRd
gRm_data_argument_label <- function(expr) {
  label <- paste(deparse(expr, width.cutoff = 60L), collapse = "")
  if (!nzchar(label)) {
    return("<unnamed>")
  }
  label
}

#' Internal new gRm analysis helper
#'
#' Supports the api constructors implementation while preserving its internal contract.
#' @param project Encoded gRm project.
#' @param data Input data for the computation.
#' @param id Internal `id` value used by this helper.
#' @param data_name Internal `data_name` value used by this helper.
#' @param score_cuts Resolved total-score cut values.
#' @param name Internal name or label.
#' @param call Captured R call.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
new_gRm_analysis <- function(project, data, id, data_name, score_cuts, name, call) {
  score_groups <- normalize_gRm_score_cuts(score_cuts, project)
  bundle <- build_item_parameters_bundle(project)
  identity <- gRm_analysis_identity_fields(
    project,
    data,
    id,
    score_groups,
    bundle = bundle
  )
  out <- list(
    data = data,
    project = project,
    name = name,
    data_name = data_name,
    items = project$items$name,
    exogenous = project$backgrounds$name,
    id = id,
    score_groups = score_groups,
    analysis_identity = identity$identity,
    analysis_fingerprint = identity$fingerprint,
    likelihood_sample = identity$likelihood_sample,
    source_trace = c(project$source_trace %||% character(), api = "gRm_analysis"),
    validation = list(status = "not_validated", corpus = NA_character_),
    unmodeled = character(),
    warnings = character(),
    call = call
  )
  class(out) <- c("gRm_analysis", "list")
  out
}

#' Internal normalize gRm score cuts helper
#'
#' Supports the api constructors implementation while preserving its internal contract.
#' @param score_cuts Resolved total-score cut values.
#' @param project Encoded gRm project.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
normalize_gRm_score_cuts <- function(score_cuts, project) {
  if (identical(score_cuts, "auto")) {
    return(resolve_auto_score_groups(project))
  }
  if (is.numeric(score_cuts) || is.integer(score_cuts)) {
    return(resolve_explicit_score_groups(project, score_cuts))
  }
  stop("`score_cuts` must be \"auto\" or an integer-like vector of score cuts.", call. = FALSE)
}

#' Internal resolve auto score groups helper
#'
#' Supports the api constructors implementation while preserving its internal contract.
#' @param project Encoded gRm project.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
resolve_auto_score_groups <- function(project) {
  values <- tryCatch(
    items_select_values(project),
    error = function(e) {
      stop(
        "Automatic `score_cuts` could not be determined: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  validate_gRm_constructor_score_cuts(
    project,
    values$score_groups$to_score,
    what = "Automatic `score_cuts`",
    require_multiple_groups = TRUE
  )
}

#' Internal resolve explicit score groups helper
#'
#' Supports the api constructors implementation while preserving its internal contract.
#' @param project Encoded gRm project.
#' @param score_cuts Resolved total-score cut values.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
resolve_explicit_score_groups <- function(project, score_cuts) {
  validate_gRm_constructor_score_cuts(
    project,
    score_cuts,
    what = "`score_cuts` cut values",
    require_multiple_groups = FALSE
  )
}

#' Internal validate gRm constructor score cuts helper
#'
#' Supports the api constructors implementation while preserving its internal contract.
#' @param project Encoded gRm project.
#' @param score_cuts Resolved total-score cut values.
#' @param what Internal `what` value used by this helper.
#' @param require_multiple_groups Internal `require_multiple_groups` value used by this helper.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
validate_gRm_constructor_score_cuts <- function(project,
                                                score_cuts,
                                                what,
                                                require_multiple_groups) {
  cuts <- normalize_public_integer_like(
    score_cuts,
    paste0(what, " must contain at least two non-missing integer-like values."),
    min_length = 2L
  )
  if (is.unsorted(cuts, strictly = TRUE)) {
    stop(what, " must be strictly increasing.", call. = FALSE)
  }
  max_score <- sum(project$items$raw_max - 1L)
  if (any(cuts < 0L | cuts > max_score)) {
    stop(what, " must lie within the possible score range 0..", max_score, ".", call. = FALSE)
  }
  if (isTRUE(require_multiple_groups)) {
    bundle <- build_item_parameters_bundle(project)
    groups <- tryCatch(
      global_homogeneity_score_groups(bundle, cuts),
      error = function(e) data.frame()
    )
    if (nrow(groups) < 2L) {
      stop(what, " must define at least two usable source score groups.", call. = FALSE)
    }
  }
  cuts
}
