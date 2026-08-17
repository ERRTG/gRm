#' Specify a GLLRM model
#'
#' @param project A `gRm_analysis`, `gRm_project`, `gRm_screen`, or fitted object.
#' @param ld One-sided formula of local-dependence item:item terms.
#' @param dif One-sided formula of DIF item:exogeneous terms.
#' @return A `gRm_model` object. This function specifies a model; it
#'   does not fit it.
#' @details
#' A GLLRM model defines a score-inclusive IRT graph. Every item is connected
#' to the score node, LD terms add item:item edges, and DIF terms add
#' item:exogenous edges. Use [model_graph()] to return that graph as an
#' `igraph` object, or `plot(model)` to draw it with the score node on the left.
#' Item and exogeneous names that are not syntactic R names can be used in
#' formulas by wrapping them in R backquote characters. For example, an item
#' named `"item one"` should be written between backquotes in the formula.
#' @seealso [model_graph()], [fit()], [update.gRm_model()]
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
#' model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I1:site)
#' summary(model)
#' graph <- model_graph(model)
#' \donttest{
#' plot(model)
#' }
gllrm <- function(project,
                  ld = NULL,
                  dif = NULL) {
  if (inherits(project, "gRm_screen")) {
    return(gllrm_from_screen(project, ld = ld, dif = dif, call = match.call()))
  }
  analysis <- as_gRm_analysis(project)
  ld_terms <- parse_ld_formula(ld, analysis, source = "user")
  dif_terms <- parse_dif_formula(dif, analysis, source = "user")
  new_gRm_model(
    analysis = analysis,
    ld = ld_terms,
    dif = dif_terms,
    screen = NULL,
    call = match.call()
  )
}

#' Extract canonical model terms for internal dispatch
#'
#' This internal generic gives package orchestration one stable representation
#' of LD, DIF, and score-effect terms across specifications, fits, and SCREEN J
#' results. Methods return source-ordered term tables; callers must not infer
#' model membership from provisional SCREEN rows outside this contract.
#'
#' @param object A `gRm_model`, `gRm_fit`, or `gRm_screen` object.
#' @param ... Reserved for S3 dispatch compatibility and must be empty.
#' @return A list with canonical `ld`, `dif`, and `score_effects` data frames
#'   plus a `source_trace` character vector.
#' @usage model_terms(object, ...)
#' @name gRm-model-terms
#' @aliases model_terms model_terms.gRm_model model_terms.gRm_fit
#'   model_terms.gRm_screen
#' @keywords internal
model_terms <- function(object, ...) {
  UseMethod("model_terms")
}

#' @export
model_terms.gRm_model <- function(object, ...) {
  reject_public_dots(...)
  list(
    ld = object$ld,
    dif = object$dif,
    score_effects = empty_score_effect_terms(),
    source_trace = object$source_trace %||% character()
  )
}

#' @export
model_terms.gRm_fit <- function(object, ...) {
  reject_public_dots(...)
  model_terms(object$model %||% object$spec)
}

#' @export
model_terms.gRm_screen <- function(object, ...) {
  reject_public_dots(...)
  values <- object$values
  list(
    ld = screen_ld_terms(values),
    dif = screen_dif_terms(values),
    score_effects = screen_score_effect_terms(values),
    source_trace = object$source_trace %||% character()
  )
}

#' Internal as gRm analysis helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
as_gRm_analysis <- function(x) {
  if (inherits(x, "gRm_analysis")) {
    return(x)
  }
  if (inherits(x, "gRm_model")) {
    return(x$analysis)
  }
  if (inherits(x, "gRm_fit")) {
    return((x$model %||% x$spec)$analysis)
  }
  if (inherits(x, "gRm_screen")) {
    return(x$analysis)
  }
  if (inherits(x, "gRm_project")) {
    return(new_gRm_analysis(
      project = x,
      data = x$source_data %||% data.frame(),
      id = x$import$idvar %||% NULL,
      data_name = x$import$project_name %||% x$paths$input_dir %||% "gRm_project",
      score_cuts = "auto",
      name = x$paths$input_dir %||% "gRm_project",
      call = match.call()
    ))
  }
  stop("Expected a gRm analysis, GLLRM object, or DIGRAM project.", call. = FALSE)
}

#' Internal new gRm model helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param analysis Prepared gRm analysis.
#' @param ld Internal `ld` value used by this helper.
#' @param dif Internal `dif` value used by this helper.
#' @param screen Internal `screen` value used by this helper.
#' @param call Captured R call.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
new_gRm_model <- function(analysis, ld, dif, screen = NULL, call) {
  terms <- source_order_model_terms(analysis, ld, dif)
  ld <- terms$ld
  dif <- terms$dif
  model_type <- if (nrow(ld) || nrow(dif)) "gllrm" else "rasch"
  out <- list(
    analysis = analysis,
    project = analysis$project,
    analysis_fingerprint = analysis$analysis_fingerprint,
    likelihood_sample = analysis$likelihood_sample,
    ld = ld,
    dif = dif,
    terms = list(ld = ld, dif = dif),
    model_type = model_type,
    screen = screen,
    source_trace = c(analysis$source_trace %||% character(), model = model_type),
    validation = list(status = "specified_not_fitted"),
    unmodeled = character(),
    warnings = character(),
    call = call
  )
  class(out) <- c("gRm_model", "list")
  out
}

#' Internal new gRm model from canonical terms helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param model gRm model specification.
#' @param ld Internal `ld` value used by this helper.
#' @param dif Internal `dif` value used by this helper.
#' @param call Captured R call.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
new_gRm_model_from_canonical_terms <- function(model, ld, dif, call = NULL) {
  if (!inherits(model, "gRm_model")) {
    stop("Expected a gRm_model object.", call. = FALSE)
  }
  ld <- ld %||% empty_ld_terms()
  dif <- dif %||% empty_dif_terms()
  terms <- source_order_model_terms(model$analysis, ld, dif)
  ld <- terms$ld
  dif <- terms$dif
  screen <- if (any(c(ld$source, dif$source) %in% "screen")) model$screen else NULL
  new_gRm_model(
    analysis = model$analysis,
    ld = ld,
    dif = dif,
    screen = screen,
    call = call %||% model$call
  )
}

#' Internal gllrm from screen helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param screen Internal `screen` value used by this helper.
#' @param ld Internal `ld` value used by this helper.
#' @param dif Internal `dif` value used by this helper.
#' @param call Captured R call.
#' @return The internal `gllrm_from_screen()` computation result.
#' @keywords internal
#' @noRd
gllrm_from_screen <- function(screen, ld = NULL, dif = NULL, call) {
  selected <- selected_screen_terms(screen)
  ld_terms <- if (is.null(ld)) selected$ld else parse_ld_formula(ld, screen$analysis, source = "user")
  dif_terms <- if (is.null(dif)) selected$dif else parse_dif_formula(dif, screen$analysis, source = "user")
  new_gRm_model(
    analysis = screen$analysis,
    ld = ld_terms,
    dif = dif_terms,
    screen = screen,
    call = call
  )
}

#' Internal selected screen terms helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param screen Internal `screen` value used by this helper.
#' @return The internal `selected_screen_terms()` computation result.
#' @keywords internal
#' @noRd
selected_screen_terms <- function(screen) {
  terms <- model_terms(screen)
  terms$ld$source <- rep("screen", nrow(terms$ld))
  terms$dif$source <- rep("screen", nrow(terms$dif))
  list(ld = terms$ld, dif = terms$dif)
}

#' Update a DIGRAM model specification
#'
#' Replace LD and/or DIF terms in an existing DIGRAM model specification.
#'
#' @param object A `gRm_model` object created by [gllrm()].
#' @param ld Replacement local-dependence formula. Use `NULL` to keep the
#'   current LD terms. Use `~ 0` to remove all LD terms.
#' @param dif Replacement DIF formula. Use `NULL` to keep the current DIF
#'   terms. Use `~ 0` to remove all DIF terms.
#' @param ... Reserved for S3 dispatch compatibility and must be empty.
#' @return An updated `gRm_model` object. The returned model is not fitted;
#'   pass it to [fit()] to estimate parameters.
#' @details
#' `update.gRm_model()` keeps the original analysis object and replaces only
#' the requested model-term families. Formula syntax is the same as in
#' [gllrm()]: LD terms are item:item interactions and DIF terms are
#' item:exogenous interactions.
#'
#' If the original model came from [screen()], screen-selected terms are
#' retained only while their source remains present in the updated model terms.
#' @examples
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0),
#'   site = c(0, 0, 1, 1, 0, 1, 0, 1)
#' )
#' analysis <- gRm(
#'   data,
#'   items = c("I1", "I2"),
#'   exogenous = "site",
#'   id = "ID",
#'   score_cuts = c(1L, 2L)
#' )
#' model <- gllrm(analysis)
#' updated <- update(model, ld = ~ I1:I2, dif = ~ I1:site)
#' summary(updated)
#' @seealso [gllrm()], [fit()]
#' @export
update.gRm_model <- function(object, ld = NULL, dif = NULL, ...) {
  reject_public_dots(...)
  ld_terms <- if (is.null(ld)) object$ld else parse_ld_formula(ld, object$analysis, source = "user")
  dif_terms <- if (is.null(dif)) object$dif else parse_dif_formula(dif, object$analysis, source = "user")
  screen <- if (any(c(ld_terms$source, dif_terms$source) %in% "screen")) object$screen else NULL
  new_gRm_model(object$analysis, ld_terms, dif_terms, screen = screen, call = match.call())
}

#' Internal parse ld formula helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param formula Internal `formula` value used by this helper.
#' @param analysis Prepared gRm analysis.
#' @param source Internal `source` value used by this helper.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
parse_ld_formula <- function(formula, analysis, source = "user") {
  parse_gRm_interaction_formula(
    formula = formula,
    items = analysis$items,
    exogenous = analysis$exogenous,
    type = "ld",
    source = source
  )
}

#' Internal parse dif formula helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param formula Internal `formula` value used by this helper.
#' @param analysis Prepared gRm analysis.
#' @param source Internal `source` value used by this helper.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
parse_dif_formula <- function(formula, analysis, source = "user") {
  parse_gRm_interaction_formula(
    formula = formula,
    items = analysis$items,
    exogenous = analysis$exogenous,
    type = "dif",
    source = source
  )
}

#' Internal parse gRm interaction formula helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param formula Internal `formula` value used by this helper.
#' @param items Item selection or item metadata.
#' @param exogenous Internal `exogenous` value used by this helper.
#' @param type Internal `type` value used by this helper.
#' @param source Internal `source` value used by this helper.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
parse_gRm_interaction_formula <- function(formula, items, exogenous, type, source = "user") {
  if (is.null(formula)) {
    return(if (identical(type, "ld")) empty_ld_terms() else empty_dif_terms())
  }
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop("`", type, "` must be a one-sided formula such as ~ item1:item2.", call. = FALSE)
  }
  if (identical(deparse(formula[[2L]]), "0")) {
    return(if (identical(type, "ld")) empty_ld_terms() else empty_dif_terms())
  }
  parsed_terms <- parse_gRm_formula_terms(formula[[2L]])
  if (length(parsed_terms) == 0L) {
    return(if (identical(type, "ld")) empty_ld_terms() else empty_dif_terms())
  }

  rows <- lapply(parsed_terms, function(term) {
    vars <- term$vars
    label <- term$label
    if (identical(type, "ld")) {
      canonical_ld_term(vars, items, label, source = source)
    } else {
      canonical_dif_term(vars, items, exogenous, label, source = source)
    }
  })
  out <- do.call(rbind, rows)
  key <- if (identical(type, "ld")) paste(out$item1, out$item2, sep = ":") else paste(out$item, out$exogenous, sep = ":")
  if (anyDuplicated(key)) {
    stop("Duplicate ", toupper(type), " terms are not allowed.", call. = FALSE)
  }
  rownames(out) <- NULL
  out
}

#' Internal parse gRm formula terms helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param expr Internal `expr` value used by this helper.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
parse_gRm_formula_terms <- function(expr) {
  if (is.numeric(expr) && length(expr) == 1L && expr %in% c(0, 1)) {
    return(list())
  }
  if (!is.call(expr)) {
    stop("DIGRAM model formulas allow only `:` interaction terms joined by `+`.", call. = FALSE)
  }

  operator <- as.character(expr[[1L]])
  if (operator %in% c("+", "-") && length(expr) == 2L) {
    value <- expr[[2L]]
    if (is.numeric(value) && length(value) == 1L && value %in% c(0, 1)) {
      return(list())
    }
    stop("DIGRAM model formulas allow only `:` interaction terms joined by `+`.", call. = FALSE)
  }
  if (identical(operator, "+") && length(expr) == 3L) {
    return(c(parse_gRm_formula_terms(expr[[2L]]), parse_gRm_formula_terms(expr[[3L]])))
  }
  if (identical(operator, ":") && length(expr) == 3L) {
    vars <- c(parse_gRm_formula_variable(expr[[2L]]), parse_gRm_formula_variable(expr[[3L]]))
    return(list(list(vars = vars, label = paste(vars, collapse = ":"))))
  }

  stop("DIGRAM model formulas allow only `:` interaction terms joined by `+`.", call. = FALSE)
}

#' Internal parse gRm formula variable helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param expr Internal `expr` value used by this helper.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
parse_gRm_formula_variable <- function(expr) {
  if (is.name(expr)) {
    return(as.character(expr))
  }
  stop("DIGRAM model formulas allow only `:` interaction terms joined by `+`.", call. = FALSE)
}

#' Internal canonical ld term helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param vars Internal `vars` value used by this helper.
#' @param items Item selection or item metadata.
#' @param label Internal `label` value used by this helper.
#' @param source Internal `source` value used by this helper.
#' @return The internal `canonical_ld_term()` computation result.
#' @keywords internal
#' @noRd
canonical_ld_term <- function(vars, items, label, source = "user") {
  if (!all(vars %in% items)) {
    stop("LD terms must be item:item interactions. Invalid term: ", label, call. = FALSE)
  }
  if (identical(vars[[1L]], vars[[2L]])) {
    stop("LD terms must contain two different items. Invalid term: ", label, call. = FALSE)
  }
  ordered <- items[sort(match(vars, items))]
  data.frame(
    type = "ld",
    item1 = ordered[[1L]],
    item2 = ordered[[2L]],
    source = source,
    status = "specified",
    stringsAsFactors = FALSE
  )
}

#' Internal canonical dif term helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param vars Internal `vars` value used by this helper.
#' @param items Item selection or item metadata.
#' @param exogenous Internal `exogenous` value used by this helper.
#' @param label Internal `label` value used by this helper.
#' @param source Internal `source` value used by this helper.
#' @return The internal `canonical_dif_term()` computation result.
#' @keywords internal
#' @noRd
canonical_dif_term <- function(vars, items, exogenous, label, source = "user") {
  item <- vars[vars %in% items]
  exo <- vars[vars %in% exogenous]
  if (length(item) != 1L || length(exo) != 1L) {
    stop("DIF terms must be item:exogenous interactions. Invalid term: ", label, call. = FALSE)
  }
  data.frame(
    type = "dif",
    item = item,
    exogenous = exo,
    source = source,
    status = "specified",
    stringsAsFactors = FALSE
  )
}

#' Internal empty ld terms helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @return The internal `empty_ld_terms()` computation result.
#' @keywords internal
#' @noRd
empty_ld_terms <- function() {
  data.frame(
    type = character(),
    item1 = character(),
    item2 = character(),
    source = character(),
    status = character(),
    stringsAsFactors = FALSE
  )
}

#' Internal empty dif terms helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @return The internal `empty_dif_terms()` computation result.
#' @keywords internal
#' @noRd
empty_dif_terms <- function() {
  data.frame(
    type = character(),
    item = character(),
    exogenous = character(),
    source = character(),
    status = character(),
    stringsAsFactors = FALSE
  )
}

#' Internal empty score effect terms helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @return The internal `empty_score_effect_terms()` computation result.
#' @keywords internal
#' @noRd
empty_score_effect_terms <- function() {
  data.frame(
    type = character(),
    exogenous = character(),
    source = character(),
    status = character(),
    stringsAsFactors = FALSE
  )
}

#' Internal screen ld terms helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @return The internal `screen_ld_terms()` computation result.
#' @keywords internal
#' @noRd
screen_ld_terms <- function(values) {
  rows <- values$model$local_dependence$rows %||% data.frame()
  if (!nrow(rows)) {
    return(empty_ld_terms())
  }
  items <- values$items
  data.frame(
    type = "ld",
    item1 = items$name[rows$row],
    item2 = items$name[rows$col],
    source = "screen",
    # `rows` is the post-filter model table. The provisional SCREEN J stage
    # remains in `values$model$local_dependence$stepwise_rows`; it is evidence
    # metadata, not a model-term status.
    status = "selected",
    stringsAsFactors = FALSE
  )
}

#' Internal screen dif terms helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @return The internal `screen_dif_terms()` computation result.
#' @keywords internal
#' @noRd
screen_dif_terms <- function(values) {
  mat <- values$model$item_bias
  if (is.null(mat) || !length(mat) || !any(mat, na.rm = TRUE)) {
    return(empty_dif_terms())
  }
  idx <- which(mat, arr.ind = TRUE)
  data.frame(
    type = "dif",
    item = values$items$name[idx[, 1L]],
    exogenous = values$backgrounds$name[idx[, 2L]],
    source = "screen",
    status = "selected",
    stringsAsFactors = FALSE
  )
}

#' Internal screen score effect terms helper
#'
#' Supports the api model spec implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @return The internal `screen_score_effect_terms()` computation result.
#' @keywords internal
#' @noRd
screen_score_effect_terms <- function(values) {
  rows <- values$model$score_effects$rows %||% data.frame()
  if (!nrow(rows)) {
    return(empty_score_effect_terms())
  }
  selected <- rows$selected %||% rep(FALSE, nrow(rows))
  data.frame(
    type = "score_effect",
    exogenous = rows$name,
    source = "screen",
    status = ifelse(selected, "selected", "not_selected"),
    stringsAsFactors = FALSE
  )
}
