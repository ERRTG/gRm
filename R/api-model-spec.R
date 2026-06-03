#' Specify a GLLRM model
#'
#' @param project A `gRm_analysis`, `gRm_project`, `gRm_screen`, or fitted object.
#' @param ld One-sided formula of local-dependence item:item terms.
#' @param dif One-sided formula of DIF item:exogeneous terms.
#' @return A `gRm_model` object. This function specifies a model; it
#'   does not fit it.
#' @export
#' @examples
#' data <- data.frame(
#'   ID = 1:6,
#'   I1 = c(0, 1, 0, 1, 0, 1),
#'   I2 = c(1, 0, 1, 0, 1, 0),
#'   group = c(0, 0, 1, 1, 0, 1)
#' )
#' analysis <- gRm(data, items = c("I1", "I2"), exogenous = "group", id = "ID")
#' model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I1:group)
#' summary(model)
gllrm <- function(project,
                  ld = NULL,
                  dif = NULL) {
  if (inherits(project, "gRm_screen")) {
    return(gllrm_from_screen(project, ld = ld, dif = dif, call = match.call()))
  }
  analysis <- as_gRm_analysis(project)
  new_gRm_model(
    analysis = analysis,
    ld = parse_ld_formula(ld, analysis, source = "user"),
    dif = parse_dif_formula(dif, analysis, source = "user"),
    screen = NULL,
    call = match.call()
  )
}

#' Extract LD, DIF, and score-effect terms
#'
#' @param object A clean gRm model, fit, or screen object.
#' @param ... Reserved for S3 dispatch compatibility; ignored.
#' @return A list of canonical term tables.
#' @noRd
model_terms <- function(object, ...) {
  UseMethod("model_terms")
}

model_terms.gRm_model <- function(object, ...) {
  list(
    ld = object$ld,
    dif = object$dif,
    score_effects = empty_score_effect_terms(),
    source_trace = object$source_trace %||% character()
  )
}

model_terms.gRm_gllrm_spec <- function(object, ...) {
  model_terms.gRm_model(object, ...)
}

model_terms.gRm_fit <- function(object, ...) {
  model_terms(object$model %||% object$spec)
}

model_terms.gRm_gllrm_fit <- function(object, ...) {
  model_terms(object$model %||% object$spec)
}

model_terms.gRm_screen <- function(object, ...) {
  values <- object$values
  list(
    ld = screen_ld_terms(values),
    dif = screen_dif_terms(values),
    score_effects = screen_score_effect_terms(values),
    source_trace = object$source_trace %||% character()
  )
}

as_gRm_analysis <- function(x) {
  if (inherits(x, "gRm_analysis") || inherits(x, "gRm_item_analysis")) {
    return(x)
  }
  if (inherits(x, "gRm_model") || inherits(x, "gRm_gllrm_spec")) {
    return(x$analysis)
  }
  if (inherits(x, "gRm_fit") || inherits(x, "gRm_gllrm_fit")) {
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
      groups = "auto",
      name = x$paths$input_dir %||% "gRm_project",
      call = match.call()
    ))
  }
  stop("Expected a gRm analysis, GLLRM object, or DIGRAM project.", call. = FALSE)
}

new_gRm_model <- function(analysis, ld, dif, screen = NULL, call) {
  model_type <- if (nrow(ld) || nrow(dif)) "gllrm" else "rasch"
  out <- list(
    analysis = analysis,
    project = analysis$project,
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
#' @param ... Reserved for S3 dispatch compatibility; ignored.
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
#'   group = c(0, 0, 1, 1, 0, 1, 0, 1)
#' )
#' analysis <- gRm(data, items = c("I1", "I2"), exogenous = "group", id = "ID")
#' model <- gllrm(analysis)
#' updated <- update(model, ld = ~ I1:I2, dif = ~ I1:group)
#' summary(updated, which = "model")
#' @seealso [gllrm()], [fit()]
#' @export
update.gRm_model <- function(object, ld = NULL, dif = NULL, ...) {
  ld_terms <- if (is.null(ld)) object$ld else parse_ld_formula(ld, object$analysis, source = "user")
  dif_terms <- if (is.null(dif)) object$dif else parse_dif_formula(dif, object$analysis, source = "user")
  screen <- if (any(c(ld_terms$source, dif_terms$source) %in% "screen")) object$screen else NULL
  new_gRm_model(object$analysis, ld_terms, dif_terms, screen = screen, call = match.call())
}

parse_ld_formula <- function(formula, analysis, source = "user") {
  parse_gRm_interaction_formula(
    formula = formula,
    items = analysis$items,
    exogenous = analysis$exogenous,
    type = "ld",
    source = source
  )
}

parse_dif_formula <- function(formula, analysis, source = "user") {
  parse_gRm_interaction_formula(
    formula = formula,
    items = analysis$items,
    exogenous = analysis$exogenous,
    type = "dif",
    source = source
  )
}

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
  term_labels <- attr(stats::terms(formula, specials = NULL), "term.labels")
  if (length(term_labels) == 0L) {
    return(if (identical(type, "ld")) empty_ld_terms() else empty_dif_terms())
  }
  if (any(grepl("[*/^()]| I\\(", term_labels))) {
    stop("DIGRAM model formulas allow only `:` interaction terms joined by `+`.", call. = FALSE)
  }

  rows <- lapply(term_labels, function(label) {
    vars <- strsplit(label, ":", fixed = TRUE)[[1L]]
    if (length(vars) != 2L) {
      stop("DIGRAM model terms must be two-way interactions: ", label, call. = FALSE)
    }
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

empty_score_effect_terms <- function() {
  data.frame(
    type = character(),
    exogenous = character(),
    source = character(),
    status = character(),
    stringsAsFactors = FALSE
  )
}

gRm_terms_formula <- function(terms) {
  if (length(terms) == 0L) {
    return(NULL)
  }
  stats::as.formula(
    paste("~", paste(terms, collapse = " + ")),
    env = parent.frame()
  )
}

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
    status = rows$stage %||% "selected",
    stringsAsFactors = FALSE
  )
}

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

#' @export
print.gRm_model <- function(x, ...) {
  label <- if (identical(x$model_type, "rasch")) "DIGRAM Rasch specification" else "DIGRAM GLLRM specification"
  cat(label, "\n", sep = "")
  cat("  items: ", paste(x$analysis$items, collapse = ", "), "\n", sep = "")
  cat("  LD terms: ", nrow(x$ld), "\n", sep = "")
  cat("  DIF terms: ", nrow(x$dif), "\n", sep = "")
  invisible(x)
}

#' @export
print.gRm_gllrm_spec <- function(x, ...) {
  print.gRm_model(x, ...)
}
