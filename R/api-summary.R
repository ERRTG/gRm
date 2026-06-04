summary_table <- function(x) {
  UseMethod("summary_table")
}

#' @export
print.gRm_analysis <- function(x, ...) {
  cat("<gRm_analysis>\n")
  cat("  rows: ", nrow(x$project$raw_data), "\n", sep = "")
  cat("  items: ", paste(x$items, collapse = ", "), "\n", sep = "")
  cat("  exogenous: ", paste(x$exogenous, collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' @export
print.gRm_model <- function(x, ...) {
  cat(public_model_label(x, "specification"), "\n", sep = "")
  cat("  items: ", paste(x$analysis$items, collapse = ", "), "\n", sep = "")
  cat("  LD terms: ", nrow(x$ld), "\n", sep = "")
  cat("  DIF terms: ", nrow(x$dif), "\n", sep = "")
  invisible(x)
}

#' @export
print.gRm_fit <- function(x, ...) {
  cat(public_model_label(x$spec, "fit"), "\n", sep = "")
  cat("  converged: ", as.character(x$convergence$converged %||% NA), "\n", sep = "")
  invisible(x)
}

#' @export
print.gRm_score_effects <- function(x, ...) print_gRm_result(x)
#' @export
print.gRm_item_parameters <- function(x, ...) print_gRm_result(x)
#' @export
print.gRm_item_fit <- function(x, ...) print_gRm_result(x)
#' @export
print.gRm_local_dependence <- function(x, ...) print_gRm_result(x)
#' @export
print.gRm_dif <- function(x, ...) print_gRm_result(x)
#' @export
print.gRm_global_homogeneity <- function(x, ...) print_gRm_result(x)

print_gRm_result <- function(x) {
  cat("<", class(x)[[1L]], ">\n", sep = "")
  cat("  items: ", paste(x$analysis$items, collapse = ", "), "\n", sep = "")
  if (length(x$analysis$exogenous)) {
    cat("  exogenous: ", paste(x$analysis$exogenous, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.gRm_analysis <- function(object, which = "data", ...) {
  new_gRm_summary(
    object,
    title = "DIGRAM analysis",
    which = validate_summary_which(which, "data"),
    tables = list(data = analysis_summary_table(object))
  )
}

#' @export
summary.gRm_item_analysis <- summary.gRm_analysis

#' @export
summary.gRm_model <- function(object, which = c("model", "ld", "dif"), ...) {
  which <- validate_summary_which(which, c("model", "ld", "dif"))
  new_gRm_summary(
    object,
    title = public_model_label(object, "specification"),
    which = which,
    tables = list(
      model = model_summary_table(object),
      ld = object$ld,
      dif = object$dif
    )
  )
}

#' @export
summary.gRm_gllrm_spec <- summary.gRm_model

#' @export
summary.gRm_fit <- function(object, which = c("fit", "parameters", "terms"), ...) {
  which <- validate_summary_which(which, c("fit", "parameters", "terms"))
  param <- item_parameter_detail_tables(object$values)
  new_gRm_summary(
    object,
    title = public_model_label(object$spec, "fit"),
    which = which,
    tables = list(
      fit = fit_summary_table(object),
      parameters = param$item_statistics %||% data.frame(),
      terms = rbind_fill(object$spec$ld %||% data.frame(), object$spec$dif %||% data.frame())
    )
  )
}

#' @export
summary.gRm_gllrm_fit <- summary.gRm_fit

#' @export
summary.gRm_screen <- function(object, which = c("tests", "selected", "all", "score_effects", "bh"), ...) {
  which <- validate_summary_which(which, c("tests", "selected", "all", "score_effects", "bh"))
  terms <- public_screen_terms(object)
  all_terms <- rbind_fill(terms$ld, terms$dif)
  selected_terms <- public_selected_rows(all_terms)
  new_gRm_summary(
    object,
    title = "DIGRAM screen",
    which = which,
    tables = list(
      tests = data.frame(
        selection = "Benjamini-Hochberg",
        n_ld = nrow(terms$ld),
        n_dif = nrow(terms$dif),
        n_score_effects = nrow(terms$score_effects),
        stringsAsFactors = FALSE
      ),
      selected = selected_terms,
      all = all_terms,
      score_effects = terms$score_effects,
      bh = public_screen_bh_table(object)
    )
  )
}

#' @export
summary.gRm_score_effects <- function(object, which = c("selected", "tests", "bh"), ...) {
  which <- validate_summary_which(which, c("selected", "tests", "bh"))
  tables <- public_value_tables(object$values)
  new_gRm_summary(
    object,
    title = "DIGRAM score effects",
    which = which,
    tables = list(
      selected = tables$score_effect_selected %||% data.frame(),
      tests = tables$score_effect_tests %||% data.frame(),
      bh = tables$score_effect_bh_thresholds %||% tables$bh_thresholds %||% data.frame()
    )
  )
}

#' @export
summary.gRm_item_parameters <- function(object, which = c("tests", "coefficients", "items", "thresholds", "fit"), ...) {
  which <- validate_summary_which(which, c("tests", "coefficients", "items", "thresholds", "fit"))
  tables <- item_parameter_detail_tables(object$values)
  coefficients <- tables$item_statistics %||% data.frame()
  new_gRm_summary(
    object,
    title = "DIGRAM item parameters",
    which = which,
    tables = list(
      tests = coefficients,
      coefficients = coefficients,
      items = coefficients,
      thresholds = tables$thresholds %||% data.frame(),
      fit = tables$fit_summary %||% data.frame()
    )
  )
}

item_fit_summary_items_table <- function(object, public_tables) {
  public_tables$item_fit_summaries %||%
    object$values$extended$summaries %||%
    data.frame()
}

#' @export
summary.gRm_item_fit <- function(object, which = c("tests", "items", "bh"), ...) {
  which <- validate_summary_which(which, c("tests", "items", "bh"))
  tables <- public_value_tables(object$values)
  item_tests <- tables$statistics %||% object$values$items %||% data.frame()
  item_summaries <- item_fit_summary_items_table(object, tables)
  new_gRm_summary(
    object,
    title = "DIGRAM item fit",
    which = which,
    tables = list(
      tests = item_tests,
      items = item_summaries,
      bh = tables$bh_thresholds %||% data.frame()
    )
  )
}

#' @export
summary.gRm_local_dependence <- function(object, which = c("selected", "tests", "bh"), ...) {
  which <- validate_summary_which(which, c("selected", "tests", "bh"))
  tables <- public_value_tables(object$values)
  new_gRm_summary(
    object,
    title = "DIGRAM local dependence",
    which = which,
    tables = list(
      selected = tables$selected %||% public_selected_by_bh(object$values$tests, object$values$bh_critical_p),
      tests = object$values$tests %||% data.frame(),
      bh = tables$bh_thresholds %||% data.frame()
    )
  )
}

#' @export
summary.gRm_dif <- function(object, which = c("selected", "tests", "active", "bh"), ...) {
  which <- validate_summary_which(which, c("selected", "tests", "active", "bh"))
  tables <- public_value_tables(object$values)
  new_gRm_summary(
    object,
    title = "DIGRAM DIF",
    which = which,
    tables = list(
      selected = tables$selected %||% public_selected_by_bh(object$values$tests, object$values$bh_critical_p),
      tests = object$values$tests %||% data.frame(),
      active = object$values$active_tests %||% data.frame(),
      bh = tables$bh_thresholds %||% data.frame()
    )
  )
}

#' @export
summary.gRm_global_homogeneity <- function(object, which = c("tests", "summary", "groups", "items"), ...) {
  which <- validate_summary_which(which, c("tests", "summary", "groups", "items"))
  gh_summary <- list_to_one_row(object$values$summary %||% list())
  new_gRm_summary(
    object,
    title = "DIGRAM global homogeneity",
    which = which,
    tables = list(
      tests = gh_summary,
      summary = gh_summary,
      groups = object$values$score_groups %||% data.frame(),
      items = object$values$items %||% data.frame()
    )
  )
}

new_gRm_summary <- function(object, title, which, tables) {
  tables <- tables[which]
  out <- c(list(
    title = title,
    which = which,
    tables = tables,
    object_class = class(object)
  ), tables)
  class(out) <- c(paste0("summary.", class(object)[[1L]]), "summary.gRm", "list")
  out
}

#' @export
print.summary.gRm <- function(x, ...) {
  cat(x$title, "\n", sep = "")
  if ("tests" %in% names(x$tables)) {
    cat("Benjamini-Hochberg\n")
  }
  for (name in names(x$tables)) {
    table <- x$tables[[name]]
    cat("\n", name, "\n", sep = "")
    if (is.data.frame(table) && nrow(table)) {
      print(public_format_table(table), row.names = FALSE)
    } else {
      cat("  <none>\n")
    }
  }
  invisible(x)
}

validate_summary_which <- function(which, allowed) {
  if (missing(which) || is.null(which)) {
    return(allowed)
  }
  if (!is.character(which) || anyNA(which) || !length(which)) {
    stop("`which` must name one or more summary sections.", call. = FALSE)
  }
  unknown <- setdiff(which, allowed)
  if (length(unknown)) {
    stop(
      "Unknown `which` summary section",
      if (length(unknown) > 1L) "s" else "",
      ": ",
      paste(unknown, collapse = ", "),
      ". Available sections: ",
      paste(allowed, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  unique(which)
}

analysis_summary_table <- function(object) {
  data.frame(
    n_rows = nrow(object$project$raw_data),
    n_items = length(object$items),
    n_exogenous = length(object$exogenous),
    items = paste(object$items, collapse = ", "),
    exogenous = paste(object$exogenous, collapse = ", "),
    score_groups = paste(object$score_groups %||% integer(), collapse = ", "),
    stringsAsFactors = FALSE
  )
}

model_summary_table <- function(object) {
  data.frame(
    model_type = public_model_type(object),
    n_items = length(object$analysis$items),
    n_exogenous = length(object$analysis$exogenous),
    n_ld = nrow(object$ld %||% data.frame()),
    n_dif = nrow(object$dif %||% data.frame()),
    stringsAsFactors = FALSE
  )
}

fit_summary_table <- function(object) {
  values <- object$values %||% list()
  data.frame(
    model_type = public_model_type(object$spec),
    log_likelihood = values$log_likelihood %||% NA_real_,
    n_parameters = values$n_parameters %||% NA_integer_,
    likelihood_n = values$likelihood_n %||% NA_integer_,
    converged = object$convergence$converged %||% NA,
    iterations = object$convergence$iterations %||% values$n_step %||% NA_integer_,
    delta = object$convergence$delta %||% values$delta %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

public_model_type <- function(object) {
  object$model_type %||% if (
    nrow(object$ld %||% data.frame()) > 0L ||
      nrow(object$dif %||% data.frame()) > 0L
  ) {
    "gllrm"
  } else {
    "rasch"
  }
}

public_model_label <- function(object, noun) {
  paste("DIGRAM", if (identical(public_model_type(object), "rasch")) "Rasch" else "GLLRM", noun)
}

public_screen_terms <- function(object) {
  if (!is.null(object$terms)) {
    return(object$terms)
  }
  model_terms(object)
}

public_screen_bh_table <- function(object) {
  bh <- object$values$bh %||% list()
  data.frame(
    fdr = c("0.05", "0.01", "0.001"),
    p_value = c(bh$fdr_05 %||% NA_real_, bh$fdr_01 %||% NA_real_, bh$fdr_001 %||% NA_real_),
    n_tests = bh$n_tests %||% NA_integer_,
    stringsAsFactors = FALSE
  )
}

public_value_tables <- function(values) {
  if (inherits(values, "gRm_item_parameters_values") || inherits(values, "gRm_active_gllrm_values")) {
    return(item_parameter_detail_tables(values))
  }
  if (exists("details", mode = "function")) {
    return(details(values)$tables %||% list())
  }
  list()
}

public_selected_by_bh <- function(tests, threshold) {
  if (!is.data.frame(tests) || !"p_value" %in% names(tests) || is.na(threshold)) {
    return(data.frame())
  }
  out <- tests[tests$p_value <= threshold, , drop = FALSE]
  rownames(out) <- NULL
  out
}

public_selected_rows <- function(rows) {
  if (!is.data.frame(rows) || !nrow(rows)) {
    return(data.frame())
  }
  if ("selected" %in% names(rows)) {
    rows <- rows[rows$selected %in% TRUE, , drop = FALSE]
  } else if ("status" %in% names(rows)) {
    rows <- rows[rows$status %in% "selected", , drop = FALSE]
  }
  rownames(rows) <- NULL
  rows
}

public_format_table <- function(table) {
  out <- table
  numeric_cols <- vapply(out, is.numeric, logical(1L))
  out[numeric_cols] <- lapply(out[numeric_cols], function(x) signif(x, 6L))
  out
}
