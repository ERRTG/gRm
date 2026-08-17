#' Internal summary table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @return The internal `summary_table()` computation result.
#' @keywords internal
#' @noRd
summary_table <- function(x) {
  UseMethod("summary_table")
}
#' @export
print.gRm_analysis <- function(x, ...) {
  reject_public_dots(...)
  data_name <- x$data_name %||% x$name %||% "<unnamed>"
  id <- x$id %||% "none"
  cat("gRm: Graphical Log-Linear Rasch Model analysis\n\n")
  cat("  Source: ", data_name, "\n", sep = "")
  cat("  Rows: ", nrow(x$project$raw_data), "\n", sep = "")
  cat("  ID: ", id, "\n", sep = "")
  cat("  Items: ", length(x$items), "\n", sep = "")
  cat("  Exogenous: ", length(x$exogenous), "\n", sep = "")
  cat("  Score groups: ", length(x$score_groups %||% integer()), "\n", sep = "")
  cat("\n")
  cat("Use summary(x) to show variable names, level counts, and score-group distribution.\n")
  invisible(x)
}

#' @export
print.gRm_model <- function(x, ...) {
  reject_public_dots(...)
  print(summary(x))
  invisible(x)
}

#' @export
print.gRm_fit <- function(x, ...) {
  reject_public_dots(...)
  values <- x$values %||% list()
  cat(public_model_label(x$spec, "fit"), "\n\n", sep = "")
  cat("  Converged: ", summary_scalar(x$convergence$converged %||% NA), "\n", sep = "")
  cat("  Iterations: ", summary_scalar(x$convergence$iterations %||% values$n_step %||% NA_integer_), "\n", sep = "")
  cat("  Delta: ", summary_scalar(x$convergence$report_delta %||% x$convergence$delta %||% NA_real_), "\n", sep = "")
  cat(
    "  Negative log likelihood (DIGRAM): ",
    summary_scalar(values$negative_log_likelihood %||% values$log_likelihood %||% NA_real_),
    "\n",
    sep = ""
  )
  if (!is.null(values$n_parameters)) {
    cat("  Parameters: ", summary_scalar(values$n_parameters), "\n", sep = "")
  }
  cat("\n")
  cat("Use summary(x) to show the fitted-model details.\n")
  invisible(x)
}

#' @export
print.gRm_direct_table <- function(x, ...) {
  reject_public_dots(...)
  display <- x
  class(display) <- "data.frame"
  title <- attr(x, "title", exact = TRUE)
  if (!is.null(title) && nzchar(title)) {
    cat(title, "\n\n", sep = "")
  }
  if (inherits(x, "gRm_item_fit") && identical(attr(x, "which", exact = TRUE), "tests")) {
    print_item_fit_tests_table(display)
  } else if (inherits(x, "gRm_item_fit") && identical(attr(x, "which", exact = TRUE), "items")) {
    print_item_fit_items_table(display)
  } else if (inherits(x, "gRm_score_effects")) {
    print_score_effects_tests_table(display)
  } else {
    print_summary_table(display)
  }
  print_summary_table_note(attr(x, "table_note", exact = TRUE))
  invisible(x)
}
#' @export
print.gRm_local_dependence <- function(x, ...) {
  reject_public_dots(...)
  print_diagnostic_status(
    x,
    title = "gRm: Local dependence tests",
    candidate_label = "Candidate pairs"
  )
}
#' @export
print.gRm_dif <- function(x, ...) {
  reject_public_dots(...)
  print_diagnostic_status(
    x,
    title = "gRm: DIF tests",
    candidate_label = "Candidate tests",
    before_selection = c(
      `Included DIF terms` = nrow(diagnostic_included_tests(x))
    ),
    after_convergence = c(
      `Output-stable fits` = diagnostic_true_count(diagnostic_tests(x), "output_stable")
    )
  )
}
#' @export
print.gRm_global_homogeneity <- function(x, ...) {
  reject_public_dots(...)
  summary <- x$values$summary %||% list()
  groups <- x$values$score_groups %||% data.frame()
  uniform_ld_n <- global_homogeneity_table_nrow(x$values$uniform_ld %||% data.frame())
  uniform_dif_n <- global_homogeneity_table_nrow(x$values$uniform_dif %||% data.frame())
  cat("gRm: Global homogeneity test\n\n")
  cat("  Score groups: ", summary_scalar(summary$n_groups %||% nrow(groups)), "\n", sep = "")
  cat("  Parameters: ", summary_scalar(summary$n_parameters %||% NA_integer_), "\n", sep = "")
  cat("  CLR: ", summary_scalar(summary$clr %||% NA_real_), "\n", sep = "")
  cat("  Df: ", summary_scalar(summary$df %||% NA_integer_), "\n", sep = "")
  cat("  Pr(>CLR): ", summary_p_value(summary$p_value %||% NA_real_), "\n", sep = "")
  cat("  Non-converged group fits: ", diagnostic_false_count(groups, "converged"), "\n", sep = "")
  if (uniform_ld_n > 0L) {
    cat("  Uniform LD tests: ", uniform_ld_n, "\n", sep = "")
  }
  if (uniform_dif_n > 0L) {
    cat("  Uniform DIF tests: ", uniform_dif_n, "\n", sep = "")
  }
  cat("\n")
  if (uniform_ld_n > 0L || uniform_dif_n > 0L) {
    cat("Use summary(x) to show the global test, score groups, item means, and uniform interaction tests.\n")
  } else {
    cat("Use summary(x) to show the global test, score groups, and item means.\n")
  }
  invisible(x)
}

#' Internal print diagnostic status helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @param title Internal `title` value used by this helper.
#' @param candidate_label Internal `candidate_label` value used by this helper.
#' @param before_selection Internal `before_selection` value used by this helper.
#' @param after_convergence Internal `after_convergence` value used by this helper.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
print_diagnostic_status <- function(x,
                                    title,
                                    candidate_label,
                                    before_selection = integer(),
                                    after_convergence = integer()) {
  tests <- diagnostic_tests(x)
  bh_critical_p <- x$values$bh_critical_p %||% NA_real_
  selected <- public_selected_by_bh(tests, bh_critical_p)
  cat(title, "\n\n", sep = "")
  cat("  ", candidate_label, ": ", nrow(tests), "\n", sep = "")
  for (name in names(before_selection)) {
    cat("  ", name, ": ", before_selection[[name]], "\n", sep = "")
  }
  cat("  Selected by BH FDR 0.05: ", nrow(selected), "\n", sep = "")
  cat("  Non-converged fits: ", diagnostic_false_count(tests, "converged"), "\n", sep = "")
  for (name in names(after_convergence)) {
    cat("  ", name, ": ", after_convergence[[name]], "\n", sep = "")
  }
  cat("  BH threshold: ", diagnostic_threshold_label(bh_critical_p), "\n\n", sep = "")
  cat("Use summary(x) to show the complete test table.\n")
  invisible(x)
}

#' Internal diagnostic tests helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @return The internal `diagnostic_tests()` computation result.
#' @keywords internal
#' @noRd
diagnostic_tests <- function(x) {
  tests <- x$values$tests %||% data.frame()
  if (is.data.frame(tests)) {
    tests
  } else {
    data.frame()
  }
}

#' Internal diagnostic included tests helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @return The internal `diagnostic_included_tests()` computation result.
#' @keywords internal
#' @noRd
diagnostic_included_tests <- function(x) {
  included <- x$values$included_tests %||% data.frame()
  if (is.data.frame(included)) {
    included
  } else {
    data.frame()
  }
}

#' Internal diagnostic true count helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param tests Diagnostic test rows.
#' @param column Internal `column` value used by this helper.
#' @return The internal `diagnostic_true_count()` computation result.
#' @keywords internal
#' @noRd
diagnostic_true_count <- function(tests, column) {
  if (!is.data.frame(tests) || !column %in% names(tests)) {
    return(0L)
  }
  sum(!is.na(tests[[column]]) & tests[[column]])
}

#' Internal diagnostic false count helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param tests Diagnostic test rows.
#' @param column Internal `column` value used by this helper.
#' @return The internal `diagnostic_false_count()` computation result.
#' @keywords internal
#' @noRd
diagnostic_false_count <- function(tests, column) {
  if (!is.data.frame(tests) || !column %in% names(tests)) {
    return(0L)
  }
  sum(!is.na(tests[[column]]) & !tests[[column]])
}

#' Internal diagnostic threshold label helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param threshold Internal `threshold` value used by this helper.
#' @return The internal `diagnostic_threshold_label()` computation result.
#' @keywords internal
#' @noRd
diagnostic_threshold_label <- function(threshold) {
  threshold <- threshold %||% NA_real_
  threshold <- threshold[[1L]]
  if (is.na(threshold)) {
    return("NA")
  }
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  format.pval(threshold, digits = dig_tst, eps = .Machine$double.eps)
}

#' Internal print gRm result helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
print_gRm_result <- function(x) {
  cat("<", class(x)[[1L]], ">\n", sep = "")
  cat("  items: ", paste(x$analysis$items, collapse = ", "), "\n", sep = "")
  if (length(x$analysis$exogenous)) {
    cat("  exogenous: ", paste(x$analysis$exogenous, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' Internal make gRm summary helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @param title Internal `title` value used by this helper.
#' @param which Internal `which` value used by this helper.
#' @param allowed Internal `allowed` value used by this helper.
#' @param tables Internal `tables` value used by this helper.
#' @param bh_tests Internal `bh_tests` value used by this helper.
#' @param header Internal `header` value used by this helper.
#' @param print_tables Internal `print_tables` value used by this helper.
#' @param table_names Internal `table_names` value used by this helper.
#' @param table_notes Internal `table_notes` value used by this helper.
#' @param extra_tables Internal `extra_tables` value used by this helper.
#' @param summary_attributes Internal `summary_attributes` value used by this helper.
#' @param print_table_names Internal `print_table_names` value used by this helper.
#' @param remark_tables Internal `remark_tables` value used by this helper.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
make_gRm_summary <- function(object,
                             title,
                             which,
                             allowed,
                             tables,
                             bh_tests = FALSE,
                             header = character(),
                             print_tables = TRUE,
                             table_names = NULL,
                             table_notes = list(),
                             extra_tables = list(),
                             summary_attributes = list(),
                             print_table_names = TRUE,
                             remark_tables = NULL) {
  new_gRm_summary(
    object,
    title = title,
    which = validate_summary_which(which, allowed),
    tables = tables,
    bh_tests = bh_tests,
    header = header,
    print_tables = print_tables,
    table_names = table_names,
    table_notes = table_notes,
    extra_tables = extra_tables,
    summary_attributes = summary_attributes,
    print_table_names = print_table_names,
    remark_tables = remark_tables
  )
}

#' @export
summary.gRm_analysis <- function(object, ...) {
  reject_summary_which(...)
  make_gRm_summary(
    object,
    title = "gRm: Graphical Log-Linear Rasch Model analysis",
    which = c("data", "score_groups"),
    allowed = c("data", "score_groups"),
    tables = list(
      data = analysis_summary_table(object),
      score_groups = analysis_score_group_distribution(object)
    ),
    header = analysis_summary_header(object),
    print_tables = FALSE
  )
}

#' @export
summary.gRm_model <- function(object, ...) {
  reject_summary_which(...)
  make_gRm_summary(
    object,
    title = public_model_label(object, "specification"),
    which = "model",
    allowed = "model",
    tables = list(
      model = model_summary_table(object),
      ld = object$ld,
      dif = object$dif
    ),
    table_names = c("model", "ld", "dif")
  )
}

#' @export
summary.gRm_fit <- function(object, which = NULL, ...) {
  reject_public_dots(...)
  param <- item_parameter_result_tables(object$values)
  tables <- list(
    parameters = param$item_statistics %||% data.frame(),
    thresholds = param$thresholds %||% data.frame(),
    terms = rbind_fill(object$spec$ld %||% data.frame(), object$spec$dif %||% data.frame())
  )
  fit <- fit_summary_table(object)
  requested <- validate_summary_which(which, c("parameters", "thresholds"))
  table_names <- if (is.null(which)) {
    c("parameters", "thresholds", "terms")
  } else {
    requested
  }
  new_gRm_summary(
    object,
    title = public_model_label(object$spec, "fit"),
    which = requested,
    tables = tables,
    table_names = table_names,
    extra_tables = list(fit = fit)
  )
}

#' @export
summary.gRm_screen <- function(object, ...) {
  reject_summary_which(...)
  tables <- public_screen_summary_tables(object)
  make_gRm_summary(
    object,
    title = "gRm: SCREEN J tests",
    which = c("local_dependence", "dif", "score_effects"),
    allowed = c("local_dependence", "dif", "score_effects"),
    header = screen_summary_header(object, tables),
    tables = list(
      local_dependence = tables$local_dependence,
      dif = tables$dif,
      score_effects = tables$score_effects
    ),
    table_notes = list(
      score_effects = public_screen_bh_marker_note(tables$bh)
    ),
    extra_tables = list(
      selected = tables$selected
    ),
    summary_attributes = list(
      bh = tables$bh,
      model_terms = tables$model_terms
    ),
    print_table_names = FALSE
  )
}

#' @export
summary.gRm_local_dependence <- function(object, ...) {
  reject_summary_which(...)
  tables <- public_value_tables(object$values)
  make_gRm_summary(
    object,
    title = "gRm: Local dependence tests",
    which = "tests",
    allowed = "tests",
    bh_tests = FALSE,
    tables = list(
      tests = public_local_dependence_tests(
        object$values$tests %||% data.frame(),
        object$values$bh_critical_p
      )
    ),
    table_notes = list(
      tests = public_bh_marker_note(tables$bh_thresholds %||% data.frame())
    ),
    extra_tables = list(
      selected = tables$selected %||% public_selected_by_bh(object$values$tests, object$values$bh_critical_p)
    ),
    summary_attributes = list(
      bh = tables$bh_thresholds %||% data.frame()
    ),
    print_table_names = FALSE,
    remark_tables = list(
      tests = object$values$tests %||% data.frame()
    )
  )
}

#' @export
summary.gRm_dif <- function(object, ...) {
  reject_summary_which(...)
  tables <- public_value_tables(object$values)
  make_gRm_summary(
    object,
    title = "gRm: DIF tests",
    which = "tests",
    allowed = "tests",
    bh_tests = FALSE,
    tables = list(
      tests = public_dif_tests(
        object$values$tests %||% data.frame(),
        object$values$bh_critical_p
      )
    ),
    table_notes = list(
      tests = public_bh_marker_note(tables$bh_thresholds %||% data.frame())
    ),
    extra_tables = list(
      selected = tables$selected %||% public_selected_by_bh(object$values$tests, object$values$bh_critical_p),
      included = object$values$included_tests %||% data.frame()
    ),
    summary_attributes = list(
      bh = tables$bh_thresholds %||% data.frame()
    ),
    print_table_names = FALSE,
    remark_tables = list(
      tests = object$values$tests %||% data.frame()
    )
  )
}

#' @export
summary.gRm_global_homogeneity <- function(object, which = NULL, ...) {
  reject_public_dots(...)
  source_status <- global_homogeneity_source_status(object$values$items %||% data.frame())
  values <- object$values %||% list()
  default_which <- global_homogeneity_default_summary_sections(values)
  allowed <- c("test", "score_groups", "item_means", "uniform_ld", "uniform_dif")
  requested <- if (is.null(which)) {
    default_which
  } else {
    validate_summary_which(which, allowed)
  }
  table_names <- if (is.null(which)) {
    default_which
  } else {
    requested
  }
  make_gRm_summary(
    object,
    title = "gRm: Global homogeneity test",
    which = requested,
    allowed = allowed,
    tables = list(
      test = public_global_homogeneity_test(values$summary %||% list()),
      score_groups = public_global_homogeneity_score_groups(values$score_groups %||% data.frame()),
      item_means = public_global_homogeneity_item_means(
        values$items %||% data.frame(),
        values$score_groups %||% data.frame()
      ),
      uniform_ld = public_global_homogeneity_uniform_ld(
        values$uniform_ld %||% data.frame(),
        values
      ),
      uniform_dif = public_global_homogeneity_uniform_dif(
        values$uniform_dif %||% data.frame(),
        values
      )
    ),
    table_names = table_names,
    table_notes = list(
      item_means = "Note: DIGRAM item residual and marker cells are not source-backed in gRm and are not printed."
    ),
    extra_tables = list(
      source_status = source_status
    ),
    summary_attributes = list(
      source_status = source_status
    ),
    print_table_names = FALSE,
    remark_tables = list(
      score_groups = values$score_groups %||% data.frame()
    )
  )
}

#' Internal global homogeneity source status helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param items Item selection or item metadata.
#' @return The internal `global_homogeneity_source_status()` computation result.
#' @keywords internal
#' @noRd
global_homogeneity_source_status <- function(items) {
  residual <- runtime_source_backed_summary(items, "residual_runtime_source_backed")
  marker <- runtime_source_backed_summary(items, "marker_runtime_source_backed")
  data.frame(
    item_residual_runtime_source_backed = residual$source_backed,
    item_residual_status = residual$status,
    item_marker_runtime_source_backed = marker$source_backed,
    item_marker_status = marker$status,
    stringsAsFactors = FALSE
  )
}

#' Internal runtime source backed summary helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param items Item selection or item metadata.
#' @param column Internal `column` value used by this helper.
#' @return The internal `runtime_source_backed_summary()` computation result.
#' @keywords internal
#' @noRd
runtime_source_backed_summary <- function(items, column) {
  if (!is.data.frame(items) || !column %in% names(items) || length(items[[column]]) == 0L) {
    return(list(source_backed = NA, status = "not_available"))
  }
  known <- items[[column]][!is.na(items[[column]])]
  if (length(known) == 0L) {
    return(list(source_backed = NA, status = "not_available"))
  }
  if (all(known)) {
    return(list(source_backed = TRUE, status = "source_backed"))
  }
  if (any(known)) {
    return(list(source_backed = FALSE, status = "partially_source_backed"))
  }
  list(source_backed = FALSE, status = "not_source_backed")
}

#' Internal new gRm summary helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @param title Internal `title` value used by this helper.
#' @param which Internal `which` value used by this helper.
#' @param tables Internal `tables` value used by this helper.
#' @param bh_tests Internal `bh_tests` value used by this helper.
#' @param header Internal `header` value used by this helper.
#' @param print_tables Internal `print_tables` value used by this helper.
#' @param table_names Internal `table_names` value used by this helper.
#' @param table_notes Internal `table_notes` value used by this helper.
#' @param extra_tables Internal `extra_tables` value used by this helper.
#' @param summary_attributes Internal `summary_attributes` value used by this helper.
#' @param print_table_names Internal `print_table_names` value used by this helper.
#' @param remark_tables Internal `remark_tables` value used by this helper.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
new_gRm_summary <- function(object,
                            title,
                            which,
                            tables,
                            bh_tests = FALSE,
                            header = character(),
                            print_tables = TRUE,
                            table_names = NULL,
                            table_notes = list(),
                            extra_tables = list(),
                            summary_attributes = list(),
                            print_table_names = TRUE,
                            remark_tables = NULL) {
  table_names <- table_names %||% which
  tables <- tables[table_names]
  remarks <- summary_convergence_remarks(remark_tables %||% tables)
  out <- c(list(
    title = title,
    which = which,
    tables = tables,
    header = header,
    remarks = remarks,
    bh_tests = isTRUE(bh_tests),
    print_tables = isTRUE(print_tables),
    print_table_names = isTRUE(print_table_names),
    table_notes = table_notes,
    object_class = class(object)
  ), tables, extra_tables)
  class(out) <- c(paste0("summary.", class(object)[[1L]]), "summary.gRm", "list")
  for (name in names(summary_attributes)) {
    attr(out, name) <- summary_attributes[[name]]
  }
  out
}
