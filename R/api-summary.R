summary_table <- function(x) {
  UseMethod("summary_table")
}

#' @export
print.gRm_analysis <- function(x, ...) {
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
  print(summary(x))
  invisible(x)
}

#' @export
print.gRm_fit <- function(x, ...) {
  values <- x$values %||% list()
  cat(public_model_label(x$spec, "fit"), "\n\n", sep = "")
  cat("  Converged: ", summary_scalar(x$convergence$converged %||% NA), "\n", sep = "")
  cat("  Iterations: ", summary_scalar(x$convergence$iterations %||% values$n_step %||% NA_integer_), "\n", sep = "")
  cat("  Delta: ", summary_scalar(x$convergence$delta %||% values$delta %||% NA_real_), "\n", sep = "")
  cat("  Log likelihood: ", summary_scalar(values$log_likelihood %||% NA_real_), "\n", sep = "")
  if (!is.null(values$n_parameters)) {
    cat("  Parameters: ", summary_scalar(values$n_parameters), "\n", sep = "")
  }
  cat("\n")
  cat("Use summary(x) to show the fitted-model details.\n")
  invisible(x)
}

#' @export
print.gRm_direct_table <- function(x, ...) {
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
  print_diagnostic_status(
    x,
    title = "gRm: Local dependence tests",
    candidate_label = "Candidate pairs"
  )
}
#' @export
print.gRm_dif <- function(x, ...) {
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

diagnostic_tests <- function(x) {
  tests <- x$values$tests %||% data.frame()
  if (is.data.frame(tests)) {
    tests
  } else {
    data.frame()
  }
}

diagnostic_included_tests <- function(x) {
  included <- x$values$included_tests %||% data.frame()
  if (is.data.frame(included)) {
    included
  } else {
    data.frame()
  }
}

diagnostic_true_count <- function(tests, column) {
  if (!is.data.frame(tests) || !column %in% names(tests)) {
    return(0L)
  }
  sum(!is.na(tests[[column]]) & tests[[column]])
}

diagnostic_false_count <- function(tests, column) {
  if (!is.data.frame(tests) || !column %in% names(tests)) {
    return(0L)
  }
  sum(!is.na(tests[[column]]) & !tests[[column]])
}

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

print_gRm_result <- function(x) {
  cat("<", class(x)[[1L]], ">\n", sep = "")
  cat("  items: ", paste(x$analysis$items, collapse = ", "), "\n", sep = "")
  if (length(x$analysis$exogenous)) {
    cat("  exogenous: ", paste(x$analysis$exogenous, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

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

#' @export
print.summary.gRm <- function(x, ...) {
  cat(x$title, "\n", sep = "")
  if (length(x$header)) {
    cat(paste(x$header, collapse = "\n"), "\n", sep = "")
  }
  if (isTRUE(x$bh_tests) && "tests" %in% names(x$tables)) {
    cat("Benjamini-Hochberg\n")
  }
  if (is.data.frame(x$remarks) && nrow(x$remarks)) {
    cat("\nRemarks\n")
    for (row in seq_len(nrow(x$remarks))) {
      cat("  ", x$remarks$message[[row]], "\n", sep = "")
    }
  }
  if (inherits(x, "summary.gRm_analysis") && "score_groups" %in% names(x$tables)) {
    print_analysis_score_group_distribution(x$tables$score_groups)
  }
  if (!isTRUE(x$print_tables)) {
    return(invisible(x))
  }
  for (name in names(x$tables)) {
    table <- x$tables[[name]]
    if (inherits(x, "summary.gRm_model")) {
      if (identical(name, "model")) {
        cat("\nModel\n")
        print_model_summary_table(table)
        next
      }
      if (identical(name, "ld")) {
        print_model_term_table("Local dependence terms", table, "ld")
        next
      }
      if (identical(name, "dif")) {
        print_model_term_table("DIF terms", table, "dif")
        next
      }
    }
    if (inherits(x, "summary.gRm_fit")) {
      if (identical(name, "parameters")) {
        cat("\nItem parameters\n")
        print_fit_parameter_table(table)
        next
      }
      if (identical(name, "thresholds")) {
        cat("\nThresholds\n")
        print_fit_threshold_table(table)
        next
      }
      if (identical(name, "terms")) {
        cat("\nModel terms\n")
        print_fit_term_table(table)
        next
      }
    }
    if (inherits(x, "summary.gRm_screen")) {
      if (identical(name, "local_dependence")) {
        cat("\nLocal dependence\n")
      } else if (identical(name, "dif")) {
        cat("\nDIF\n")
      } else if (identical(name, "score_effects")) {
        cat("\nScore effects\n")
      } else if (isTRUE(x$print_table_names)) {
        cat("\n", name, "\n", sep = "")
      } else {
        cat("\n")
      }
      print_screen_tests_table(table)
      print_summary_table_note(x$table_notes[[name]])
      next
    }
    if (isTRUE(x$print_table_names)) {
      cat("\n", name, "\n", sep = "")
    } else {
      cat("\n")
    }
    if (inherits(x, "summary.gRm_local_dependence") && identical(name, "tests")) {
      print_diagnostic_tests_table(table, gamma_columns = "WPG")
    } else if (inherits(x, "summary.gRm_dif") && identical(name, "tests")) {
      print_diagnostic_tests_table(table, gamma_columns = "Gamma")
    } else if (inherits(x, "summary.gRm_m2") || inherits(x, "summary.gRm_m3")) {
      print_diagnostic_tests_table(table)
    } else if (inherits(x, "summary.gRm_global_homogeneity")) {
      if (identical(name, "test")) {
        cat("Global test\n")
        print_global_homogeneity_test_table(table)
      } else if (identical(name, "score_groups")) {
        cat("Score groups\n")
        print_global_homogeneity_score_group_table(table)
      } else if (identical(name, "item_means")) {
        cat("Item means\n")
        print_global_homogeneity_item_means_table(table)
      } else if (identical(name, "uniform_ld")) {
        cat("Uniform local dependence\n")
        print_global_homogeneity_uniform_table(table)
      } else if (identical(name, "uniform_dif")) {
        cat("Uniform DIF\n")
        print_global_homogeneity_uniform_table(table)
      } else {
        print_summary_table(table)
      }
    } else {
      print_summary_table(table)
    }
    print_summary_table_note(x$table_notes[[name]])
  }
  invisible(x)
}

print_summary_table_note <- function(note) {
  if (!length(note) || all(is.na(note)) || !nzchar(note[[1L]])) {
    return(invisible(NULL))
  }
  cat("---\n")
  cat(note[[1L]], "\n", sep = "")
  invisible(NULL)
}

print_summary_table <- function(table) {
  if (is.data.frame(table) && nrow(table)) {
    print(public_format_table(table), row.names = FALSE)
  } else {
    cat("  <none>\n")
  }
  invisible(NULL)
}

print_diagnostic_tests_table <- function(table, gamma_columns = character()) {
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  <none>\n")
    return(invisible(NULL))
  }
  display <- table
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  if ("Chisq" %in% names(display)) {
    display$Chisq <- format(signif(display$Chisq, dig_tst), digits = dig_tst)
  }
  if ("Pr(>Chisq)" %in% names(display)) {
    display[["Pr(>Chisq)"]] <- format.pval(
      display[["Pr(>Chisq)"]],
      digits = dig_tst,
      eps = .Machine$double.eps
    )
  }
  for (column in gamma_columns) {
    if (column %in% names(display)) {
      display[[column]] <- format(signif(display[[column]], digits), digits = digits)
    }
  }
  if ("delta" %in% names(display)) {
    display$delta <- format(signif(display$delta, digits), digits = digits)
  }
  print(
    display,
    row.names = FALSE,
    quote = FALSE,
    right = TRUE,
    na.print = "NA",
    max = diagnostic_print_max(display)
  )
  invisible(NULL)
}

print_screen_tests_table <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  <none>\n")
    return(invisible(NULL))
  }
  display <- table
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  statistic_columns <- c("Chisq", "Gamma", "WPG")
  p_value_columns <- c("Pr(>Chisq)", "Pr(>|Gamma|)", "Pr(>|WPG|)")
  for (column in statistic_columns) {
    if (column %in% names(display)) {
      missing <- is.na(display[[column]])
      display[[column]] <- format(signif(display[[column]], dig_tst), digits = dig_tst)
      display[[column]][missing] <- ""
    }
  }
  for (column in p_value_columns) {
    if (column %in% names(display)) {
      missing <- is.na(display[[column]])
      display[[column]] <- format.pval(
        display[[column]],
        digits = dig_tst,
        eps = .Machine$double.eps
      )
      display[[column]][missing] <- ""
    }
  }
  if ("Df" %in% names(display)) {
    missing <- is.na(display$Df)
    display$Df <- as.character(display$Df)
    display$Df[missing] <- ""
  }
  print(
    display,
    row.names = FALSE,
    quote = FALSE,
    right = TRUE,
    na.print = "",
    max = diagnostic_print_max(display)
  )
  invisible(NULL)
}

print_global_homogeneity_test_table <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  <none>\n")
    return(invisible(NULL))
  }
  display <- table
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  statistic_columns <- c("LogLik full", "LogLik groups", "CLR")
  for (column in statistic_columns) {
    if (column %in% names(display)) {
      display[[column]] <- format(signif(display[[column]], dig_tst), digits = dig_tst)
    }
  }
  if ("Pr(>CLR)" %in% names(display)) {
    display[["Pr(>CLR)"]] <- format.pval(
      display[["Pr(>CLR)"]],
      digits = dig_tst,
      eps = .Machine$double.eps
    )
  }
  print(
    display,
    row.names = FALSE,
    quote = FALSE,
    right = TRUE,
    na.print = "NA"
  )
  invisible(NULL)
}

print_global_homogeneity_score_group_table <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  <none>\n")
    return(invisible(NULL))
  }
  display <- table
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  if ("LogLik" %in% names(display)) {
    display$LogLik <- format(signif(display$LogLik, dig_tst), digits = dig_tst)
  }
  if ("delta" %in% names(display)) {
    display$delta <- format(signif(display$delta, digits), digits = digits)
  }
  print(
    display,
    row.names = FALSE,
    quote = FALSE,
    right = TRUE,
    na.print = "NA"
  )
  invisible(NULL)
}

print_global_homogeneity_item_means_table <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  <none>\n")
    return(invisible(NULL))
  }
  display <- table
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  statistic_columns <- c("Observed mean", "Expected mean")
  for (column in statistic_columns) {
    if (column %in% names(display)) {
      display[[column]] <- format(signif(display[[column]], dig_tst), digits = dig_tst)
    }
  }
  print(
    display,
    row.names = FALSE,
    quote = FALSE,
    right = TRUE,
    na.print = "NA"
  )
  invisible(NULL)
}

print_score_effects_tests_table <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  <none>\n")
    return(invisible(NULL))
  }
  display <- table
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  statistic_columns <- c("Chisq", "Gamma")
  p_value_columns <- c(
    "Pr(>Chisq)",
    "Pr(Gamma+)",
    "Pr(|Gamma|)",
    "Exact Pr(>Chisq)",
    "Exact Pr(Gamma+)",
    "Exact Pr(|Gamma|)"
  )
  for (column in statistic_columns) {
    if (column %in% names(display)) {
      display[[column]] <- format(signif(display[[column]], dig_tst), digits = dig_tst)
    }
  }
  for (column in p_value_columns) {
    if (column %in% names(display)) {
      display[[column]] <- format.pval(
        display[[column]],
        digits = dig_tst,
        eps = .Machine$double.eps
      )
    }
  }
  print(
    display,
    row.names = FALSE,
    quote = FALSE,
    right = TRUE,
    na.print = "NA",
    max = diagnostic_print_max(display)
  )
  invisible(NULL)
}

print_global_homogeneity_uniform_table <- function(table) {
  gamma_columns <- grep("^(Obs|Exp) gamma ", names(table), value = TRUE)
  print_diagnostic_tests_table(table, gamma_columns = gamma_columns)
}

print_item_fit_tests_table <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  <none>\n")
    return(invisible(NULL))
  }
  display <- table
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  statistic_columns <- c("Outfit", "Outfit SE", "Infit", "Infit SE")
  p_value_columns <- c("Pr(>Outfit)", "Pr(>Infit)", "Pr(>Gamma)")
  gamma_columns <- c("Observed gamma", "Expected gamma", "Gamma SE")
  for (column in statistic_columns) {
    if (column %in% names(display)) {
      display[[column]] <- format(signif(display[[column]], dig_tst), digits = dig_tst)
    }
  }
  for (column in p_value_columns) {
    if (column %in% names(display)) {
      display[[column]] <- format.pval(
        display[[column]],
        digits = dig_tst,
        eps = .Machine$double.eps
      )
    }
  }
  for (column in gamma_columns) {
    if (column %in% names(display)) {
      display[[column]] <- format(signif(display[[column]], digits), digits = digits)
    }
  }
  print(
    display,
    row.names = FALSE,
    quote = FALSE,
    right = TRUE,
    na.print = "NA",
    max = diagnostic_print_max(display)
  )
  invisible(NULL)
}

print_item_fit_items_table <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  <none>\n")
    return(invisible(NULL))
  }
  display <- table
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  statistic_columns <- c(
    "Outfit observed",
    "Outfit expected",
    "Outfit total",
    "Infit observed",
    "Infit expected",
    "Infit variance",
    "Infit ratio"
  )
  for (column in statistic_columns) {
    if (column %in% names(display)) {
      display[[column]] <- format(signif(display[[column]], dig_tst), digits = dig_tst)
    }
  }
  print(
    display,
    row.names = FALSE,
    quote = FALSE,
    right = TRUE,
    na.print = "NA",
    max = diagnostic_print_max(display)
  )
  invisible(NULL)
}

diagnostic_print_max <- function(table) {
  cells <- as.double(nrow(table)) * as.double(ncol(table)) + 1
  as.integer(min(cells, .Machine$integer.max))
}

print_model_summary_table <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  <none>\n")
    return(invisible(NULL))
  }
  row <- table[1L, , drop = FALSE]
  cat("  Model type: ", public_model_type_label(row$model_type[[1L]]), "\n", sep = "")
  cat("  Items: ", row$n_items[[1L]], "\n", sep = "")
  cat("  Exogenous variables: ", row$n_exogenous[[1L]], "\n", sep = "")
  cat("  Local dependence terms: ", row$n_ld[[1L]], "\n", sep = "")
  cat("  DIF terms: ", row$n_dif[[1L]], "\n", sep = "")
  invisible(NULL)
}

print_fit_summary_table <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  <none>\n")
    return(invisible(NULL))
  }
  row <- table[1L, , drop = FALSE]
  cat("  Model type: ", public_model_type_label(row$model_type[[1L]]), "\n", sep = "")
  cat("  Converged: ", summary_scalar(row$converged[[1L]]), "\n", sep = "")
  cat("  Iterations: ", summary_scalar(row$iterations[[1L]]), "\n", sep = "")
  cat("  Delta: ", summary_scalar(row$delta[[1L]]), "\n", sep = "")
  cat("  Log likelihood: ", summary_scalar(row$log_likelihood[[1L]]), "\n", sep = "")
  cat("  Parameters: ", summary_scalar(row$n_parameters[[1L]]), "\n", sep = "")
  cat("  Likelihood rows: ", summary_scalar(row$likelihood_n[[1L]]), "\n", sep = "")
  invisible(NULL)
}

summary_scalar <- function(x) {
  if (length(x) == 0L || is.na(x)) {
    return("NA")
  }
  if (is.logical(x)) {
    return(if (isTRUE(x)) "yes" else "no")
  }
  if (is.numeric(x)) {
    return(as.character(signif(x, 6L)))
  }
  as.character(x)
}

summary_p_value <- function(x) {
  if (length(x) == 0L || is.na(x)) {
    return("NA")
  }
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  format.pval(x, digits = dig_tst, eps = .Machine$double.eps)
}

print_fit_parameter_table <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  <none>\n")
    return(invisible(NULL))
  }
  display <- public_format_table(table)
  names(display) <- fit_parameter_display_names(names(display))
  print(display, row.names = FALSE)
  invisible(NULL)
}

print_fit_threshold_table <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  <none>\n")
    return(invisible(NULL))
  }
  display <- table
  if ("score" %in% names(display)) {
    display$step <- paste(display$score - 1L, display$score, sep = " -> ")
    preferred <- c("item", "step", "score", "threshold")
    display <- display[c(intersect(preferred, names(display)), setdiff(names(display), preferred))]
  }
  display <- public_format_table(display)
  names(display) <- fit_threshold_display_names(names(display))
  print(display, row.names = FALSE)
  invisible(NULL)
}

fit_threshold_display_names <- function(names) {
  replacements <- c(
    item = "Item",
    step = "Step",
    score = "Score",
    threshold = "Threshold"
  )
  out <- replacements[names]
  missing <- is.na(out)
  out[missing] <- names[missing]
  unname(out)
}

fit_parameter_display_names <- function(names) {
  replacements <- c(
    item = "Item",
    location = "Location",
    midpoint = "Midpoint",
    target = "Target score",
    info_at_target = "Info at target",
    info_per_step = "Info per step"
  )
  out <- replacements[names]
  missing <- is.na(out)
  out[missing] <- names[missing]
  unname(out)
}

print_model_term_table <- function(title, table, type) {
  cat("\n", title, "\n", sep = "")
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  None\n")
    return(invisible(NULL))
  }
  for (row in seq_len(nrow(table))) {
    term <- if (identical(type, "ld")) {
      paste(table$item1[[row]], table$item2[[row]], sep = " -- ")
    } else {
      paste(table$item[[row]], table$exogenous[[row]], sep = " by ")
    }
    annotation <- model_term_annotation(table[row, , drop = FALSE])
    cat("  ", term, annotation, "\n", sep = "")
  }
  invisible(NULL)
}

print_fit_term_table <- function(table) {
  ld <- term_rows(table, "ld")
  dif <- term_rows(table, "dif")
  print_model_term_table("Local dependence terms", ld, "ld")
  print_model_term_table("DIF terms", dif, "dif")
  invisible(NULL)
}

term_rows <- function(table, type) {
  if (!is.data.frame(table) || !nrow(table) || !"type" %in% names(table)) {
    return(data.frame())
  }
  out <- table[table$type %in% type, , drop = FALSE]
  rownames(out) <- NULL
  out
}

model_term_annotation <- function(row) {
  source <- row$source[[1L]] %||% NA_character_
  status <- row$status[[1L]] %||% NA_character_
  if (identical(source, "user") && identical(status, "specified")) {
    return("")
  }
  details <- c()
  if (!is.na(status) && nzchar(status)) {
    details <- c(details, status)
  }
  if (!is.na(source) && nzchar(source)) {
    details <- c(details, paste("by", source))
  }
  if (!length(details)) {
    return("")
  }
  paste0(" (", paste(details, collapse = " "), ")")
}

summary_convergence_remarks <- function(tables) {
  rows <- list()
  for (name in names(tables)) {
    table <- tables[[name]]
    if (!is.data.frame(table) || !"converged" %in% names(table) || !nrow(table)) {
      next
    }
    nonconverged <- !is.na(table$converged) & !table$converged
    n_nonconverged <- sum(nonconverged)
    if (n_nonconverged == 0L) {
      next
    }
    reason_text <- ""
    if ("stop_reason" %in% names(table)) {
      reasons <- table$stop_reason[nonconverged & !is.na(table$stop_reason)]
      if (length(reasons)) {
        counts <- base::table(reasons)
        reason_text <- paste0(
          "; stop reasons: ",
          paste(paste0(names(counts), "=", as.integer(counts)), collapse = ", ")
        )
      }
    }
    rows[[length(rows) + 1L]] <- data.frame(
      table = name,
      n_nonconverged = as.integer(n_nonconverged),
      message = paste0(
        "Non-converged candidate fits: ",
        n_nonconverged,
        " candidate fit",
        if (n_nonconverged == 1L) "" else "s",
        " did not converge",
        " (",
        name,
        reason_text,
        "; see rows where converged == FALSE)."
      ),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) {
    return(data.frame(
      table = character(),
      n_nonconverged = integer(),
      message = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
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

reject_summary_which <- function(...) {
  dots <- list(...)
  if ("which" %in% names(dots)) {
    stop("This summary has one public view and does not accept `which`.", call. = FALSE)
  }
  invisible(NULL)
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

analysis_summary_header <- function(object) {
  data_name <- object$data_name %||% object$name %||% "<unnamed>"
  items <- object$items %||% character()
  exogenous <- object$exogenous %||% character()
  item_levels <- analysis_summary_level_counts(object$project$items)
  exogenous_levels <- analysis_summary_level_counts(object$project$backgrounds)
  id <- object$id %||% "none"
  score_groups <- analysis_summary_score_group_label(object)
  c(
    "",
    "Data",
    paste0("  Source: ", data_name),
    paste0("  Rows: ", nrow(object$project$raw_data)),
    paste0("  ID: ", id),
    "",
    "Variables",
    paste0("  Items: ", length(items), " (", summary_header_names(items), ")"),
    paste0("  Item levels: ", summary_header_names(item_levels, empty = "none")),
    paste0(
      "  Exogenous: ",
      length(exogenous),
      " (",
      summary_header_names(exogenous, empty = "none"),
      ")"
    ),
    paste0("  Exogenous levels: ", summary_header_names(exogenous_levels, empty = "none")),
    "",
    "Score groups",
    paste0("  Groups: ", score_groups)
  )
}

analysis_summary_score_group_label <- function(object) {
  distribution <- tryCatch(
    analysis_score_group_distribution(object),
    error = function(e) data.frame()
  )
  if (is.data.frame(distribution) && nrow(distribution)) {
    groups <- distribution[distribution$score != "Total", , drop = FALSE]
    if (nrow(groups)) {
      return(paste0(nrow(groups), " (", paste(groups$score, collapse = ", "), ")"))
    }
  }
  summary_header_names(object$score_groups %||% integer(), empty = "none")
}

analysis_score_group_distribution <- function(object) {
  values <- items_select_values(object$project)
  score_summary <- values$score_summary
  score_distribution <- values$score_distribution
  cuts <- as.integer(object$score_groups %||% integer())
  if (!length(cuts)) {
    return(data.frame())
  }

  n_complete <- as.integer(score_summary$n)
  rows <- list()
  from_score <- 0L
  cumulative_cases <- 0L
  for (group_index in seq_along(cuts)) {
    to_score <- cuts[[group_index]]
    if (group_index == length(cuts)) {
      to_score <- min(to_score, score_summary$observed_max)
    }
    if (from_score <= to_score) {
      in_group <- score_distribution$score >= from_score & score_distribution$score <= to_score
      cases <- as.integer(sum(score_distribution$count[in_group]))
      cumulative_cases <- cumulative_cases + cases
      rows[[length(rows) + 1L]] <- data.frame(
        score = score_range_label(from_score, to_score),
        count = cases,
        percent = 100 * cases / n_complete,
        cumulative = 100 * cumulative_cases / n_complete,
        stringsAsFactors = FALSE
      )
    }
    from_score <- cuts[[group_index]] + 1L
  }

  rows[[length(rows) + 1L]] <- data.frame(
    score = "Total",
    count = n_complete,
    percent = 100,
    cumulative = 100,
    stringsAsFactors = FALSE
  )
  out <- do.call(rbind, rows)
  attr(out, "observed_min") <- as.integer(score_summary$observed_min)
  attr(out, "observed_max") <- as.integer(score_summary$observed_max)
  attr(out, "missing_item_score_rows") <- as.integer(score_summary$missing)
  out
}

score_range_label <- function(from_score, to_score) {
  if (identical(as.integer(from_score), as.integer(to_score))) {
    return(as.character(as.integer(from_score)))
  }
  paste0(as.integer(from_score), "-", as.integer(to_score))
}

print_analysis_score_group_distribution <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    return(invisible(NULL))
  }
  observed_min <- attr(table, "observed_min", exact = TRUE)
  observed_max <- attr(table, "observed_max", exact = TRUE)
  missing <- attr(table, "missing_item_score_rows", exact = TRUE)

  cat("\nScore group distribution\n")
  cat("Observed score range: ", observed_min, "-", observed_max, "\n", sep = "")
  cat("Score-group cases: ", sum(table$count[table$score != "Total"]), "\n", sep = "")
  cat("Missing item-score rows: ", missing, "\n\n", sep = "")

  display <- table
  display$percent <- sprintf("%.1f", display$percent)
  display$cumulative <- sprintf("%.1f", display$cumulative)
  names(display) <- c("Score", "Count", "Percent", "Cumulative")
  print(display, row.names = FALSE)
  invisible(NULL)
}

analysis_summary_level_counts <- function(variables) {
  if (!is.data.frame(variables) || !nrow(variables)) {
    return(character())
  }
  paste0(variables$name, "=", variables$raw_max)
}

summary_header_names <- function(x, empty = "none") {
  if (!length(x)) {
    return(empty)
  }
  paste(as.character(x), collapse = ", ")
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
  paste("gRm:", public_model_type_label(public_model_type(object)), noun)
}

public_model_type_label <- function(model_type) {
  if (identical(model_type, "rasch")) {
    "Rasch model"
  } else {
    "Graphical Log-Linear Rasch Model"
  }
}

public_screen_terms <- function(object) {
  if (!is.null(object$terms)) {
    return(object$terms)
  }
  model_terms(object)
}

screen_summary_header <- function(object, tables) {
  exact_state <- object$exact_state %||% list()
  bh <- object$values$bh %||% list()
  header <- c(
    "",
    "Screening",
    paste0("  Inference: ", object$inference %||% "asymptotic")
  )
  if (isTRUE(exact_state$exact)) {
    header <- c(
      header,
      paste0("  Simulations: ", summary_scalar(object$nsim %||% exact_state$nsim %||% NA_integer_))
    )
    if (isTRUE(exact_state$sequential)) {
      header <- c(
        header,
        paste0("  Sequential limit: ", summary_scalar(exact_state$seq_limit %||% NA_integer_))
      )
    }
    header <- c(
      header,
      paste0("  Seed: ", summary_scalar(object$seed %||% exact_state$seed %||% NA_integer_))
    )
  }
  c(
    header,
    paste0("  Tested relations: ", summary_scalar(bh$n_tests %||% NA_integer_)),
    paste0("  Tested local-dependence pairs: ", nrow(tables$local_dependence)),
    paste0("  Tested DIF relations: ", nrow(tables$dif)),
    paste0("  Tested score effects: ", nrow(tables$score_effects))
  )
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

public_screen_summary_tables <- function(object) {
  values <- object$values %||% list()
  terms <- public_screen_terms(object)
  model_terms <- rbind_fill(terms$ld, terms$dif)
  selected_ld <- public_selected_rows(terms$ld)
  selected_dif <- public_selected_rows(terms$dif)
  list(
    local_dependence = public_screen_local_dependence_tests(values),
    dif = public_screen_dif_tests(values),
    score_effects = public_screen_score_effect_tests(values),
    selected = rbind_fill(selected_ld, selected_dif),
    selected_ld = selected_ld,
    selected_dif = selected_dif,
    model_terms = model_terms,
    bh = public_screen_bh_table(object)
  )
}

public_screen_local_dependence_tests <- function(values) {
  p <- values$partial$item_p %||% matrix(numeric(), nrow = 0L, ncol = 0L)
  wpg <- values$partial$weighted_gamma %||% matrix(numeric(), nrow = 0L, ncol = 0L)
  items <- values$items %||% data.frame()
  n_items <- nrow(items)
  if (!is.matrix(p) || n_items < 2L) {
    return(data.frame(
      `Item 1` = character(),
      `Item 2` = character(),
      WPG = numeric(),
      `Pr(>|WPG|)` = numeric(),
      ` ` = character(),
      check.names = FALSE
    ))
  }
  selected <- values$model$local_dependence$matrix %||% matrix(FALSE, nrow = n_items, ncol = n_items)
  rows <- which(upper.tri(matrix(FALSE, nrow = n_items, ncol = n_items)), arr.ind = TRUE)
  item_names <- screen_item_names(items, n_items)
  out <- data.frame(
    `Item 1` = item_names[rows[, "row"]],
    `Item 2` = item_names[rows[, "col"]],
    WPG = screen_matrix_value(wpg, rows),
    `Pr(>|WPG|)` = screen_matrix_value(p, rows),
    ` ` = ifelse(screen_matrix_value(selected, rows, default = FALSE) %in% TRUE, "*", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out[!is.na(out[["Pr(>|WPG|)"]]) & out[["Pr(>|WPG|)"]] > 1, "Pr(>|WPG|)"] <- NA_real_
  out
}

public_screen_dif_tests <- function(values) {
  p <- values$partial$exo_p %||% matrix(numeric(), nrow = 0L, ncol = 0L)
  stat <- values$partial$exo_stat %||% matrix(numeric(), nrow = 0L, ncol = 0L)
  kinds <- values$partial$exo_kind %||% character()
  items <- values$items %||% data.frame()
  backgrounds <- values$backgrounds %||% data.frame()
  n_items <- nrow(items)
  n_exo <- nrow(backgrounds)
  if (!is.matrix(p) || n_items == 0L || n_exo == 0L) {
    return(data.frame(
      Item = character(),
      Exogenous = character(),
      Chisq = numeric(),
      Df = integer(),
      `Pr(>Chisq)` = numeric(),
      Gamma = numeric(),
      `Pr(>|Gamma|)` = numeric(),
      ` ` = character(),
      check.names = FALSE
    ))
  }
  rows <- expand.grid(item = seq_len(n_items), exogenous = seq_len(n_exo))
  item_names <- screen_item_names(items, n_items)
  exo_names <- screen_item_names(backgrounds, n_exo)
  selected <- values$model$item_bias %||% matrix(FALSE, nrow = n_items, ncol = n_exo)
  kind <- rep_len(as.character(kinds), n_exo)
  row_kind <- kind[rows$exogenous]
  use_gamma <- row_kind == "Gamma"
  row_index <- cbind(rows$item, rows$exogenous)
  statistic <- stat[row_index]
  p_value <- p[row_index]
  data.frame(
    Item = item_names[rows$item],
    Exogenous = exo_names[rows$exogenous],
    Chisq = ifelse(use_gamma, NA_real_, statistic),
    Df = NA_integer_,
    `Pr(>Chisq)` = ifelse(use_gamma, NA_real_, p_value),
    Gamma = ifelse(use_gamma, statistic, NA_real_),
    `Pr(>|Gamma|)` = ifelse(use_gamma, p_value, NA_real_),
    ` ` = ifelse(screen_matrix_value(selected, row_index, default = FALSE) %in% TRUE, "*", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

public_screen_score_effect_tests <- function(values) {
  rows <- values$model$score_effects$rows %||% data.frame()
  if (!is.data.frame(rows) || !nrow(rows)) {
    return(data.frame(
      Exogenous = character(),
      Chisq = numeric(),
      Df = integer(),
      `Pr(>Chisq)` = numeric(),
      Gamma = numeric(),
      `Pr(>|Gamma|)` = numeric(),
      ` ` = character(),
      check.names = FALSE
    ))
  }
  out <- data.frame(
    Exogenous = diagnostic_column(rows, "name", NA_character_),
    Chisq = diagnostic_column(rows, "chi_square", NA_real_),
    Df = diagnostic_column(rows, "df", NA_integer_),
    `Pr(>Chisq)` = diagnostic_column(rows, "p_chi", NA_real_),
    Gamma = diagnostic_column(rows, "gamma", NA_real_),
    `Pr(>|Gamma|)` = diagnostic_column(rows, "p_gamma", NA_real_),
    ` ` = ifelse(diagnostic_column(rows, "selected", FALSE) %in% TRUE, "*", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  exact_nsim <- diagnostic_column(rows, "exact_nsim", NA_integer_)
  if (any(!is.na(exact_nsim) & exact_nsim > 0L)) {
    out$Simulations <- exact_nsim
  }
  out
}

screen_item_names <- function(table, n) {
  names <- if (is.data.frame(table) && "name" %in% names(table)) {
    as.character(table$name)
  } else {
    character()
  }
  if (length(names) < n) {
    names <- c(names, as.character(seq.int(length(names) + 1L, n)))
  }
  names
}

screen_matrix_value <- function(matrix, index, default = NA_real_) {
  if (!is.matrix(matrix) || !nrow(index)) {
    return(rep(default, nrow(index)))
  }
  valid <- index[, 1L] >= 1L & index[, 1L] <= nrow(matrix) &
    index[, 2L] >= 1L & index[, 2L] <= ncol(matrix)
  out <- rep(default, nrow(index))
  out[valid] <- matrix[index[valid, , drop = FALSE]]
  out
}

public_local_dependence_tests <- function(tests, bh_critical_p = NA_real_) {
  if (!is.data.frame(tests)) {
    tests <- data.frame()
  }
  converged <- diagnostic_column(tests, "converged", NA)
  p_value <- diagnostic_column(tests, "p_value", NA_real_)
  threshold <- bh_critical_p %||% NA_real_
  threshold <- threshold[[1L]]
  data.frame(
    `Item 1` = diagnostic_column(tests, "item1_name", NA_character_),
    `Item 2` = diagnostic_column(tests, "item2_name", NA_character_),
    Chisq = diagnostic_column(tests, "chi_square", NA_real_),
    Df = diagnostic_column(tests, "degrees_of_freedom", NA_integer_),
    `Pr(>Chisq)` = p_value,
    WPG = diagnostic_column(tests, "wpg_gamma", NA_real_),
    Converged = yes_no_display(converged),
    delta = diagnostic_column(tests, "delta", NA_real_),
    ` ` = ifelse(!is.na(p_value) & !is.na(threshold) & p_value <= threshold, "*", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

public_dif_tests <- function(tests, bh_critical_p = NA_real_) {
  if (!is.data.frame(tests)) {
    tests <- data.frame()
  }
  converged <- diagnostic_column(tests, "converged", NA)
  stable <- diagnostic_column(tests, "output_stable", FALSE)
  p_value <- diagnostic_column(tests, "p_value", NA_real_)
  threshold <- bh_critical_p %||% NA_real_
  threshold <- threshold[[1L]]
  data.frame(
    Item = diagnostic_column(tests, "item_name", NA_character_),
    Exogenous = diagnostic_column(tests, "background_name", NA_character_),
    Chisq = diagnostic_column(tests, "chi_square", NA_real_),
    Df = diagnostic_column(tests, "degrees_of_freedom", NA_integer_),
    `Pr(>Chisq)` = p_value,
    Gamma = diagnostic_column(tests, "gamma", NA_real_),
    Converged = yes_no_display(converged),
    Stable = yes_no_display(stable),
    delta = diagnostic_column(tests, "delta", NA_real_),
    ` ` = ifelse(!is.na(p_value) & !is.na(threshold) & p_value <= threshold, "*", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

public_global_homogeneity_test <- function(summary) {
  summary <- list_to_one_row(summary %||% list())
  data.frame(
    `Score groups` = diagnostic_column(summary, "n_groups", NA_integer_),
    Parameters = diagnostic_column(summary, "n_parameters", NA_integer_),
    `LogLik full` = diagnostic_column(summary, "full_log_likelihood", NA_real_),
    `LogLik groups` = diagnostic_column(summary, "subgroup_log_likelihood_sum", NA_real_),
    CLR = diagnostic_column(summary, "clr", NA_real_),
    Df = diagnostic_column(summary, "df", NA_integer_),
    `Pr(>CLR)` = diagnostic_column(summary, "p_value", NA_real_),
    check.names = FALSE
  )
}

public_global_homogeneity_score_groups <- function(groups) {
  if (!is.data.frame(groups)) {
    groups <- data.frame()
  }
  data.frame(
    `Score group` = global_homogeneity_score_group_labels(groups),
    Cases = diagnostic_column(groups, "n", NA_integer_),
    LogLik = diagnostic_column(groups, "log_likelihood", NA_real_),
    Converged = yes_no_display(diagnostic_column(groups, "converged", NA)),
    delta = diagnostic_column(groups, "delta", NA_real_),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

public_global_homogeneity_item_means <- function(items, groups) {
  if (!is.data.frame(items)) {
    items <- data.frame()
  }
  data.frame(
    `Score group` = global_homogeneity_item_score_group_labels(items, groups),
    Item = diagnostic_column(items, "item_name", NA_character_),
    Cases = diagnostic_column(items, "n", NA_integer_),
    `Observed mean` = diagnostic_column(items, "observed_mean", NA_real_),
    `Expected mean` = diagnostic_column(items, "expected_mean", NA_real_),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

public_global_homogeneity_uniform_ld <- function(table, values) {
  if (!is.data.frame(table) || !nrow(table)) {
    return(data.frame())
  }
  out <- data.frame(
    `Item 1` = global_homogeneity_label_names(
      diagnostic_column(table, "item1_label", NA_character_),
      global_homogeneity_context_variables(values, "items")
    ),
    `Item 2` = global_homogeneity_label_names(
      diagnostic_column(table, "item2_label", NA_character_),
      global_homogeneity_context_variables(values, "items")
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- cbind(
    out,
    global_homogeneity_uniform_gamma_columns(table, values),
    data.frame(
      Chisq = diagnostic_column(table, "chi_square", NA_real_),
      Df = diagnostic_column(table, "df", NA_integer_),
      `Pr(>Chisq)` = diagnostic_column(table, "p_value", NA_real_),
      check.names = FALSE
    )
  )
  rownames(out) <- NULL
  out
}

public_global_homogeneity_uniform_dif <- function(table, values) {
  if (!is.data.frame(table) || !nrow(table)) {
    return(data.frame())
  }
  out <- data.frame(
    Item = global_homogeneity_label_names(
      diagnostic_column(table, "item_label", NA_character_),
      global_homogeneity_context_variables(values, "items")
    ),
    Exogenous = global_homogeneity_label_names(
      diagnostic_column(table, "background_label", NA_character_),
      global_homogeneity_context_variables(values, "backgrounds")
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- cbind(
    out,
    global_homogeneity_uniform_gamma_columns(table, values),
    data.frame(
      Chisq = diagnostic_column(table, "chi_square", NA_real_),
      Df = diagnostic_column(table, "df", NA_integer_),
      `Pr(>Chisq)` = diagnostic_column(table, "p_value", NA_real_),
      check.names = FALSE
    )
  )
  rownames(out) <- NULL
  out
}

global_homogeneity_default_summary_sections <- function(values) {
  sections <- c("test", "score_groups", "item_means")
  if (global_homogeneity_table_nrow(values$uniform_ld %||% data.frame()) > 0L) {
    sections <- c(sections, "uniform_ld")
  }
  if (global_homogeneity_table_nrow(values$uniform_dif %||% data.frame()) > 0L) {
    sections <- c(sections, "uniform_dif")
  }
  sections
}

global_homogeneity_table_nrow <- function(table) {
  if (is.data.frame(table)) {
    nrow(table)
  } else {
    0L
  }
}

global_homogeneity_uniform_gamma_columns <- function(table, values) {
  labels <- global_homogeneity_uniform_score_group_labels(table, values$score_groups %||% data.frame())
  n_groups <- length(labels)
  if (n_groups == 0L) {
    return(data.frame(row.names = seq_len(nrow(table))))
  }
  observed <- global_homogeneity_uniform_gamma_matrix(table, "observed_gamma", n_groups)
  expected <- global_homogeneity_uniform_gamma_matrix(table, "expected_gamma", n_groups)
  columns <- vector("list", 2L * n_groups)
  column_names <- character(2L * n_groups)
  for (group_index in seq_len(n_groups)) {
    observed_column <- (group_index - 1L) * 2L + 1L
    expected_column <- observed_column + 1L
    columns[[observed_column]] <- observed[, group_index]
    columns[[expected_column]] <- expected[, group_index]
    column_names[[observed_column]] <- paste("Obs gamma", labels[[group_index]])
    column_names[[expected_column]] <- paste("Exp gamma", labels[[group_index]])
  }
  out <- as.data.frame(columns, stringsAsFactors = FALSE, check.names = FALSE)
  names(out) <- column_names
  out
}

global_homogeneity_uniform_score_group_labels <- function(table, groups) {
  labels <- global_homogeneity_score_group_labels(groups)
  n_groups <- max(
    length(labels),
    global_homogeneity_uniform_n_groups(table),
    0L
  )
  if (n_groups == 0L) {
    return(character())
  }
  if (length(labels) < n_groups) {
    labels <- c(labels, as.character(seq.int(length(labels) + 1L, n_groups)))
  }
  labels[seq_len(n_groups)]
}

global_homogeneity_uniform_n_groups <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    return(0L)
  }
  lengths <- integer()
  for (column in c("observed_gamma", "expected_gamma")) {
    if (column %in% names(table)) {
      lengths <- c(lengths, vapply(table[[column]], length, integer(1L)))
    }
  }
  if (length(lengths)) {
    max(lengths, na.rm = TRUE)
  } else {
    0L
  }
}

global_homogeneity_uniform_gamma_matrix <- function(table, column, n_groups) {
  out <- matrix(NA_real_, nrow = nrow(table), ncol = n_groups)
  if (!column %in% names(table) || n_groups == 0L) {
    return(out)
  }
  for (row_index in seq_len(nrow(table))) {
    values <- as.numeric(table[[column]][[row_index]])
    n <- min(length(values), n_groups)
    if (n > 0L) {
      out[row_index, seq_len(n)] <- values[seq_len(n)]
    }
  }
  out
}

global_homogeneity_context_variables <- function(values, kind = c("items", "backgrounds")) {
  kind <- match.arg(kind)
  context <- values$fit$context %||% list()
  variables <- context[[kind]] %||% NULL
  if (is.data.frame(variables)) {
    return(variables)
  }
  model <- values$bundle$model %||% list()
  variables <- model[[kind]] %||% NULL
  if (is.data.frame(variables)) {
    return(variables)
  }
  data.frame()
}

global_homogeneity_label_names <- function(labels, variables) {
  labels <- as.character(labels)
  if (!is.data.frame(variables) || !all(c("label_code", "name") %in% names(variables))) {
    return(labels)
  }
  index <- match(labels, as.character(variables$label_code))
  out <- as.character(variables$name[index])
  missing <- is.na(out)
  out[missing] <- labels[missing]
  out
}

global_homogeneity_score_group_labels <- function(groups) {
  if (!is.data.frame(groups) || !nrow(groups)) {
    return(character())
  }
  if (all(c("from_score", "to_score") %in% names(groups))) {
    return(mapply(score_range_label, groups$from_score, groups$to_score, USE.NAMES = FALSE))
  }
  if ("label" %in% names(groups)) {
    return(gsub("[[:space:]]+-[[:space:]]+", "-", groups$label))
  }
  as.character(diagnostic_column(groups, "group", NA_integer_))
}

global_homogeneity_item_score_group_labels <- function(items, groups) {
  if (!is.data.frame(items) || !nrow(items)) {
    return(character())
  }
  item_groups <- diagnostic_column(items, "group", NA_integer_)
  if (!is.data.frame(groups) || !nrow(groups) || !"group" %in% names(groups)) {
    return(as.character(item_groups))
  }
  group_labels <- global_homogeneity_score_group_labels(groups)
  names(group_labels) <- as.character(groups$group)
  out <- unname(group_labels[as.character(item_groups)])
  missing <- is.na(out)
  out[missing] <- as.character(item_groups[missing])
  out
}

diagnostic_column <- function(tests, column, default) {
  if (column %in% names(tests)) {
    return(tests[[column]])
  }
  rep(default, nrow(tests))
}

yes_no_display <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    ifelse(x, "yes", "no")
  )
}

public_bh_marker_note <- function(bh, digits = max(3L, getOption("digits") - 3L)) {
  if (!is.data.frame(bh) || !nrow(bh)) {
    return(character())
  }
  fdr <- bh$fdr[[1L]]
  threshold <- bh$p_value[[1L]]
  if (is.na(fdr) || is.na(threshold)) {
    return(character())
  }
  dig_tst <- max(1L, min(5L, digits - 1L))
  paste0(
    "*: p <= Benjamini-Hochberg threshold for FDR = ",
    fdr,
    " (threshold = ",
    format.pval(threshold, digits = dig_tst, eps = .Machine$double.eps),
    ")"
  )
}

public_screen_bh_marker_note <- function(bh, digits = max(3L, getOption("digits") - 3L)) {
  if (!is.data.frame(bh) || !nrow(bh)) {
    return(character())
  }
  fdr05 <- which(as.character(bh$fdr) == "0.05")
  if (!length(fdr05)) {
    fdr05 <- 1L
  }
  fdr <- bh$fdr[[fdr05[[1L]]]]
  threshold <- bh$p_value[[fdr05[[1L]]]]
  if (is.na(fdr) || is.na(threshold)) {
    return(character())
  }
  dig_tst <- max(1L, min(5L, digits - 1L))
  paste0(
    "*: selected by the SCREEN J source decision path at the 5% level. ",
    "Local-dependence and DIF candidate evidence uses the global ",
    "Benjamini-Hochberg threshold for FDR = ",
    fdr,
    " (threshold = ",
    format.pval(threshold, digits = dig_tst, eps = .Machine$double.eps),
    "); stricter global cutoffs are available in attr(x, \"bh\")."
  )
}

public_value_tables <- function(values) {
  if (inherits(values, "gRm_item_parameters_values") || inherits(values, "gRm_gllrm_values")) {
    return(item_parameter_result_tables(values))
  }

  if (inherits(values, "gRm_item_fits_values")) {
    extended <- values$extended %||% list()
    return(list(
      statistics = normalize_summary_table(values$items %||% data.frame()),
      compact = normalize_summary_table(values$side_file %||% data.frame()),
      bh_thresholds = result_named_numeric_table(values$bh_limits %||% numeric(), "threshold", "p_value"),
      score_n = normalize_summary_table(extended$score_n %||% data.frame()),
      score_level_fit = normalize_summary_table(extended$scores %||% data.frame()),
      item_fit_summaries = normalize_summary_table(extended$summaries %||% data.frame())
    ))
  }

  if (inherits(values, "gRm_local_independence_values")) {
    return(list(
      tests = normalize_summary_table(values$tests %||% data.frame()),
      selected = normalize_summary_table(public_selected_by_bh(values$tests, values$bh_critical_p)),
      bh_thresholds = data.frame(
        result = "missing_ld",
        fdr = "0.05",
        p_value = values$bh_critical_p %||% NA_real_,
        stringsAsFactors = FALSE
      )
    ))
  }

  if (inherits(values, "gRm_dif_tests_values")) {
    tables <- list(
      tests = normalize_summary_table(values$tests %||% data.frame()),
      selected = normalize_summary_table(public_selected_by_bh(values$tests, values$bh_critical_p)),
      bh_thresholds = data.frame(
        result = "missing_dif",
        fdr = "0.05",
        p_value = values$bh_critical_p %||% NA_real_,
        stringsAsFactors = FALSE
      )
    )
    if (!is.null(values$included_tests)) {
      tables$included_tests <- normalize_summary_table(values$included_tests)
    }
    return(tables)
  }

  if (inherits(values, "gRm_exo_select_values")) {
    return(list(
      score_effect_tests = normalize_summary_table(values$screen %||% data.frame()),
      score_effect_selected = normalize_summary_table(values$selected %||% data.frame()),
      bh_thresholds = list_to_one_row(values$bh %||% list())
    ))
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
