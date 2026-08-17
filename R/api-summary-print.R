#' @export
print.summary.gRm <- function(x, ...) {
  reject_public_dots(...)
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
  if (inherits(x, "summary.gRm_screen")) {
    cat("\nSelected model terms\n")
    print_fit_term_table(x$selected)
  }
  invisible(x)
}

#' Internal print summary table note helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param note Internal `note` value used by this helper.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
print_summary_table_note <- function(note) {
  if (!length(note) || all(is.na(note)) || !nzchar(note[[1L]])) {
    return(invisible(NULL))
  }
  cat("---\n")
  cat(note[[1L]], "\n", sep = "")
  invisible(NULL)
}

#' Internal print summary table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
print_summary_table <- function(table) {
  if (is.data.frame(table) && nrow(table)) {
    print(public_format_table(table), row.names = FALSE)
  } else {
    cat("  <none>\n")
  }
  invisible(NULL)
}

#' Internal print diagnostic tests table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @param gamma_columns Internal `gamma_columns` value used by this helper.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
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

#' Internal print screen tests table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
print_screen_tests_table <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  <none>\n")
    return(invisible(NULL))
  }
  display <- table
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  statistic_columns <- c(
    "Chisq", "Gamma", "WPG", "Gamma 1->2", "Gamma 2->1", "Gamma sum"
  )
  p_value_columns <- c(
    "Pr(>Chisq)", "Pr(>|Gamma|)",
    "Pr(>|Gamma 1->2|)", "Pr(>|Gamma 2->1|)"
  )
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

#' Internal print global homogeneity test table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
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

#' Internal print global homogeneity score group table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
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

#' Internal print global homogeneity item means table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
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

#' Internal print score effects tests table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
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

#' Internal print global homogeneity uniform table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
print_global_homogeneity_uniform_table <- function(table) {
  gamma_columns <- grep("^(Obs|Exp) gamma ", names(table), value = TRUE)
  print_diagnostic_tests_table(table, gamma_columns = gamma_columns)
}

#' Internal print item fit tests table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
print_item_fit_tests_table <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    cat("  <none>\n")
    return(invisible(NULL))
  }
  display <- table
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  statistic_columns <- c("Outfit", "Outfit SE", "Infit", "Infit SE")
  p_value_markers <- c(
    `Pr(>Outfit)` = "Outfit FDR",
    `Pr(>Infit)` = "Infit FDR",
    `Pr(>Gamma)` = "Gamma FDR"
  )
  gamma_columns <- c("Observed gamma", "Expected gamma", "Gamma SE")
  for (column in statistic_columns) {
    if (column %in% names(display)) {
      display[[column]] <- format(signif(display[[column]], dig_tst), digits = dig_tst)
    }
  }
  for (column in names(p_value_markers)) {
    if (column %in% names(display)) {
      formatted <- format.pval(
        display[[column]],
        digits = dig_tst,
        eps = .Machine$double.eps
      )
      marker_column <- unname(p_value_markers[[column]])
      marker <- if (marker_column %in% names(display)) {
        as.character(display[[marker_column]])
      } else {
        rep("", nrow(display))
      }
      marker[is.na(marker)] <- ""
      # DIGRAM writes every FDR marker into a three-character suffix field
      # (`*  `, `** `, `***`, or three blanks). Reserving that field keeps the
      # numeric p-values aligned independently of their marker grade.
      display[[column]] <- paste0(formatted, sprintf("%-3s", marker))
    }
  }
  display[intersect(unname(p_value_markers), names(display))] <- NULL
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

#' Internal print item fit items table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
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

#' Internal diagnostic print max helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The internal `diagnostic_print_max()` computation result.
#' @keywords internal
#' @noRd
diagnostic_print_max <- function(table) {
  cells <- as.double(nrow(table)) * as.double(ncol(table)) + 1
  as.integer(min(cells, .Machine$integer.max))
}

#' Internal print model summary table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
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

#' Internal print fit summary table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
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
  cat(
    "  Negative log likelihood (DIGRAM): ",
    summary_scalar(row$negative_log_likelihood %||% row$log_likelihood[[1L]]),
    "\n",
    sep = ""
  )
  cat("  Parameters: ", summary_scalar(row$n_parameters[[1L]]), "\n", sep = "")
  cat("  Likelihood rows: ", summary_scalar(row$likelihood_n[[1L]]), "\n", sep = "")
  invisible(NULL)
}

#' Internal summary scalar helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @return The internal `summary_scalar()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal summary p value helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @return The internal `summary_p_value()` computation result.
#' @keywords internal
#' @noRd
summary_p_value <- function(x) {
  if (length(x) == 0L || is.na(x)) {
    return("NA")
  }
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  format.pval(x, digits = dig_tst, eps = .Machine$double.eps)
}

#' Internal print fit parameter table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
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

#' Internal print fit threshold table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
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

#' Internal fit threshold display names helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param names Internal `names` value used by this helper.
#' @return The internal `fit_threshold_display_names()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal fit parameter display names helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param names Internal `names` value used by this helper.
#' @return The internal `fit_parameter_display_names()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal print model term table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param title Internal `title` value used by this helper.
#' @param table Numeric contingency or result table.
#' @param type Internal `type` value used by this helper.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
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

#' Internal print fit term table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
print_fit_term_table <- function(table) {
  ld <- term_rows(table, "ld")
  dif <- term_rows(table, "dif")
  print_model_term_table("Local dependence terms", ld, "ld")
  print_model_term_table("DIF terms", dif, "dif")
  invisible(NULL)
}

#' Internal term rows helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @param type Internal `type` value used by this helper.
#' @return The internal `term_rows()` computation result.
#' @keywords internal
#' @noRd
term_rows <- function(table, type) {
  if (!is.data.frame(table) || !nrow(table) || !"type" %in% names(table)) {
    return(data.frame())
  }
  out <- table[table$type %in% type, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Internal model term annotation helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param row Internal `row` value used by this helper.
#' @return The internal `model_term_annotation()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal summary convergence remarks helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param tables Internal `tables` value used by this helper.
#' @return The internal `summary_convergence_remarks()` computation result.
#' @keywords internal
#' @noRd
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
