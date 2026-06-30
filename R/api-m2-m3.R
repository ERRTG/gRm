#' M2 fit diagnostic
#'
#' Compute DIGRAM's source-backed M2 diagnostic for a fitted gRm model.
#'
#' @param fit A fitted `gRm_fit` object returned by [fit()].
#' @param items Optional item selection for this diagnostic only. `NULL`
#'   selects all items in source item order. A non-null value may be R-facing
#'   item names or one-based item indices. The selection is normalized to
#'   source order because DIGRAM stores item selection as a set and then loops
#'   over source item indices.
#' @param score_cuts Optional integer-like upper total-score cut values for
#'   this diagnostic result. `NULL` uses the score cuts stored on the
#'   `gRm_analysis`. A non-null value supplies the temporary DIGRAM `scorecuts`
#'   state used for score-group rows in this result only; it does not mutate the
#'   analysis object and does not refit the model.
#' @param ... Reserved for future extensions and must be empty.
#' @return A `gRm_m2` result object. Programmatic tables are available in
#'   `result$values`.
#' @details
#' `m2()` is a diagnostic for the current fitted model. It does not search for,
#' select, add, or refit local-dependence or DIF terms.
#'
#' In the DIGRAM source, M2 is not a separately registered command. It is the
#' two-way-margin part of DIGRAM's CM3 command. The source `Prepare_CM3tests`
#' routine prepares item-item margins, item-exogenous margins, and
#' item-score-group margins for the selected items. Item-item margins that are
#' already included as local-dependence terms are skipped. Item-exogenous rows
#' follow the source `ItemBias(.i2,i1.)` check in `Prepare_CM3tests`: because
#' the stored DIF matrix is otherwise item-first and exogenous-second, this
#' branch uses the DIF matrix in transposed index order. Active DIF rows can
#' therefore remain in the M2 table as zero rows, while transposed source
#' artifact rows can be suppressed. Item-score-group rows are included for
#' selected items.
#'
#' The M2 aggregate is the sum of the prepared two-way Pearson chi-square
#' contributions. Degrees of freedom and p-values follow the DIGRAM source,
#' using the package's `source_pfchi()` implementation for chi-square tail
#' probabilities. Public output uses R variable names. DIGRAM's compact
#' alphabetic labels are retained internally for source tracing and oracle
#' validation.
#'
#' Score-group rows use the resolved diagnostic score cuts. Supplying
#' `score_cuts` changes only this diagnostic result. Bootstrap calibration is
#' not implemented in the initial deterministic/asymptotic API.
#' @export
#' @examples
#' \donttest{
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(1, 2, 1, 2, 1, 2, 2, 1),
#'   I2 = c(2, 1, 2, 1, 2, 1, 2, 1),
#'   I3 = c(1, 1, 2, 1, 2, 2, 1, 2)
#' )
#' analysis <- gRm(data, items = c("I1", "I2", "I3"), id = "ID")
#' fit0 <- fit(gllrm(analysis))
#' m2_result <- m2(fit0)
#' summary(m2_result)
#' }
m2 <- function(fit, items = NULL, score_cuts = NULL, ...) {
  reject_public_dots(...)
  fit <- as_public_gRm_fit(fit)
  values <- m2_values(fit, items = items, score_cuts = score_cuts)
  new_m2_m3_result(fit, values, "gRm_m2", "m2", match.call())
}

#' M3 fit diagnostic
#'
#' Compute DIGRAM's source-backed M3 diagnostic for a fitted gRm model.
#'
#' @inheritParams m2
#' @return A `gRm_m3` result object. Programmatic tables are available in
#'   `result$values`.
#' @details
#' `m3()` is a diagnostic for the current fitted model. It does not search for,
#' select, add, or refit local-dependence or DIF terms.
#'
#' `m3()` implements the deterministic/asymptotic output of DIGRAM's CM3
#' command for the selected items. DIGRAM's CM3 output includes M2 first; the
#' M3 aggregate is therefore not a three-way-only statistic. It includes the
#' complete M2 two-way aggregate and then adds the prepared three-way margin
#' contributions.
#'
#' The three-way rows are prepared after the M2 rows in Pascal source order:
#' item-item-item rows, item-item-exogenous and item-item-score-group rows, and
#' then item-exogenous-exogenous and item-exogenous-score-group rows. Included
#' two-way local-dependence or DIF terms do not remove related three-way rows;
#' this is source behavior from `Prepare_CM3tests`.
#'
#' Item selection, public variable names, score-group handling, and bootstrap
#' scope are the same as for [m2()]. Public output uses R variable names, while
#' DIGRAM alphabetic labels remain internal for source/oracle matching.
#' @export
#' @examples
#' \donttest{
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(1, 2, 1, 2, 1, 2, 2, 1),
#'   I2 = c(2, 1, 2, 1, 2, 1, 2, 1),
#'   I3 = c(1, 1, 2, 1, 2, 2, 1, 2)
#' )
#' analysis <- gRm(data, items = c("I1", "I2", "I3"), id = "ID")
#' fit0 <- fit(gllrm(analysis))
#' m3_result <- m3(fit0)
#' summary(m3_result)
#' }
m3 <- function(fit, items = NULL, score_cuts = NULL, ...) {
  reject_public_dots(...)
  fit <- as_public_gRm_fit(fit)
  values <- m3_values(fit, items = items, score_cuts = score_cuts)
  new_m2_m3_result(fit, values, "gRm_m3", "m3", match.call())
}

new_m2_m3_result <- function(fit, values, class, result, call) {
  analysis <- fit$analysis %||% fit$spec$analysis
  selected <- values$selected_items %||% data.frame()
  metadata <- list(
    items = selected$item_name %||% character(),
    item_indices = selected$item_index %||% integer(),
    score_cuts = values$score_cuts %||% integer()
  )
  new_gRm_result(
    class = class,
    analysis = analysis,
    fit = fit,
    values = values,
    result = result,
    metadata = metadata,
    call = call
  )
}

#' @export
print.gRm_m2 <- function(x, ...) {
  print_m2_m3_result(x, include_m3 = FALSE)
}

#' @export
print.gRm_m3 <- function(x, ...) {
  print_m2_m3_result(x, include_m3 = TRUE)
}

print_m2_m3_result <- function(x, include_m3) {
  values <- x$values %||% list()
  title <- if (isTRUE(include_m3)) "gRm: M3 fit diagnostic" else "gRm: M2 fit diagnostic"
  cat(title, "\n\n", sep = "")
  cat("  Selection: ", m2_m3_selection_label(values), "\n", sep = "")
  cat("  Selected items: ", summary_scalar(nrow(values$selected_items %||% data.frame())), "\n", sep = "")
  cat("  Score cuts: ", m2_m3_score_cuts_label(values$score_cuts %||% integer()), "\n", sep = "")
  cat("  Two-way margins: ", summary_scalar(values$n_two_way_margins %||% 0L), "\n", sep = "")
  if (isTRUE(include_m3)) {
    cat("  Three-way margins: ", summary_scalar(values$n_three_way_margins %||% 0L), "\n", sep = "")
  }

  m2_row <- if (isTRUE(include_m3)) values$m2 else values$aggregate
  print_m2_m3_stat_block("M2", m2_row)
  if (isTRUE(include_m3)) {
    print_m2_m3_stat_block("M3", values$m3)
  }
  print_m2_m3_stat_block("Item-trait", values$item_trait)
  cat("  Invariance terms: ", summary_scalar(nrow(values$invariance %||% data.frame())), "\n\n", sep = "")
  if (isTRUE(include_m3)) {
    cat("Use summary(x) to show two-way margins, three-way margins, and decompositions.\n")
  } else {
    cat("Use summary(x) to show two-way margins and invariance by exogenous variable.\n")
  }
  invisible(x)
}

print_m2_m3_stat_block <- function(label, row) {
  row <- row %||% data.frame()
  if (!is.data.frame(row) || nrow(row) == 0L) {
    chi <- df <- p <- NA
  } else {
    chi <- row$chi_square[[1L]]
    df <- row$degrees_of_freedom[[1L]]
    p <- row$p_value[[1L]]
  }
  cat("\n")
  cat("  ", label, " Chisq: ", m2_m3_summary_chisq(chi), "\n", sep = "")
  cat("  ", label, " Df: ", summary_scalar(df), "\n", sep = "")
  cat("  ", label, " Pr(>Chisq): ", summary_p_value(p), "\n", sep = "")
}

m2_m3_summary_chisq <- function(x) {
  if (length(x) == 0L || is.na(x)) {
    return("NA")
  }
  digits <- max(3L, getOption("digits") - 3L)
  dig_tst <- max(1L, min(5L, digits - 1L))
  format(signif(x, dig_tst), digits = dig_tst)
}

#' @export
summary.gRm_m2 <- function(object, ...) {
  reject_summary_which(...)
  new_m2_m3_summary(object, include_m3 = FALSE)
}

#' @export
summary.gRm_m3 <- function(object, ...) {
  reject_summary_which(...)
  new_m2_m3_summary(object, include_m3 = TRUE)
}

new_m2_m3_summary <- function(object, include_m3) {
  values <- object$values %||% list()
  title <- if (isTRUE(include_m3)) "gRm: M3 fit diagnostic" else "gRm: M2 fit diagnostic"
  tables <- list(
    aggregates = m2_m3_aggregate_table(values, include_m3 = include_m3),
    item_trait = m2_m3_item_trait_table(values),
    invariance = m2_m3_invariance_table(values),
    margins = m2_m3_margin_summary_table(values, include_m3 = include_m3)
  )
  new_gRm_summary(
    object,
    title = title,
    which = names(tables),
    tables = tables,
    header = paste0("Score cuts: ", m2_m3_score_cuts_label(values$score_cuts %||% integer())),
    table_names = names(tables),
    print_table_names = FALSE
  )
}

m2_m3_selection_label <- function(values) {
  selected <- values$selected_items %||% data.frame()
  items <- selected$item_name %||% character()
  exogenous <- values$exogenous_names %||% character()
  paste(
    if (length(items)) paste(items, collapse = ", ") else "<none>",
    if (length(exogenous)) paste(exogenous, collapse = ", ") else "<none>",
    sep = " | "
  )
}

m2_m3_score_cuts_label <- function(score_cuts) {
  if (!length(score_cuts)) {
    return("<none>")
  }
  paste(as.integer(score_cuts), collapse = ", ")
}

m2_m3_aggregate_table <- function(values, include_m3) {
  rows <- if (isTRUE(include_m3)) {
    rbind_fill(values$m2 %||% data.frame(), values$m3 %||% data.frame())
  } else {
    values$aggregate %||% data.frame()
  }
  m2_m3_public_stat_table(rows, first_column = "Diagnostic")
}

m2_m3_item_trait_table <- function(values) {
  m2_m3_public_stat_table(values$item_trait %||% data.frame(), first_column = "Diagnostic")
}

m2_m3_invariance_table <- function(values) {
  table <- values$invariance %||% data.frame()
  if (!is.data.frame(table) || !nrow(table)) {
    return(data.frame())
  }
  data.frame(
    Exogenous = table$background_name,
    Chisq = table$chi_square,
    Df = as.integer(table$degrees_of_freedom),
    `Pr(>Chisq)` = table$p_value,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

m2_m3_margin_summary_table <- function(values, include_m3) {
  tests <- values$tests %||% data.frame()
  if (!isTRUE(include_m3) && is.data.frame(tests) && nrow(tests)) {
    tests <- tests[tests$is_m2, , drop = FALSE]
  }
  if (!is.data.frame(tests) || !nrow(tests)) {
    return(data.frame())
  }
  data.frame(
    Margin = tests$margin,
    Type = tests$margin_type,
    Item = m2_m3_blank_missing(tests$variable1),
    `Item/Exogenous` = m2_m3_blank_missing(tests$variable2),
    `Item/Exogenous/Score` = m2_m3_blank_missing(tests$variable3),
    Chisq = tests$chi_square,
    Df = as.integer(tests$degrees_of_freedom),
    `Pr(>Chisq)` = tests$p_value,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

m2_m3_blank_missing <- function(x) {
  x <- as.character(x %||% character())
  x[is.na(x)] <- ""
  x
}

m2_m3_public_stat_table <- function(table, first_column) {
  if (!is.data.frame(table) || !nrow(table)) {
    return(data.frame())
  }
  out <- data.frame(
    Diagnostic = table$diagnostic,
    Chisq = table$chi_square,
    Df = as.integer(table$degrees_of_freedom),
    `Pr(>Chisq)` = table$p_value,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(out)[[1L]] <- first_column
  out
}
