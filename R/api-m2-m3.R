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
#' @param bootstrap Single non-missing logical. If `TRUE`, run DIGRAM's
#'   source-shaped parametric bootstrap after the observed diagnostic. The
#'   default is `FALSE` so an ordinary diagnostic remains deterministic and
#'   inexpensive; DIGRAM's prompted runtime default is reproduced by the
#'   default `nsim = 100L` when bootstrap is enabled.
#' @param nsim Positive integer-like number of bootstrap samples. DIGRAM prompts
#'   with a default of 100.
#' @param seed `NULL` for a source-style nondeterministic initialization, or a
#'   positive integer-like seed in `1..2147483646` for the private reproducible
#'   validation stream. Supplying a seed does not modify R's global RNG state.
#' @param reestimate Single non-missing logical. If `TRUE` (the DIGRAM CM3
#'   command behavior), start a fresh fit of the same GLLRM in every bootstrap
#'   sample. If `FALSE`, evaluate samples under the supplied fitted parameters.
#' @param bootstrap_max_step Positive integer-like maximum number of GLLRM
#'   fitting steps for each bootstrap refit. The source value is 500.
#' @param keep_bootstrap_samples Single non-missing logical. Retain every
#'   generated zero-based response matrix and one-based exogenous matrix under
#'   `result$values$bootstrap$samples`. Raw per-replicate fit and statistic
#'   tables are always retained; response matrices default to being omitted.
#' @param resample_score_distribution Single non-missing logical corresponding
#'   to Pascal's `ResampleScoreDist`. The preserved source routine has an empty
#'   resampling branch, so both values keep the observed score distribution;
#'   the requested value and no-op mode are recorded in bootstrap metadata.
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
#' `score_cuts` changes only this diagnostic result.
#'
#' With `bootstrap = TRUE`, score counts are fixed within complete exogenous
#' strata and response patterns are generated conditionally from the supplied
#' fitted GLLRM. Active LD and DIF parameters participate in the same component
#' probabilities as fitting. Optional refits start from fresh source parameter
#' values and use `bootstrap_max_step`; only samples whose final sufficient-count
#' discrepancy is strictly less than 0.1 contribute. Bootstrap p-values are the
#' proportion of accepted simulated asymptotic p-values less than or equal to
#' the observed p-value. `nsim`, `nused`, every generated statistic, convergence
#' metadata, RNG state, and the individual exceedance counts are retained in
#' `result$values$bootstrap`.
#'
#' CM2/CM3 record eligibility is independent of the CML estimation score
#' window. Complete score-zero records contribute to observed and expected
#' diagnostic margins. The historical Pascal `Count_IJK` path also retains an
#' item-complete record when exogenous data are missing, while the corresponding
#' expected table and every other margin require complete exogenous values.
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
m2 <- function(fit,
               items = NULL,
               score_cuts = NULL,
               bootstrap = FALSE,
               nsim = 100L,
               seed = NULL,
               reestimate = TRUE,
               bootstrap_max_step = 500L,
               keep_bootstrap_samples = FALSE,
               resample_score_distribution = FALSE,
               ...) {
  reject_public_dots(...)
  fit <- as_public_gRm_fit(fit)
  bootstrap_control <- normalize_m2_m3_bootstrap_control(
    bootstrap,
    nsim,
    seed,
    reestimate,
    bootstrap_max_step,
    keep_bootstrap_samples,
    resample_score_distribution
  )
  values <- m2_values(
    fit,
    items = items,
    score_cuts = score_cuts,
    bootstrap_control = bootstrap_control
  )
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
m3 <- function(fit,
               items = NULL,
               score_cuts = NULL,
               bootstrap = FALSE,
               nsim = 100L,
               seed = NULL,
               reestimate = TRUE,
               bootstrap_max_step = 500L,
               keep_bootstrap_samples = FALSE,
               resample_score_distribution = FALSE,
               ...) {
  reject_public_dots(...)
  fit <- as_public_gRm_fit(fit)
  bootstrap_control <- normalize_m2_m3_bootstrap_control(
    bootstrap,
    nsim,
    seed,
    reestimate,
    bootstrap_max_step,
    keep_bootstrap_samples,
    resample_score_distribution
  )
  values <- m3_values(
    fit,
    items = items,
    score_cuts = score_cuts,
    bootstrap_control = bootstrap_control
  )
  new_m2_m3_result(fit, values, "gRm_m3", "m3", match.call())
}

#' Internal new m2 m3 result helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param fit Fitted gRm model.
#' @param values Values to validate or transform.
#' @param class Internal `class` value used by this helper.
#' @param result Result value to assemble or transform.
#' @param call Captured R call.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
new_m2_m3_result <- function(fit, values, class, result, call) {
  analysis <- fit$analysis %||% fit$spec$analysis
  selected <- values$selected_items %||% data.frame()
  metadata <- list(
    items = selected$item_name %||% character(),
    item_indices = selected$item_index %||% integer(),
    score_cuts = values$score_cuts %||% integer(),
    bootstrap = values$bootstrap[c(
      "enabled", "possible", "nsim", "nused", "seed", "reestimate",
      "acceptance_delta"
    )]
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
  reject_public_dots(...)
  print_m2_m3_result(x, include_m3 = FALSE)
}

#' @export
print.gRm_m3 <- function(x, ...) {
  reject_public_dots(...)
  print_m2_m3_result(x, include_m3 = TRUE)
}

#' Internal print m2 m3 result helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @param include_m3 Internal `include_m3` value used by this helper.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
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
  print_m2_m3_bootstrap_status(values$bootstrap)
  if (isTRUE(include_m3)) {
    cat("Use summary(x) to show two-way margins, three-way margins, and decompositions.\n")
  } else {
    cat("Use summary(x) to show two-way margins and invariance by exogenous variable.\n")
  }
  invisible(x)
}

#' Internal print m2 m3 stat block helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param label Internal `label` value used by this helper.
#' @param row Internal `row` value used by this helper.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
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
  if (is.data.frame(row) && nrow(row) && "bootstrap_p_value" %in% names(row)) {
    cat(
      "  ", label, " Bootstrap Pr: ",
      summary_p_value(row$bootstrap_p_value[[1L]]),
      "\n",
      sep = ""
    )
  }
}

#' Internal print m2 m3 bootstrap status helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param bootstrap Whether to run bootstrap calibration.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
print_m2_m3_bootstrap_status <- function(bootstrap) {
  if (!is.list(bootstrap) || !isTRUE(bootstrap$enabled)) {
    cat("  Parametric bootstrap: not requested\n\n")
    return(invisible(NULL))
  }
  if (!isTRUE(bootstrap$possible)) {
    cat("  Parametric bootstrap: not possible\n")
    if (length(bootstrap$reasons)) {
      cat("  Reason: ", paste(bootstrap$reasons, collapse = "; "), "\n", sep = "")
    }
    cat("\n")
    return(invisible(NULL))
  }
  cat(
    "  Parametric bootstrap: ", bootstrap$nused, "/", bootstrap$nsim,
    " samples accepted (delta < ", bootstrap$acceptance_delta, ")\n\n",
    sep = ""
  )
  invisible(NULL)
}

#' Internal m2 m3 summary chisq helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @return The internal `m2_m3_summary_chisq()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal new m2 m3 summary helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @param include_m3 Internal `include_m3` value used by this helper.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
new_m2_m3_summary <- function(object, include_m3) {
  values <- object$values %||% list()
  title <- if (isTRUE(include_m3)) "gRm: M3 fit diagnostic" else "gRm: M2 fit diagnostic"
  tables <- list(
    aggregates = m2_m3_aggregate_table(values, include_m3 = include_m3),
    item_trait = m2_m3_item_trait_table(values),
    invariance = m2_m3_invariance_table(values),
    margins = m2_m3_margin_summary_table(values, include_m3 = include_m3)
  )
  if (isTRUE(values$bootstrap$enabled)) {
    tables$bootstrap <- m2_m3_bootstrap_summary_table(values$bootstrap)
  }
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

#' Internal m2 m3 selection label helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @return The internal `m2_m3_selection_label()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal m2 m3 score cuts label helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param score_cuts Resolved total-score cut values.
#' @return The internal `m2_m3_score_cuts_label()` computation result.
#' @keywords internal
#' @noRd
m2_m3_score_cuts_label <- function(score_cuts) {
  if (!length(score_cuts)) {
    return("<none>")
  }
  paste(as.integer(score_cuts), collapse = ", ")
}

#' Internal m2 m3 aggregate table helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @param include_m3 Internal `include_m3` value used by this helper.
#' @return The internal `m2_m3_aggregate_table()` computation result.
#' @keywords internal
#' @noRd
m2_m3_aggregate_table <- function(values, include_m3) {
  rows <- if (isTRUE(include_m3)) {
    rbind_fill(values$m2 %||% data.frame(), values$m3 %||% data.frame())
  } else {
    values$aggregate %||% data.frame()
  }
  m2_m3_public_stat_table(rows, first_column = "Diagnostic")
}

#' Internal m2 m3 item trait table helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @return The internal `m2_m3_item_trait_table()` computation result.
#' @keywords internal
#' @noRd
m2_m3_item_trait_table <- function(values) {
  m2_m3_public_stat_table(values$item_trait %||% data.frame(), first_column = "Diagnostic")
}

#' Internal m2 m3 invariance table helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @return The internal `m2_m3_invariance_table()` computation result.
#' @keywords internal
#' @noRd
m2_m3_invariance_table <- function(values) {
  table <- values$invariance %||% data.frame()
  if (!is.data.frame(table) || !nrow(table)) {
    return(data.frame())
  }
  out <- data.frame(
    Exogenous = table$background_name,
    Chisq = table$chi_square,
    Df = as.integer(table$degrees_of_freedom),
    `Pr(>Chisq)` = table$p_value,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if ("bootstrap_p_value" %in% names(table)) {
    out$`Bootstrap Pr` <- table$bootstrap_p_value
    out$`Bootstrap n` <- as.integer(table$bootstrap_nused)
  }
  out
}

#' Internal m2 m3 margin summary table helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @param include_m3 Internal `include_m3` value used by this helper.
#' @return The internal `m2_m3_margin_summary_table()` computation result.
#' @keywords internal
#' @noRd
m2_m3_margin_summary_table <- function(values, include_m3) {
  tests <- values$tests %||% data.frame()
  if (!isTRUE(include_m3) && is.data.frame(tests) && nrow(tests)) {
    tests <- tests[tests$is_m2, , drop = FALSE]
  }
  if (!is.data.frame(tests) || !nrow(tests)) {
    return(data.frame())
  }
  out <- data.frame(
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
  if ("bootstrap_p_value" %in% names(tests)) {
    out$`Bootstrap Pr` <- tests$bootstrap_p_value
    out$`Bootstrap n` <- as.integer(tests$bootstrap_nused)
  }
  out
}

#' Internal m2 m3 blank missing helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @return The internal `m2_m3_blank_missing()` computation result.
#' @keywords internal
#' @noRd
m2_m3_blank_missing <- function(x) {
  x <- as.character(x %||% character())
  x[is.na(x)] <- ""
  x
}

#' Internal m2 m3 public stat table helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @param first_column Internal `first_column` value used by this helper.
#' @return The internal `m2_m3_public_stat_table()` computation result.
#' @keywords internal
#' @noRd
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
  if ("bootstrap_p_value" %in% names(table)) {
    out$`Bootstrap Pr` <- table$bootstrap_p_value
    out$`Bootstrap n` <- as.integer(table$bootstrap_nused)
  }
  names(out)[[1L]] <- first_column
  out
}

#' Internal m2 m3 bootstrap summary table helper
#'
#' Supports the api m2 m3 implementation while preserving its internal contract.
#' @param bootstrap Whether to run bootstrap calibration.
#' @return The internal `m2_m3_bootstrap_summary_table()` computation result.
#' @keywords internal
#' @noRd
m2_m3_bootstrap_summary_table <- function(bootstrap) {
  if (!is.list(bootstrap) || !isTRUE(bootstrap$enabled)) {
    return(data.frame())
  }
  data.frame(
    Status = bootstrap$source_status %||% NA_character_,
    Possible = bootstrap$possible %||% NA,
    Simulations = as.integer(bootstrap$nsim %||% 0L),
    Accepted = as.integer(bootstrap$nused %||% 0L),
    Seed = as.integer(bootstrap$seed %||% NA_integer_),
    Reestimated = bootstrap$reestimate %||% NA,
    `Acceptance delta` = bootstrap$acceptance_delta %||% 0.1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
