#' Internal yes no display helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @return The internal `yes_no_display()` computation result.
#' @keywords internal
#' @noRd
yes_no_display <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    ifelse(x, "yes", "no")
  )
}

#' Internal public bh marker note helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param bh Internal `bh` value used by this helper.
#' @param digits Internal `digits` value used by this helper.
#' @return The internal `public_bh_marker_note()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal public screen bh marker note helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param bh Internal `bh` value used by this helper.
#' @param digits Internal `digits` value used by this helper.
#' @return The internal `public_screen_bh_marker_note()` computation result.
#' @keywords internal
#' @noRd
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
    "For local dependence, * means retained in the final screen model; ",
    "provisional negative LD is named in the Decision column and excluded. ",
    "Local-dependence and DIF candidate evidence use the global ",
    "Benjamini-Hochberg threshold for FDR = ",
    fdr,
    " (threshold = ",
    format.pval(threshold, digits = dig_tst, eps = .Machine$double.eps),
    "); score effects use their source screening routine, and stricter global ",
    "cutoffs are available in attr(x, \"bh\")."
  )
}

#' Internal public value tables helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @return The internal `public_value_tables()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal public selected by bh helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param tests Diagnostic test rows.
#' @param threshold Internal `threshold` value used by this helper.
#' @return The internal `public_selected_by_bh()` computation result.
#' @keywords internal
#' @noRd
public_selected_by_bh <- function(tests, threshold) {
  if (!is.data.frame(tests) || !"p_value" %in% names(tests) || is.na(threshold)) {
    return(data.frame())
  }
  out <- tests[tests$p_value <= threshold, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Internal public selected rows helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param rows Rows used by the computation.
#' @return The internal `public_selected_rows()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal public format table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The internal `public_format_table()` computation result.
#' @keywords internal
#' @noRd
public_format_table <- function(table) {
  out <- table
  numeric_cols <- vapply(out, is.numeric, logical(1L))
  out[numeric_cols] <- lapply(out[numeric_cols], function(x) signif(x, 6L))
  out
}
