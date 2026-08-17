#' Source support summary for public fitting
#'
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param bundle Source-shaped estimation bundle.
#' @return Item support table and observed Rasch counts.
#' @keywords internal
source_public_fit_item_support <- function(bundle) {
  counts <- rasch_counts(bundle)
  items <- bundle$model$items
  rows <- lapply(seq_len(nrow(items)), function(item) {
    declared <- seq.int(0L, items$raw_max[[item]] - 1L)
    observed <- declared[counts$item_counts[item, as.character(declared)] > 0L]
    data.frame(
      item_index = item,
      item = items$name[[item]],
      declared_categories = length(declared),
      supported_categories = length(observed),
      supported_scores = paste(observed, collapse = ","),
      complete_declared_support = length(observed) == length(declared),
      valid_item = length(observed) >= 2L,
      stringsAsFactors = FALSE
    )
  })
  list(counts = counts, items = do.call(rbind, rows))
}

#' Source reference-free cells for one interaction margin
#'
#' The source chooses the row and column having the largest number of positive
#' cells, then estimates positive cells outside that reference row and column.
#'
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param observed Observed IJ or IX margin.
#' @return Two-column integer matrix of one-based free-cell indices.
#' @keywords internal
source_public_fit_free_cells <- function(observed) {
  row_support <- rowSums(observed > 0)
  column_support <- colSums(observed > 0)
  if (!any(row_support > 0) || !any(column_support > 0)) {
    return(matrix(integer(), nrow = 0L, ncol = 2L))
  }
  row_reference <- which.max(row_support)
  column_reference <- which.max(column_support)
  cells <- which(observed > 0, arr.ind = TRUE)
  cells[
    cells[, 1L] != row_reference & cells[, 2L] != column_reference,
    ,
    drop = FALSE
  ]
}

#' Build the source-counted public fit support
#'
#' Source trace: `skbias12b.pas::Count_Margins`, `InitializeParameters`, and
#' `calculate_Nparameters`. The source counts supported item parameters and
#' positive reference-free IJ/IX cells; it does not infer estimability from the
#' affine rank of only the response patterns observed in the sample. Such an
#' empirical-pattern rank is not a valid conditional-likelihood support test:
#' a sufficient-statistic mean can be interior to the full conditional support
#' even when only a small number of distinct patterns has positive frequency.
#'
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param model Canonical `gRm_model` specification.
#' @param bundle Source-shaped estimation bundle.
#' @return A list with the source parameter count, per-term support counts, and
#'   optional GLLRM context.
#' @keywords internal
source_public_fit_support <- function(model, bundle) {
  has_dependencies <- nrow(model$ld) > 0L || nrow(model$dif) > 0L
  context <- if (has_dependencies) build_gllrm_context(model, bundle) else NULL
  counts <- rasch_counts(bundle)
  term_counts <- data.frame(
    type = character(),
    term = character(),
    n_parameters = integer(),
    stringsAsFactors = FALSE
  )
  if (!is.null(context)) {
    for (ld_index in seq_along(context$ld_specs)) {
      spec <- context$ld_specs[[ld_index]]
      free <- source_public_fit_free_cells(context$observed_ld[[ld_index]])
      term_counts <- rbind(
        term_counts,
        data.frame(type = "LD", term = spec$term, n_parameters = nrow(free))
      )
    }
    for (dif_index in seq_along(context$dif_specs)) {
      spec <- context$dif_specs[[dif_index]]
      free <- source_public_fit_free_cells(context$observed_dif[[dif_index]])
      term_counts <- rbind(
        term_counts,
        data.frame(type = "DIF", term = spec$term, n_parameters = nrow(free))
      )
    }
  }
  # Mathematical step: reproduce calculate_Nparameters exactly after the
  # source has selected positive item and reference-free interaction cells.
  expected_parameters <- if (is.null(context)) {
    calculate_source_n_parameters(counts$item_counts)
  } else {
    calculate_gllrm_n_parameters(context)
  }
  list(
    n_parameters = as.integer(expected_parameters),
    term_counts = term_counts,
    context = context
  )
}

#' Validate public Rasch/GLLRM estimability
#'
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param bundle Source-shaped estimation bundle.
#' @param model Canonical `gRm_model` specification.
#' @return `bundle`, invisibly.
#' @keywords internal
assert_public_estimable_fit_bundle <- function(bundle, model) {
  if (isTRUE(bundle$manifest$nvalid < 1L)) {
    stop(
      paste(
        "Rasch fitting requires at least one source-valid complete response pattern",
        "inside the DIGRAM score window."
      ),
      call. = FALSE
    )
  }

  support <- source_public_fit_item_support(bundle)
  invalid_items <- support$items$item[!support$items$valid_item]
  if (length(invalid_items)) {
    stop(
      "Rasch fitting requires at least two supported response categories for every item; insufficient support for: ",
      paste(invalid_items, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  incomplete_items <- support$items$item[!support$items$complete_declared_support]
  if (length(incomplete_items)) {
    stop(
      "Rasch fitting requires every declared item category to occur in the source estimation rows; missing category support for: ",
      paste(incomplete_items, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  estimability <- source_public_fit_support(model, bundle)
  zero_terms <- estimability$term_counts$term[estimability$term_counts$n_parameters < 1L]
  if (length(zero_terms)) {
    stop(
      "The requested GLLRM contains interaction terms without estimable observed support: ",
      paste(zero_terms, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!is.finite(estimability$n_parameters) || estimability$n_parameters < 1L) {
    stop(
      "Rasch fitting requires at least one non-reference parameter after source support rules.",
      call. = FALSE
    )
  }
  invisible(bundle)
}

#' Validate values before constructing a public fit
#'
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param values Source-shaped fitted value list.
#' @return `values`, invisibly.
#' @keywords internal
assert_public_fit_values <- function(values) {
  n_parameters <- values$n_parameters
  if (
    length(n_parameters) != 1L ||
      !is.numeric(n_parameters) ||
      is.na(n_parameters) ||
      !is.finite(n_parameters) ||
      n_parameters < 0 ||
      n_parameters != floor(n_parameters)
  ) {
    stop("The fitted model produced an invalid source parameter count.", call. = FALSE)
  }
  negative_log_likelihood <- values$log_likelihood
  if (
    length(negative_log_likelihood) != 1L ||
      !is.numeric(negative_log_likelihood) ||
      is.na(negative_log_likelihood) ||
      !is.finite(negative_log_likelihood)
  ) {
    stop("The fitted model produced an invalid negative conditional log likelihood.", call. = FALSE)
  }

  threshold_values <- as.numeric(values$thresholds %||% numeric())
  finite_thresholds <- threshold_values[is.finite(threshold_values)]
  if (any(finite_thresholds <= -999999)) {
    stop("The fitted model produced a source sentinel threshold and is not a valid public fit.", call. = FALSE)
  }
  item_gamma <- as.numeric(values$item_gamma %||% numeric())
  if (!length(item_gamma) || any(!is.finite(item_gamma))) {
    stop("The fitted model produced invalid item parameters.", call. = FALSE)
  }
  invisible(values)
}
