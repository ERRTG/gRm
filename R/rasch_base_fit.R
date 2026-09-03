#' Count observed item scores and total scores
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param bundle An item-parameters bundle from
#'   `build_item_parameters_bundle()`.
#' @return A list containing the number of valid rows, an item-by-score count
#'   matrix, and a total-score count vector.
#' @keywords internal
#' @noRd
rasch_counts <- function(bundle) {
  # Source trace: pascal_harness/SourceRaschCore.pas::LoadCounts.
  # Only status==1 rows contribute to the CML item margins and total-score
  # margins used by the base Rasch iteration.
  items <- bundle$model$items
  data <- bundle$data
  valid <- data$status == 1L
  max_total_score <- bundle$model$max_total_score

  item_counts <- matrix(
    0L,
    nrow = nrow(items),
    ncol = max(items$raw_max),
    dimnames = list(items$name, as.character(seq_len(max(items$raw_max)) - 1L))
  )
  for (item_index in seq_len(nrow(items))) {
    values <- data[[items$name[[item_index]]]][valid]
    tab <- table(factor(values, levels = seq.int(0L, items$raw_max[[item_index]] - 1L)))
    item_counts[item_index, seq_along(tab)] <- as.integer(tab)
  }

  score_counts <- integer(max_total_score + 1L)
  names(score_counts) <- as.character(seq.int(0L, max_total_score))
  score_tab <- table(factor(data$score[valid], levels = seq.int(0L, max_total_score)))
  score_counts[] <- as.integer(score_tab)

  list(
    n_valid = sum(valid),
    item_counts = item_counts,
    score_counts = score_counts
  )
}

#' Initialize item gamma parameters
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param bundle An item-parameters bundle from
#'   `build_item_parameters_bundle()`.
#' @return A matrix of initial multiplicative item score parameters, with valid
#'   categories initialized to one and unused cells set to zero.
#' @keywords internal
#' @noRd
initial_item_gamma <- function(bundle) {
  # Source trace: pascal_harness/SourceRaschCore.pas::InitializeRaschFit.
  # DIGRAM starts all observed item-score gamma parameters at 1 and leaves
  # unused source-bounded cells at 0.
  items <- bundle$model$items
  gamma <- matrix(
    0,
    nrow = nrow(items),
    ncol = max(items$raw_max),
    dimnames = list(items$name, as.character(seq_len(max(items$raw_max)) - 1L))
  )
  for (item_index in seq_len(nrow(items))) {
    gamma[item_index, seq_len(items$raw_max[[item_index]])] <- 1
  }
  gamma
}

#' Build the score-generating function excluding one item
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param bundle An item-parameters bundle.
#' @param item_gamma Matrix of current item gamma parameters.
#' @param excluded_item One-based item index to exclude from the convolution.
#' @return A numeric vector indexed by total score plus one.
#' @keywords internal
#' @noRd
build_gamma_excluding_item <- function(bundle, item_gamma, excluded_item) {
  # Source trace: pascal_harness/SourceRaschCore.pas::BuildGammaExcludingItem.
  items <- bundle$model$items
  max_total_score <- bundle$model$max_total_score
  gamma_values <- numeric(max_total_score + 1L)
  gamma_values[[1]] <- 1
  current_max <- 0L

  for (item_index in seq_len(nrow(items))) {
    if (item_index == excluded_item) {
      next
    }
    # Convolve the score polynomial for all included items:
    # Gamma_new[s+x] += Gamma_old[s] * item_gamma[item, x].
    next_values <- numeric(max_total_score + 1L)
    next_max <- current_max + items$raw_max[[item_index]] - 1L
    for (score in seq.int(0L, current_max)) {
      current_weight <- gamma_values[[score + 1L]]
      if (current_weight != 0) {
        for (item_score in seq.int(0L, items$raw_max[[item_index]] - 1L)) {
          next_values[[score + item_score + 1L]] <-
            next_values[[score + item_score + 1L]] +
            current_weight * item_gamma[item_index, as.character(item_score)]
        }
      }
    }
    gamma_values[seq_len(next_max + 1L)] <- next_values[seq_len(next_max + 1L)]
    current_max <- next_max
  }

  gamma_values
}

#' Build the DIGRAM source score gamma array
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param bundle An item-parameters bundle.
#' @param item_gamma Matrix of current item gamma parameters.
#' @param use_items Logical vector selecting items included in the score gamma.
#' @return A numeric vector indexed by total score plus one.
#' @keywords internal
#' @noRd
build_source_score_gamma <- function(bundle, item_gamma, use_items) {
  # Source trace: source/digram_source_20260817/skunits/skbias12.pas::
  # Inexpensive_Gamma_Calculation. NewGamma starts as Sgamma, each included
  # valid item adds shifted Sgamma terms across the full source score limits,
  # and the old Sgamma is then subtracted. The report-driving Delphi copy still
  # runs the item loop for single-category valid items; when the source scale
  # adjustment has set that sole category to one, the add/subtract pair is a
  # no-op but the loop shape stays source-faithful.
  items <- bundle$model$items
  max_total_score <- bundle$model$max_total_score
  use_items <- as.logical(use_items)
  if (length(use_items) != nrow(items)) {
    stop("`use_items` must have one entry per item.", call. = FALSE)
  }

  ifra <- integer(nrow(items))
  itil <- integer(nrow(items))
  sfra <- 0L
  stil <- 0L

  for (item_index in seq_len(nrow(items))) {
    scores <- seq.int(0L, items$raw_max[[item_index]] - 1L)
    gamma_values <- as.numeric(item_gamma[item_index, as.character(scores)])
    positive <- scores[gamma_values > 0]
    if (length(positive) == 0L) {
      ifra[[item_index]] <- items$raw_max[[item_index]]
      itil[[item_index]] <- 0L
    } else {
      ifra[[item_index]] <- min(positive)
      itil[[item_index]] <- max(positive)
      if (isTRUE(use_items[[item_index]]) && ifra[[item_index]] <= itil[[item_index]]) {
        sfra <- sfra + ifra[[item_index]]
        stil <- stil + itil[[item_index]]
      }
    }
  }

  gamma_values <- numeric(max_total_score + 1L)
  gamma_values[[sfra + 1L]] <- 1
  next_values <- gamma_values

  for (item_index in seq_len(nrow(items))) {
    if (!isTRUE(use_items[[item_index]]) || ifra[[item_index]] > itil[[item_index]]) {
      next
    }

    for (item_score in seq.int(ifra[[item_index]], itil[[item_index]])) {
      offset <- item_score - ifra[[item_index]]
      if (sfra <= stil - 1L) {
        for (score in seq.int(sfra, stil - 1L)) {
          target_score <- score + offset
          if (target_score <= max_total_score) {
            next_values[[target_score + 1L]] <-
              next_values[[target_score + 1L]] +
              item_gamma[item_index, as.character(item_score)] *
                gamma_values[[score + 1L]]
          }
        }
      }
    }

    for (score in seq.int(sfra, stil)) {
      next_values[[score + 1L]] <- next_values[[score + 1L]] - gamma_values[[score + 1L]]
    }

    gamma_values <- next_values
  }

  gamma_values
}

#' Calculate expected item score counts under the current Rasch fit
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param bundle An item-parameters bundle.
#' @param counts Count list from `rasch_counts()`.
#' @param item_gamma Matrix of current item gamma parameters.
#' @return A matrix of expected item score counts with the same shape as
#'   `item_gamma`.
#' @keywords internal
#' @noRd
calculate_rasch_expected_items <- function(bundle, counts, item_gamma) {
  # Source trace: pascal_harness/SourceRaschCore.pas::CalculateRaschExpectedItems.
  items <- bundle$model$items
  expected <- item_gamma
  expected[,] <- 0
  use_items <- rep(TRUE, nrow(items))
  full_gamma <- build_source_score_gamma(bundle, item_gamma, use_items)

  for (item_index in seq_len(nrow(items))) {
    # WithoutItem is the rest-score gamma table. The denominator is the source
    # full-score gamma, not a freshly summed item-specific denominator.
    use_items[[item_index]] <- FALSE
    without_item <- build_source_score_gamma(bundle, item_gamma, use_items)
    use_items[[item_index]] <- TRUE
    for (score in seq.int(0L, bundle$model$max_total_score)) {
      score_count <- counts$score_counts[[score + 1L]]
      if (score_count == 0L) {
        next
      }
      denominator <- full_gamma[[score + 1L]]
      if (denominator <= 0) {
        next
      }
      ratio <- score_count / denominator
      for (item_score in seq.int(0L, items$raw_max[[item_index]] - 1L)) {
        if (score >= item_score) {
          numerator <- item_gamma[item_index, as.character(item_score)] *
            without_item[[score - item_score + 1L]]
          # Multiply the conditional probability by the observed count at this
          # total score to get the expected item-score margin contribution.
          expected[item_index, as.character(item_score)] <-
            expected[item_index, as.character(item_score)] +
            numerator * ratio
        }
      }
    }
  }

  expected
}

#' Calculate IPF update ratios for item gamma parameters
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param bundle An item-parameters bundle.
#' @param counts Count list from `rasch_counts()`.
#' @param expected Expected item score counts from
#'   `calculate_rasch_expected_items()`.
#' @param item_gamma Matrix of current item gamma parameters.
#' @param apply_update Logical; if `TRUE`, multiply the current gamma values by
#'   the observed/fitted ratios.
#' @return A list containing the next gamma matrix, the update-ratio matrix, and
#'   the maximum observed/fitted count discrepancy.
#' @keywords internal
#' @noRd
calculate_rasch_update_ratios <- function(bundle, counts, expected, item_gamma, apply_update) {
  # Source trace:
  # pascal_harness/SourceRaschCore.pas::CalculateRaschUpdateRatiosFromScore
  # and ::CalculateRaschUpdateRatios.
  items <- bundle$model$items
  update <- item_gamma
  update[,] <- 1
  delta <- 0
  next_gamma <- item_gamma

  for (item_index in seq_len(nrow(items))) {
    for (item_score in seq.int(0L, items$raw_max[[item_index]] - 1L)) {
      observed <- counts$item_counts[item_index, as.character(item_score)]
      fitted <- expected[item_index, as.character(item_score)]
      if (observed > 0) {
        if (fitted > 0) {
          # IPF multiplier for the current sufficient-statistic cell.
          ratio <- observed / fitted
        } else {
          ratio <- 1
        }
        # DIGRAM's convergence delta is the maximum absolute observed/fitted
        # count discrepancy, not the maximum ratio movement.
        cell_delta <- abs(fitted - observed)
        if (cell_delta > delta) {
          delta <- cell_delta
        }
      } else {
        ratio <- 0
      }
      update[item_index, as.character(item_score)] <- ratio
      if (apply_update) {
        next_gamma[item_index, as.character(item_score)] <-
          next_gamma[item_index, as.character(item_score)] * ratio
      }
    }
  }

  list(item_gamma = next_gamma, update = update, delta = delta)
}

#' Rescale item gamma parameters to the Pascal source convention
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param bundle An item-parameters bundle.
#' @param item_gamma Matrix of item gamma parameters after an IPF update.
#' @return A rescaled item gamma matrix using the source normalization.
#' @keywords internal
#' @noRd
adjust_item_gammas_source_scale <- function(bundle, item_gamma) {
  # Source trace: pascal_harness/SourceRaschCore.pas::AdjustItemGammasSourceScale.
  items <- bundle$model$items
  ifra <- integer(nrow(items))
  itil <- integer(nrow(items))
  valid_item <- logical(nrow(items))

  for (item_index in seq_len(nrow(items))) {
    scores <- seq.int(0L, items$raw_max[[item_index]] - 1L)
    gamma_values <- as.numeric(item_gamma[item_index, as.character(scores)])
    positive <- scores[gamma_values > 0]
    if (length(positive) == 0L) {
      ifra[[item_index]] <- items$raw_max[[item_index]]
      itil[[item_index]] <- 0L
    } else {
      ifra[[item_index]] <- min(positive)
      itil[[item_index]] <- max(positive)
      valid_item[[item_index]] <- ifra[[item_index]] < itil[[item_index]]
    }
  }

  last_sgamma <- 1
  s_max <- 0L
  s_min <- 0L

  for (item_index in seq_len(nrow(items))) {
    if (!isTRUE(valid_item[[item_index]])) {
      next
    }
    fra <- ifra[[item_index]]
    til <- itil[[item_index]]
    s_max <- s_max + til
    s_min <- s_min + fra
    alpha <- item_gamma[item_index, as.character(fra)]
    if (alpha > 0) {
      # First set each item's lowest-score gamma to 1 by dividing through by
      # gamma[item, fra].
      for (score in seq.int(fra, til)) {
        item_gamma[item_index, as.character(score)] <-
          item_gamma[item_index, as.character(score)] / alpha
      }
    }
    if (item_gamma[item_index, as.character(til)] > 0) {
      last_sgamma <- last_sgamma * item_gamma[item_index, as.character(til)]
    }
  }

  alpha <- 0
  if ((s_max - s_min) > 0) {
    # Then distribute the remaining product of highest-score gammas across the
    # score scale so the source normalization is stable after each IPF step.
    alpha <- -log(last_sgamma) / (s_max - s_min)
  }

  for (item_index in seq_len(nrow(items))) {
    if (!isTRUE(valid_item[[item_index]])) {
      next
    }
    fra <- ifra[[item_index]]
    til <- itil[[item_index]]
    for (score in seq.int(fra, til)) {
      item_gamma[item_index, as.character(score)] <-
        exp((score - fra) * alpha) * item_gamma[item_index, as.character(score)]
    }
  }
  for (item_index in seq_len(nrow(items))) {
    if (!isTRUE(valid_item[[item_index]]) && ifra[[item_index]] <= itil[[item_index]]) {
      item_gamma[item_index, as.character(itil[[item_index]])] <- 1
    }
  }

  item_gamma
}

#' Fit the base Rasch item parameters
#'
#' Fits the source-shaped base Rasch model used by the current DIGRAM
#' item-parameters report slice. The implementation mirrors the Pascal IPF/CML
#' update structure and returns both fitted parameters and diagnostic matrices
#' used by downstream parity tests.
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param bundle An item-parameters bundle from
#'   `build_item_parameters_bundle()`.
#' @param max_step Maximum number of IPF iterations.
#' @param max_delta Convergence threshold for the maximum item-score count
#'   discrepancy.
#' @return A list with components:
#'   \describe{
#'     \item{`n_step`}{Number of IPF iterations performed.}
#'     \item{`delta`}{Final maximum observed/fitted count discrepancy.}
#'     \item{`report_delta`}{Delta value reported at convergence, matching the
#'       Pascal report convention.}
#'     \item{`converged`}{Logical convergence flag.}
#'     \item{`item_gamma`}{Fitted multiplicative item score parameters.}
#'     \item{`expected_items`}{Expected item score counts at the fitted values.}
#'     \item{`update_items`}{Final observed/fitted update ratios.}
#'     \item{`counts`}{Observed count summaries used for fitting.}
#'   }
#' @examples
#' \dontrun{
#' project <- read_digram_project("path/to/DIGRAM")
#' bundle <- build_item_parameters_bundle(project)
#' fit <- fit_rasch_base(bundle)
#' fit$converged
#' }
#' @keywords internal
#' @noRd
fit_rasch_base <- function(bundle,
                           max_step = 5000L,
                           max_delta = 0.0001) {
  # Source trace: pascal_harness/SourceRaschCore.pas::EstimateBaseRaschParameters
  # base_Rasch_only branch.
  counts <- rasch_counts(bundle)
  item_gamma <- initial_item_gamma(bundle)
  n_step <- 0L
  report_delta <- 0
  converged <- FALSE
  expected <- item_gamma
  update <- item_gamma
  update[,] <- 1
  delta <- 0

  if (counts$n_valid == 0L || nrow(bundle$model$items) == 0L) {
    return(list(
      n_step = n_step,
      delta = delta,
      report_delta = report_delta,
      converged = TRUE,
      max_step = max_step,
      max_delta = max_delta,
      item_gamma = item_gamma,
      expected_items = expected,
      update_items = update,
      counts = counts
    ))
  }

  while (n_step < max_step) {
    n_step <- n_step + 1L
    expected <- calculate_rasch_expected_items(bundle, counts, item_gamma)
    # The source calculates ratios once to set delta/update diagnostics, then
    # again with ApplyUpdate=TRUE before source-scale normalization.
    ratio_state <- calculate_rasch_update_ratios(bundle, counts, expected, item_gamma, FALSE)
    delta <- ratio_state$delta
    update <- ratio_state$update
    ratio_state <- calculate_rasch_update_ratios(bundle, counts, expected, item_gamma, TRUE)
    item_gamma <- adjust_item_gammas_source_scale(bundle, ratio_state$item_gamma)
    if (delta < max_delta) {
      converged <- TRUE
      report_delta <- delta
      break
    }
  }

  if (!converged || report_delta == 0) {
    report_delta <- delta
  }

  expected <- calculate_rasch_expected_items(bundle, counts, item_gamma)
  # Final expected margins and update ratios are recomputed after the last
  # parameter update, matching the diagnostics emitted by the Pascal harness.
  ratio_state <- calculate_rasch_update_ratios(bundle, counts, expected, item_gamma, FALSE)
  delta <- ratio_state$delta
  update <- ratio_state$update
  converged <- delta < max_delta

  list(
    n_step = n_step,
    delta = delta,
    report_delta = report_delta,
    converged = converged,
    max_step = max_step,
    max_delta = max_delta,
    item_gamma = item_gamma,
    expected_items = expected,
    update_items = update,
    counts = counts
  )
}
