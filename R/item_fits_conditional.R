#' Internal calculate conditional item fit values helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle Source-shaped analysis bundle.
#' @param fit Fitted gRm model.
#' @param conditional Internal `conditional` value used by this helper.
#' @param incomplete Internal `incomplete` value used by this helper.
#' @return The internal `calculate_conditional_item_fit_values()` computation result.
#' @keywords internal
#' @noRd
calculate_conditional_item_fit_values <- function(bundle, fit, conditional = NULL, incomplete = NULL) {
  if (is.null(conditional)) {
    conditional <- item_conditional_moments(bundle, fit$item_gamma, include_probabilities = TRUE)
  }
  accumulators <- conditional_item_fit_complete_accumulators(bundle, conditional)
  if (is.null(incomplete)) {
    incomplete <- empty_source_incomplete_records(bundle)
  }
  accumulators <- conditional_item_fit_add_incomplete(
    accumulators,
    bundle,
    fit$item_gamma,
    incomplete
  )

  outfit <- safe_ratio(accumulators$outfit_sum, accumulators$n_used)
  outfit_expected <- safe_ratio(accumulators$outfit_mean, accumulators$n_used)
  outfit_variance <- safe_ratio(accumulators$outfit_var, accumulators$n_used^2)
  outfit_z <- safe_z(outfit, outfit_expected, outfit_variance)

  infit <- safe_ratio(accumulators$infit_sum, accumulators$infit_weight)
  infit_expected <- safe_ratio(accumulators$infit_mean, accumulators$infit_weight)
  infit_variance <- safe_ratio(accumulators$infit_var, accumulators$infit_weight^2)
  infit_z <- safe_z(infit, infit_expected, infit_variance)

  list(
    outfit = outfit,
    outfit_sd = sqrt(pmax(outfit_variance, 0)),
    p_outfit = two_sided_source_normal_p(outfit_z),
    infit = inf_replace(infit, 1),
    infit_sd = sqrt(pmax(infit_variance, 0)),
    p_infit = two_sided_source_normal_p(infit_z)
  )
}

#' Accumulate complete-record conditional item-fit contributions
#'
#' Source trace: `source/PAS_skunits/skbias15.pas::CalculateInAndOutfits`.
#' Mathematical step: traverse source score groups then items and accumulate
#' score-count-weighted outfit and variance-weighted infit moments.
#' @param bundle Source-shaped analysis bundle.
#' @param conditional Conditional item moments by item and total score.
#' @return Named item-level outfit/infit accumulators.
#' @keywords internal
#' @noRd
conditional_item_fit_complete_accumulators <- function(bundle, conditional) {
  items <- bundle$model$items
  data <- bundle$data
  n_items <- nrow(items)
  outfit_sum <- outfit_mean <- outfit_var <- numeric(n_items)
  infit_sum <- infit_mean <- infit_var <- infit_weight <- numeric(n_items)
  n_used <- integer(n_items)

  least_score <- bundle$model$least_score
  if (is.null(least_score)) {
    least_score <- 0L
  }
  largest_score <- bundle$model$largest_score
  if (is.null(largest_score)) {
    largest_score <- bundle$model$max_total_score
  }
  from_score <- if (least_score == 0L) 1L else least_score
  to_score <- if (largest_score == bundle$model$max_total_score) {
    bundle$model$max_total_score - 1L
  } else {
    largest_score
  }

  item_matrix <- data[, items$name, drop = FALSE]
  complete_items <- apply(item_matrix >= 0L, 1L, all)
  item_scores <- rowSums(item_matrix * (item_matrix >= 0L))
  score_weights <- tabulate(
    item_scores[complete_items] + 1L,
    nbins = bundle$model$max_total_score + 1L
  )
  observed_rows <- complete_items
  if (nrow(bundle$model$backgrounds) > 0L) {
    observed_rows <- observed_rows & data$status == 1L
  }

  for (score in seq.int(from_score, to_score)) {
    score_weight <- score_weights[[score + 1L]]
    if (score_weight <= 0L) {
      next
    }
    score_rows <- observed_rows & item_scores == score
    for (item_index in seq_len(n_items)) {
      item_values <- data[[items$name[[item_index]]]][score_rows]
      observed_count <- tabulate(item_values + 1L, nbins = items$raw_max[[item_index]])
      observed_total <- sum(observed_count)
      if (observed_total <= 0) {
        # Source trace: if Count_Observed has no valid-background item margin
        # for this score, CalculateOutfit leaves the score contribution at
        # zero, but the later finalization still divides outfit by the broader
        # complete-item ScoreDistribution. Infit receives no variance weight.
        n_used[[item_index]] <- n_used[[item_index]] + score_weight
        next
      }
      moments <- conditional[[item_index]][[score + 1L]]
      probabilities <- moments$probabilities
      variance <- moments$variance
      if (variance <= 0) {
        next
      }

      item_scores_for_item <- seq.int(0L, items$raw_max[[item_index]] - 1L)
      centered <- item_scores_for_item - moments$mean
      outfit_range <- centered^2 / variance
      infit_range <- centered^2
      observed_frequency <- observed_count / observed_total

      # Source trace: skbias15.pas::Count_Observed applies GET_EXOGENE before
      # filling ItemMargTables, so observed item-score proportions can be based
      # on valid-background rows. CalculateInAndOutfits then weights those
      # proportions with the broader ScoreDistribution over complete item rows.
      # That mixes conditioning populations and is mathematically unusual, but
      # reproduces DIGRAM's ItemFits output.
      outfit_score <- sum(outfit_range * observed_frequency)
      outfit_expected <- sum(outfit_range * probabilities)
      outfit_score_variance <- sum(outfit_range^2 * probabilities) - outfit_expected^2

      infit_score <- sum(infit_range * observed_frequency)
      infit_expected <- sum(infit_range * probabilities)
      infit_score_variance <- sum(infit_range^2 * probabilities) - infit_expected^2

      outfit_sum[[item_index]] <- outfit_sum[[item_index]] + score_weight * outfit_score
      outfit_mean[[item_index]] <- outfit_mean[[item_index]] + score_weight * outfit_expected
      outfit_var[[item_index]] <- outfit_var[[item_index]] + score_weight * outfit_score_variance

      infit_sum[[item_index]] <- infit_sum[[item_index]] + score_weight * infit_score
      infit_mean[[item_index]] <- infit_mean[[item_index]] + score_weight * infit_expected
      infit_var[[item_index]] <- infit_var[[item_index]] + score_weight * infit_score_variance
      infit_weight[[item_index]] <- infit_weight[[item_index]] + score_weight * variance
      n_used[[item_index]] <- n_used[[item_index]] + score_weight
    }
  }

  list(
    outfit_sum = outfit_sum,
    outfit_mean = outfit_mean,
    outfit_var = outfit_var,
    infit_sum = infit_sum,
    infit_mean = infit_mean,
    infit_var = infit_var,
    infit_weight = infit_weight,
    n_used = n_used
  )
}

#' Add source incomplete-record conditional item-fit contributions
#'
#' Source trace: `source/PAS_skunits/skbias15.pas::CalculateInAndOutfits`.
#' Mathematical step: reconstruct the available-item conditional distribution
#' for each compressed incomplete record and add its weighted moments.
#' @param accumulators Complete-record item-fit accumulators.
#' @param bundle Source-shaped analysis bundle.
#' @param item_gamma Fitted source item gamma matrix.
#' @param incomplete Compressed source incomplete-record representation.
#' @return Updated item-level outfit/infit accumulators.
#' @keywords internal
#' @noRd
conditional_item_fit_add_incomplete <- function(accumulators,
                                                bundle,
                                                item_gamma,
                                                incomplete) {
  items <- bundle$model$items
  for (record_index in seq_len(nrow(incomplete$records))) {
    record <- as.integer(incomplete$records[record_index, items$name, drop = TRUE])
    use_item <- record <= (items$raw_max - 1L)
    partial_score <- incomplete$score[[record_index]]
    count <- incomplete$count[[record_index]]
    for (item_index in which(use_item)) {
      item_score <- record[[item_index]]
      moments <- item_conditional_moment_for_subset(
        bundle,
        item_gamma,
        item_index,
        partial_score,
        use_item
      )
      variance <- moments$variance
      if (variance <= 0) {
        next
      }
      centered <- seq.int(0L, items$raw_max[[item_index]] - 1L) - moments$mean
      outfit_range <- centered^2 / variance
      infit_range <- centered^2

      accumulators$outfit_sum[[item_index]] <- accumulators$outfit_sum[[item_index]] +
        count * outfit_range[[item_score + 1L]]
      accumulators$outfit_mean[[item_index]] <- accumulators$outfit_mean[[item_index]] +
        count * sum(outfit_range * moments$probabilities)
      accumulators$outfit_var[[item_index]] <- accumulators$outfit_var[[item_index]] +
        count * (sum(outfit_range^2 * moments$probabilities) - sum(outfit_range * moments$probabilities)^2)

      accumulators$infit_sum[[item_index]] <- accumulators$infit_sum[[item_index]] +
        count * infit_range[[item_score + 1L]]
      accumulators$infit_mean[[item_index]] <- accumulators$infit_mean[[item_index]] +
        count * sum(infit_range * moments$probabilities)
      accumulators$infit_var[[item_index]] <- accumulators$infit_var[[item_index]] +
        count * (sum(infit_range^2 * moments$probabilities) - sum(infit_range * moments$probabilities)^2)
      accumulators$infit_weight[[item_index]] <- accumulators$infit_weight[[item_index]] +
        count * variance
      accumulators$n_used[[item_index]] <- accumulators$n_used[[item_index]] + count
    }
  }
  accumulators
}

#' Internal calculate item restscore gamma values helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle Source-shaped analysis bundle.
#' @param fit Fitted gRm model.
#' @param conditional Internal `conditional` value used by this helper.
#' @param incomplete Internal `incomplete` value used by this helper.
#' @return The internal `calculate_item_restscore_gamma_values()` computation result.
#' @keywords internal
#' @noRd
calculate_item_restscore_gamma_values <- function(bundle, fit, conditional = NULL, incomplete = NULL) {
  # Source trace: skbias14.pas::Calculate_item_restscore_gamma builds observed
  # and expected item-by-restscore crosstabs, calls CalculateGamma for observed
  # gamma, CalculateFittedGamma for expected gamma and ssgam, then uses
  # pIRgamma := 2*pnormal(abs((IRgamma-EIRgamma)/sdIRgamma)).
  items <- bundle$model$items
  data <- bundle$data
  n_items <- nrow(items)
  max_score <- bundle$model$max_total_score
  if (is.null(conditional)) {
    conditional <- item_conditional_moments(bundle, fit$item_gamma, include_probabilities = TRUE)
  }
  if (is.null(incomplete)) {
    incomplete <- empty_source_incomplete_records(bundle)
  }
  incomplete_gamma <- if (nrow(incomplete$records) > 0L) {
    prepare_incomplete_item_restscore_gamma(bundle, fit$item_gamma, incomplete)
  } else {
    NULL
  }
  item_matrix <- data[, items$name, drop = FALSE]
  complete_items <- apply(item_matrix >= 0L, 1L, all)
  item_scores <- rowSums(item_matrix * (item_matrix >= 0L))
  # Source trace: skbias15.pas::Count_Observed updates Nlowscore/Nhighscore
  # before GET_EXOGENE can reject rows with missing exogenous values. These
  # boundary cells therefore use complete-item rows, not only status==1 rows.
  low_count <- sum(complete_items & item_scores < bundle$model$least_score)
  high_count_source <- sum(complete_items & (item_scores > bundle$model$largest_score |
    item_scores == max_score))
  observed_rows <- complete_items
  if (nrow(bundle$model$backgrounds) > 0L) {
    observed_rows <- observed_rows & data$status == 1L
  }
  score_weights <- tabulate(
    item_scores[observed_rows] + 1L,
    nbins = max_score + 1L
  )

  observed_gamma <- expected_gamma <- gamma_sd <- p_gamma <- numeric(n_items)
  gamma_context <- list(
    items = items,
    data = data,
    max_score = max_score,
    conditional = conditional,
    incomplete = incomplete,
    incomplete_gamma = incomplete_gamma,
    item_scores = item_scores,
    observed_rows = observed_rows,
    score_weights = score_weights,
    low_count = low_count,
    high_count_source = high_count_source
  )
  for (item_index in seq_len(n_items)) {
    item_result <- item_restscore_gamma_item(item_index, gamma_context)
    observed_gamma[[item_index]] <- item_result$observed_gamma
    expected_gamma[[item_index]] <- item_result$expected_gamma
    gamma_sd[[item_index]] <- item_result$gamma_sd
    p_gamma[[item_index]] <- item_result$p_gamma
  }

  list(
    observed_gamma = observed_gamma,
    expected_gamma = expected_gamma,
    gamma_sd = gamma_sd,
    p_gamma = p_gamma
  )
}

#' Compute one source item-restscore gamma comparison
#'
#' Source trace: `source/PAS_skunits/skbias14.pas::Calculate_item_restscore_gamma`.
#' Mathematical step: build observed and fitted item-by-restscore tables with
#' deterministic extremes and compressed incomplete-record contributions, then
#' compare Goodman--Kruskal gamma using the source fitted variance.
#' @param item_index One-based item index.
#' @param gamma_context Prepared immutable inputs for the item loop.
#' @return Observed gamma, expected gamma, standard deviation, and p-value.
#' @keywords internal
#' @noRd
item_restscore_gamma_item <- function(item_index, gamma_context) {
  items <- gamma_context$items
  data <- gamma_context$data
  max_score <- gamma_context$max_score
  item_max <- items$raw_max[[item_index]] - 1L
  rest_max <- max_score - item_max
  observed <- matrix(0, nrow = item_max + 1L, ncol = rest_max + 1L)
  expected <- observed

  # Calculate_item_restscore_gamma initializes both tables with Nlowscore and
  # Nhighscore in their deterministic corner cells.
  observed[1L, 1L] <- gamma_context$low_count
  expected[1L, 1L] <- gamma_context$low_count
  observed[item_max + 1L, rest_max + 1L] <- gamma_context$high_count_source
  expected[item_max + 1L, rest_max + 1L] <- gamma_context$high_count_source

  for (total_score in seq.int(1L, max_score - 1L)) {
    n_score <- gamma_context$score_weights[[total_score + 1L]]
    if (n_score <= 0L) {
      next
    }
    moments <- gamma_context$conditional[[item_index]][[total_score + 1L]]
    probabilities <- moments$probabilities
    if (length(probabilities) == 0L) {
      next
    }
    score_rows <- gamma_context$observed_rows & gamma_context$item_scores == total_score
    observed_count <- tabulate(
      data[[items$name[[item_index]]]][score_rows] + 1L,
      nbins = item_max + 1L
    )
    observed_total <- sum(observed_count)
    observed_frequency <- if (observed_total > 0L) observed_count / observed_total else observed_count
    for (candidate_score in seq.int(0L, item_max)) {
      candidate_rest <- total_score - candidate_score
      if (candidate_rest >= 0L && candidate_rest <= rest_max) {
        # CalculateMeans stores relative item frequencies; the gamma routine
        # restores counts by multiplying observed and expected cells by n_score.
        observed[candidate_score + 1L, candidate_rest + 1L] <-
          observed[candidate_score + 1L, candidate_rest + 1L] +
          observed_frequency[[candidate_score + 1L]] * n_score
        expected[candidate_score + 1L, candidate_rest + 1L] <-
          expected[candidate_score + 1L, candidate_rest + 1L] +
          probabilities[[candidate_score + 1L]] * n_score
      }
    }
  }

  incomplete <- gamma_context$incomplete
  incomplete_gamma <- gamma_context$incomplete_gamma
  if (!is.null(incomplete_gamma)) {
    for (record_index in seq_len(nrow(incomplete$records))) {
      record <- as.integer(incomplete$records[record_index, items$name, drop = TRUE])
      item_score <- record[[item_index]]
      if (item_score > item_max) {
        next
      }
      count <- incomplete$count[[record_index]]
      observed_rest <- incomplete$score[[record_index]] - item_score
      if (observed_rest >= 0L && observed_rest <= rest_max) {
        observed[item_score + 1L, observed_rest + 1L] <-
          observed[item_score + 1L, observed_rest + 1L] + count
      }

      expected_score <- round(incomplete_gamma$expected_total[[record_index]])
      probabilities <- incomplete_gamma$probabilities[[record_index]][[item_index]]
      if (length(probabilities) == 0L) {
        next
      }
      for (candidate_score in seq.int(0L, item_max)) {
        candidate_rest <- expected_score - candidate_score
        if (candidate_rest >= 0L && candidate_rest <= rest_max) {
          expected[candidate_score + 1L, candidate_rest + 1L] <-
            expected[candidate_score + 1L, candidate_rest + 1L] +
            probabilities[[candidate_score + 1L]] * count
        }
      }
    }
  }

  observed_gamma <- goodman_kruskal_gamma(observed)
  fitted <- fitted_gamma_stats(expected)
  gamma_sd <- sqrt(fitted$variance)
  p_gamma <- if (gamma_sd > 0) {
    two_sided_source_normal_p((observed_gamma - fitted$gamma) / gamma_sd)
  } else {
    1
  }
  list(
    observed_gamma = observed_gamma,
    expected_gamma = fitted$gamma,
    gamma_sd = gamma_sd,
    p_gamma = p_gamma
  )
}
