#' Calculate a PCM threshold from adjacent gamma values
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param gamma_values Numeric vector of item gamma values indexed by score plus
#'   one.
#' @param score Positive score category whose threshold is calculated against the
#'   previous category.
#' @return The log adjacent-category threshold, or the Pascal sentinel value for
#'   structural zero cases.
#' @keywords internal
#' @noRd
source_threshold_from_gamma <- function(gamma_values, score) {
  # Source trace:
  # pascal_harness/SourceRaschCore.pas::SourceItemThresholdFromGamma.
  threshold_epsilon <- 0.000000000001
  if (gamma_values[[score + 1L]] <= threshold_epsilon) {
    9999999
  } else if (gamma_values[[score]] <= threshold_epsilon) {
    -9999999
  } else {
    # Adjacent-category PCM threshold: tau_score = log(gamma[score-1] / gamma[score]).
    log(gamma_values[[score]] / gamma_values[[score + 1L]])
  }
}

#' Calculate the source item location from gamma values
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param gamma_values Numeric vector of item gamma values.
#' @param max_score Maximum score for the item.
#' @return Mean finite adjacent-category threshold, or zero if no finite
#'   threshold is available.
#' @keywords internal
#' @noRd
source_location_from_gamma <- function(gamma_values, max_score) {
  # Source trace:
  # pascal_harness/SourceRaschCore.pas::SourceItemLocationFromGamma.
  total <- 0
  n_finite <- 0L
  for (score in seq.int(1L, max_score)) {
    threshold <- source_threshold_from_gamma(gamma_values, score)
    if (abs(threshold) < 9999999) {
      # DIGRAM averages only finite adjacent-category thresholds.
      total <- total + threshold
      n_finite <- n_finite + 1L
    }
  }
  if (n_finite > 0L) {
    total / n_finite
  } else {
    0
  }
}

#' Source-compatible log with zero sentinel
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param value Numeric scalar.
#' @return `log(value)`, or `-9999999` when `value` is zero.
#' @keywords internal
#' @noRd
ln_zero <- function(value) {
  # Source trace: pascal_harness/SourceRaschCore.pas::LnZero.
  if (value == 0) {
    -9999999
  } else {
    log(value)
  }
}

#' Calculate the expected item score at a theta value
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param theta Positive person parameter on the multiplicative scale.
#' @param max_score Maximum score for the item.
#' @param gamma_values Numeric vector of item gamma values.
#' @return Expected item score.
#' @keywords internal
#' @noRd
true_score_from_gamma <- function(theta, max_score, gamma_values) {
  # Source trace: pascal_harness/SourceRaschCore.pas::CalculateTrueScoreFromGamma.
  if (theta <= 0) {
    return(0)
  }
  numerator <- 0
  denominator <- 0
  power_value <- 1
  for (score in seq.int(0L, max_score)) {
    if (score == 0L) {
      power_value <- 1
    } else {
      power_value <- power_value * theta
    }
    # E(X | theta) = sum_x x * theta^x * gamma[x] /
    #                sum_x     theta^x * gamma[x].
    numerator <- numerator + score * power_value * gamma_values[[score + 1L]]
    denominator <- denominator + power_value * gamma_values[[score + 1L]]
  }
  if (denominator > 0) {
    numerator / denominator
  } else {
    0
  }
}

#' Initialize a person parameter search
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param target_score Target expected score.
#' @param max_score Maximum score for the item.
#' @param gamma_values Numeric vector of item gamma values.
#' @return Initial theta value on the multiplicative scale.
#' @keywords internal
#' @noRd
initialize_person_parameter2 <- function(target_score, max_score, gamma_values) {
  # Source trace:
  # pascal_harness/SourceRaschCore.pas::InitializePersonParameter2.
  max_diff <- 10 * max_score
  theta <- 1
  # The source uses a coarse exp(-5)..exp(5) grid to seed Newton iteration.
  for (index in seq.int(0L, 100L)) {
    trial_theta <- exp(-5 + index / 10)
    true_score <- true_score_from_gamma(trial_theta, max_score, gamma_values)
    delta <- abs(true_score - target_score)
    if (delta < max_diff) {
      theta <- trial_theta
      max_diff <- delta
    }
  }
  theta
}

#' Estimate a person parameter for a target score
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param target_score Target expected score.
#' @param max_score Maximum score for the item.
#' @param gamma_values Numeric vector of item gamma values.
#' @param max_iterations Maximum Newton iterations.
#' @param max_delta Lower bound used when theta becomes non-positive.
#' @return Estimated theta on the multiplicative scale, using source sentinel
#'   behavior for boundary scores.
#' @keywords internal
#' @noRd
estimate_person_parameter <- function(target_score, max_score, gamma_values, max_iterations = 1000L, max_delta = 0.0001) {
  # Source trace: pascal_harness/SourceRaschCore.pas::EstimatePersonParameter.
  theta <- initialize_person_parameter2(target_score, max_score, gamma_values)
  if (target_score == 0) {
    return(0)
  }
  if (target_score == max_score) {
    return(999999)
  }
  if (theta <= 0) {
    theta <- max_delta
  }

  for (iter in seq_len(max_iterations)) {
    h0 <- 0
    h1 <- 0
    h2 <- 0
    factor <- 1 / theta
    for (score in seq.int(0L, max_score)) {
      factor <- factor * theta
      product0 <- gamma_values[[score + 1L]] * factor
      product1 <- product0 * score
      product2 <- product1 * score
      # H0, H1, and H2 are the normalizing sum and first two score moments
      # on the multiplicative theta scale.
      h0 <- h0 + product0
      h1 <- h1 + product1
      h2 <- h2 + product2
    }
    # Newton step copied from the Pascal source for solving E(X | theta)=target.
    step <- -theta * (target_score * h0 - h1)
    if (h0 > 0) {
      if (h1 == h0) {
        ratio <- 1
      } else {
        ratio <- h1 / h0
      }
      denominator <- ratio * h1 - h2
      if (denominator != 0) {
        step <- step / denominator
      }
      theta <- theta + step
    }
    if (theta <= 0) {
      theta <- max_delta
    }
    if (abs(step) < 0.001 || theta <= 0 || theta > 999999) {
      break
    }
  }
  theta
}

#' Calculate the item midpoint difficulty
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param gamma_values Numeric vector of item gamma values.
#' @param max_score Maximum score for the item.
#' @return Source midpoint difficulty on the log scale.
#' @keywords internal
#' @noRd
source_difficulty_from_gamma <- function(gamma_values, max_score) {
  # Source trace:
  # pascal_harness/SourceRaschCore.pas::SourceItemDifficultyFromGamma.
  if (max_score == 1L) {
    return(-ln_zero(gamma_values[[2L]]))
  }
  # For polytomous items DIGRAM defines the midpoint as the log theta where
  # expected item score equals half the maximum item score.
  theta <- estimate_person_parameter(max_score / 2, max_score, gamma_values, 1000L)
  ln_zero(theta)
}

#' Calculate score probabilities at a log-theta value
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param log_theta Person parameter on the log scale.
#' @param max_score Maximum score for the item.
#' @param gamma_values Numeric vector of item gamma values.
#' @return A list containing score probabilities, mean score, and score
#'   variance.
#' @keywords internal
#' @noRd
score_probabilities_for_log_theta <- function(log_theta, max_score, gamma_values) {
  # Source trace:
  # pascal_harness/SourceRaschCore.pas::ScoreProbabilitiesForLogTheta.
  theta <- exp(log_theta)
  power_value <- 1
  probabilities <- numeric(max_score + 1L)
  total <- 0
  for (score in seq.int(0L, max_score)) {
    if (score == 0L) {
      power_value <- 1
    } else {
      power_value <- power_value * theta
    }
    # Unnormalized PCM score mass at log-theta: exp(log_theta * x) * gamma[x].
    probabilities[[score + 1L]] <- power_value * gamma_values[[score + 1L]]
    total <- total + probabilities[[score + 1L]]
  }

  mean_score <- 0
  sum_sq <- 0
  if (total > 0) {
    for (score in seq.int(0L, max_score)) {
      probabilities[[score + 1L]] <- probabilities[[score + 1L]] / total
      mean_score <- mean_score + score * probabilities[[score + 1L]]
      sum_sq <- sum_sq + score * score * probabilities[[score + 1L]]
    }
  }

  # In the source this variance is also the item information curve being
  # maximized by CalculateSourceItemTarget.
  list(
    probabilities = probabilities,
    mean_score = mean_score,
    score_variance = sum_sq - mean_score * mean_score
  )
}

#' Find the item information target
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param gamma_values Numeric vector of item gamma values.
#' @param max_score Maximum score for the item.
#' @return A list containing the log-theta target and information value where
#'   item information is maximized.
#' @keywords internal
#' @noRd
source_item_target <- function(gamma_values, max_score) {
  # Source trace:
  # pascal_harness/SourceRaschCore.pas::CalculateSourceItemTarget.
  modes <- source_item_target_modes(gamma_values, max_score)
  list(target = modes$target, info = modes$info)
}

#' Find all source item-information modes
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param gamma_values Numeric vector of item gamma values.
#' @param max_score Maximum score for the item.
#' @return A list with the global target and all coarse local modes.
#' @keywords internal
#' @noRd
source_item_target_modes <- function(gamma_values, max_score) {
  # Source trace:
  # pascal_harness/SourceRaschCore.pas::CalculateSourceItemTarget1Modes.
  target <- 0
  max_inf <- 0
  super_target <- target
  inf_before <- 0
  going_up <- TRUE
  mode_targets <- numeric()
  mode_info <- numeric()

  for (ii in seq.int(-10000L, 10000L)) {
    theta <- ii * 0.001
    props <- score_probabilities_for_log_theta(theta, max_score, gamma_values)
    if (going_up && props$score_variance < inf_before) {
      mode_targets <- c(mode_targets, theta - 0.001)
      mode_info <- c(mode_info, inf_before)
      going_up <- FALSE
    } else if (!going_up && props$score_variance > inf_before) {
      going_up <- TRUE
    }
    inf_before <- props$score_variance
    if (props$score_variance > max_inf) {
      max_inf <- props$score_variance
      super_target <- theta
    }
  }

  target <- super_target
  d <- 0.001
  for (ii in seq.int(-25L, 25L)) {
    theta <- target + ii * d / 25
    props <- score_probabilities_for_log_theta(theta, max_score, gamma_values)
    if (props$score_variance > max_inf) {
      max_inf <- props$score_variance
      super_target <- theta
    }
  }

  list(
    target = super_target,
    info = max_inf,
    modes = data.frame(target = mode_targets, info = mode_info)
  )
}

#' Calculate item information at a log-theta value
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param gamma_values Numeric vector of item gamma values.
#' @param max_score Maximum score for the item.
#' @param log_theta Person parameter on the log scale.
#' @return Item information at `log_theta`.
#' @keywords internal
#' @noRd
source_item_information_from_gamma <- function(gamma_values, max_score, log_theta) {
  # Source trace: pascal_harness/SourceRaschCore.pas::SourceItemInformation.
  fx <- numeric(max_score + 1L)
  f1x <- numeric(max_score + 1L)
  f2x <- numeric(max_score + 1L)
  g <- 0
  g1 <- 0
  g2 <- 0
  for (score in seq.int(0L, max_score)) {
    fx[[score + 1L]] <- exp(log_theta * score) * gamma_values[[score + 1L]]
    g <- g + fx[[score + 1L]]
    f1x[[score + 1L]] <- score * fx[[score + 1L]]
    g1 <- g1 + f1x[[score + 1L]]
    f2x[[score + 1L]] <- score * f1x[[score + 1L]]
    g2 <- g2 + f2x[[score + 1L]]
  }
  if (g <= 0) {
    return(0)
  }

  g_sqr <- g * g
  g_cube <- g_sqr * g
  result <- 0
  for (score in seq.int(0L, max_score)) {
    # PX is the normalized category probability; P1X and P2X are the first and
    # second derivatives of that probability with respect to log theta.
    px <- fx[[score + 1L]] / g
    p1x <- f1x[[score + 1L]] / g - fx[[score + 1L]] * g1 / g_sqr
    p2x <- f2x[[score + 1L]] / g - 2 * f1x[[score + 1L]] * g1 / g_sqr -
      fx[[score + 1L]] * g2 / g_cube + 2 * fx[[score + 1L]] * g1 * g1 / g_cube
    if (px > 0) {
      result <- result + p1x * p1x / px - p2x
    }
  }
  result
}

#' Preserve the source signed-zero display for top ICE cancellation
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param value Numeric ICE value for the highest item category.
#' @return Numeric value with source-shaped signed zero when it formats to zero.
#' @keywords internal
#' @noRd
source_top_ice_cancellation_zero <- function(value) {
  # Source trace: source/digram_source_20260817/skunits/skbias22.pas::write_iteminformation1 writes the highest
  # ICE category as ln(x[i,k]) - k * (ln(x[i,k]) / k). The value is
  # mathematically zero, but Pascal Extended arithmetic preserves a display
  # sign. R's double calculation often collapses that cancellation to +0, so
  # preserve the source-shaped signed zero for values that print as zero.
  if (is.na(value) || abs(value) >= 0.0005) {
    return(value)
  }
  if (value < 0) {
    -0
  } else {
    0
  }
}

#' Calculate the observed possible score range
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param item_counts Item-by-score observed count matrix.
#' @return Integer vector of length two containing the summed observed minimum
#'   and maximum item scores.
#' @keywords internal
#' @noRd
calculate_observed_score_range <- function(item_counts) {
  # Source trace:
  # pascal_harness/item_parameters_report/BIRT_ITEM_PARAMETERS_REPORT.pas::
  # CalculateObservedScoreRange.
  score_min <- 0L
  score_max <- 0L
  for (item_index in seq_len(nrow(item_counts))) {
    positive <- which(item_counts[item_index, ] > 0L)
    if (length(positive) > 0L) {
      scores <- as.integer(colnames(item_counts)[positive])
      # DIGRAM sums each item's lowest and highest observed supported category,
      # rather than using the empirical minimum and maximum total scores.
      score_min <- score_min + min(scores)
      score_max <- score_max + max(scores)
    }
  }
  c(score_min, score_max)
}

#' Count source-style estimated parameters
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param item_counts Item-by-score observed count matrix.
#' @return Integer number of estimated parameters using the source report
#'   convention.
#' @keywords internal
#' @noRd
calculate_source_n_parameters <- function(item_counts) {
  # Source trace:
  # pascal_harness/item_parameters_report/BIRT_ITEM_PARAMETERS_REPORT.pas::
  # CalculateSourceNParameters.
  result <- -1L
  for (item_index in seq_len(nrow(item_counts))) {
    # Source convention: each item contributes positive-category-count minus
    # one, and the model as a whole subtracts one additional reference degree.
    result <- result + sum(item_counts[item_index, ] > 0L) - 1L
  }
  result
}

#' Report native Pascal Extended compatibility
#'
#' The package can reproduce Pascal Extended fixed-field cancellation only
#' when the compiled platform provides a 64-bit significand, the Extended
#' exponent range, and the x87 control/logarithm path used by this native unit.
#' Other platforms still return ordinary platform-long-double formatting, but
#' this capability flag is false and exact Extended parity is not claimed.
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @return A schema-versioned list containing long-double mantissa bits,
#'   maximum exponent, storage size, x87 feature flags, and the final
#'   `pascal_extended_fixed_field_supported` flag.
#' @keywords internal
source_extended_native_capability <- function() {
  .Call("gRm_item_parameters_extended_capabilities", PACKAGE = "gRm")
}

#' Replay one source Extended top-category ICE field
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param gamma_top Numeric highest-category item gamma.
#' @param max_score Integer highest item score.
#' @return Character scalar formatted as a nine-character source field.
#' @keywords internal
#' @noRd
source_extended_top_ice_cancellation_field <- function(gamma_top, max_score) {
  .Call(
    "gRm_item_parameters_top_ice_field",
    as.numeric(gamma_top),
    as.integer(max_score),
    PACKAGE = "gRm"
  )
}

#' Replay source Extended item-parameter ICE fields
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param fit A fitted model list returned by `fit_rasch_base()`.
#' @param bundle The item-parameters bundle used to produce `fit`.
#' @return Character matrix of fixed-width ICE fields from the native
#'   `skbias12.pas` `CalculateICEandMICE` arithmetic path.
#' @keywords internal
#' @noRd
source_item_parameters_extended_ice_fields <- function(fit, bundle) {
  .Call(
    "gRm_item_parameters_ice_fields_from_gamma",
    fit$item_gamma,
    as.integer(bundle$model$items$raw_max),
    PACKAGE = "gRm"
  )
}

#' Calculate source-shaped item-parameters input statistics
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param bundle Source-shaped source item-parameters bundle.
#' @return A list with the input counters printed by the DIGRAM report.
#' @param gllrm_context Internal `gllrm_context` value used by this helper.
#' @keywords internal
#' @noRd
item_parameters_input_stats <- function(bundle, gllrm_context = NULL) {
  # Source trace:
  # source/digram_source_20260817/skunits/skbias12b.pas::Count_Margins, with report counters
  # printed by the WriteParameters branch for complete item responses.
  data <- bundle$data
  items <- bundle$model$items
  backgrounds <- bundle$model$backgrounds
  max_total_score <- sum(items$raw_max - 1L)
  source_case_columns <- c(items$name, backgrounds$name)
  source_case_order <- if (length(source_case_columns) > 0L) {
    do.call(order, as.data.frame(data[, source_case_columns, drop = FALSE]))
  } else {
    seq_len(nrow(data))
  }
  ld_components <- gllrm_context$ld_components_items %||% list()

  n_incomplete <- 0L
  n_responses <- 0L
  n_useful <- 0L
  n_useless <- 0L
  for (row_index in source_case_order) {
    item_values <- unlist(data[row_index, items$name, drop = FALSE], use.names = FALSE)
    background_values <- if (nrow(backgrounds) > 0L) {
      unlist(data[row_index, backgrounds$name, drop = FALSE], use.names = FALSE)
    } else {
      integer()
    }
    complete_items <- all(item_values >= 0L)
    complete_backgrounds <- all(background_values >= 0L)
    score <- sum(item_values[item_values >= 0L])

    if (complete_items) {
      # Source trace: Count_Margins in source/digram_source_20260817/skunits/skbias12b.pas
      # computes and filters the complete item score before GET_EXOGENE can
      # mark missing/invalid background values as Nuseless.
      if (score < bundle$model$least_score || score > max_total_score - 1L) {
        next
      }
      if (!complete_backgrounds) {
        n_useless <- n_useless + 1L
      }
      next
    } else if (!complete_items && complete_backgrounds) {
      n_incomplete <- n_incomplete + 1L
      pure_items <- item_values
      # Count_Margins.DealWithMissingInformation removes every response in an
      # LD component when any member of that component is missing. The
      # remaining one-based Pascal responses correspond here to nonnegative
      # zero-based values; a zero score is still a response.
      for (component in ld_components) {
        component <- as.integer(component)
        if (length(component) && any(pure_items[component] < 0L)) {
          pure_items[component] <- -1L
        }
      }
      retained <- pure_items >= 0L
      responses <- sum(retained)
      pure_score <- sum(pure_items[retained])
      pure_max_score <- sum(items$raw_max[retained] - 1L)
      # DIGRAM stores Nresponses from the latest useful incomplete TABDATA cell
      # in source lexicographic variable order rather than summing responses.
      if (
        responses > 1L &&
          pure_score > 0L &&
          pure_score < pure_max_score
      ) {
        n_useful <- n_useful + 1L
        n_responses <- responses
      }
    } else {
      n_useless <- n_useless + 1L
    }
  }

  list(
    valid_score_min = bundle$model$least_score,
    valid_score_max = max_total_score - 1L,
    n_valid = sum(data$status == 1L),
    n_incomplete = n_incomplete,
    n_responses = n_responses,
    n_useful = n_useful,
    n_useless = n_useless
  )
}

#' Derive item-parameters report values
#'
#' Calculates the derived values rendered in the DIGRAM `ItemParameters.txt`
#' report from a fitted base Rasch model. This includes PCM thresholds and
#' locations, ICE/MICE effects, item effect terms, midpoint difficulties,
#' information targets, observed score range, and source-style parameter counts.
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::GLLRM_output`.
#' @param fit A fitted model list returned by `fit_rasch_base()`.
#' @param bundle The item-parameters bundle used to produce `fit`.
#' @return A list containing report-ready item names, labels, fitted gamma
#'   parameters, threshold/location matrices, ICE/MICE effects, item statistics,
#'   convergence diagnostics, observed score range, and parameter count.
#' @examples
#' \dontrun{
#' project <- read_digram_project("path/to/DIGRAM")
#' bundle <- build_item_parameters_bundle(project)
#' fit <- fit_rasch_base(bundle)
#' values <- item_parameters_values(fit, bundle)
#' values$item_statistics
#' }
#' @keywords internal
#' @noRd
item_parameters_values <- function(fit, bundle) {
  # Source trace:
  # pascal_harness/SourceRaschCore.pas::EmitGLLRMOutputRows output(4) item
  # parameter/PCM/statistics rows, plus BIRT_ITEM_PARAMETERS_REPORT.pas footer
  # helpers for observed score range and parameter count.
  items <- bundle$model$items
  max_category_count <- ncol(fit$item_gamma)
  score_names <- as.character(seq.int(0L, max_category_count - 1L))
  item_names <- items$name

  thresholds <- matrix(
    NA_real_,
    nrow = nrow(items),
    ncol = max_category_count - 1L,
    dimnames = list(item_names, as.character(seq.int(1L, max_category_count - 1L)))
  )
  locations <- stats::setNames(numeric(nrow(items)), item_names)
  midpoints <- locations
  targets <- locations
  info_at_target <- locations
  info_per_step <- locations
  item_effect <- locations
  mice_item_effect <- locations
  ice <- matrix(NA_real_, nrow = nrow(items), ncol = max_category_count, dimnames = list(item_names, score_names))
  mice <- ice

  for (item_index in seq_len(nrow(items))) {
    max_score <- items$raw_max[[item_index]] - 1L
    gamma_values <- as.numeric(fit$item_gamma[item_index, seq_len(items$raw_max[[item_index]])])
    # PCM rows are derived directly from fitted item gammas.
    for (score in seq.int(1L, max_score)) {
      thresholds[item_index, as.character(score)] <- source_threshold_from_gamma(gamma_values, score)
    }
    locations[[item_index]] <- source_location_from_gamma(gamma_values, max_score)
    midpoints[[item_index]] <- source_difficulty_from_gamma(gamma_values, max_score)
    target <- source_item_target(gamma_values, max_score)
    targets[[item_index]] <- target$target
    info_at_target[[item_index]] <- target$info
    info_per_step[[item_index]] <- if (max_score >= 2L) target$info / max_score else NA_real_

    if (max_score > 0L && gamma_values[[max_score + 1L]] > 0) {
      # The report ICE/MICE rows use the output(4) block:
      # z = log(gamma[max_score]) / max_score,
      # ICE[x] = log(gamma[x]) - x*z, MICE[x] = exp(ICE[x]).
      # This is distinct from the lower-level CalculateSourceICEAndMICEFromGamma
      # helper, which centers by PCM location for other source paths.
      z <- log(gamma_values[[max_score + 1L]]) / max_score
      for (score in seq.int(0L, max_score)) {
        if (gamma_values[[score + 1L]] > 0) {
          ice[item_index, as.character(score)] <- log(gamma_values[[score + 1L]]) - score * z
          if (score == max_score) {
            ice[item_index, as.character(score)] <-
              source_top_ice_cancellation_zero(ice[item_index, as.character(score)])
          }
          mice[item_index, as.character(score)] <- exp(ice[item_index, as.character(score)])
        } else {
          ice[item_index, as.character(score)] <- 0
          mice[item_index, as.character(score)] <- 0
        }
      }
      item_effect[[item_index]] <- z
      mice_item_effect[[item_index]] <- exp(z)
    }
  }

  extended_capability <- source_extended_native_capability()
  result <- list(
    n_step = fit$n_step,
    delta = fit$report_delta,
    log_likelihood = base_rasch_loglike(bundle, fit$item_gamma),
    likelihood_n = fit$counts$n_valid,
    input_stats = item_parameters_input_stats(bundle),
    item_names = item_names,
    item_labels = items$label_code,
    background_labels = bundle$model$backgrounds$label_code,
    background_max = bundle$model$backgrounds$raw_max,
    n_parameters = calculate_source_n_parameters(fit$counts$item_counts),
    observed_score_range = calculate_observed_score_range(fit$counts$item_counts),
    item_gamma = fit$item_gamma,
    ice_fields = source_item_parameters_extended_ice_fields(fit, bundle),
    extended_capability = extended_capability,
    thresholds = thresholds,
    locations = locations,
    ice = ice,
    mice = mice,
    ice_item_effect = item_effect,
    mice_item_effect = mice_item_effect,
    item_statistics = data.frame(
      item = item_names,
      location = unname(locations),
      midpoint = unname(midpoints),
      target = unname(targets),
      info_at_target = unname(info_at_target),
      info_per_step = unname(info_per_step),
      stringsAsFactors = FALSE,
      row.names = item_names
    )
  )
  class(result) <- c("gRm_item_parameters_values", class(result))
  result
}
