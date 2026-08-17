#' Internal calculate extended item fit values helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle Source-shaped analysis bundle.
#' @param fit Fitted gRm model.
#' @param conditional Internal `conditional` value used by this helper.
#' @return The internal `calculate_extended_item_fit_values()` computation result.
#' @keywords internal
#' @noRd
calculate_extended_item_fit_values <- function(bundle, fit, conditional = NULL) {
  # Source trace: skbias15.pas::Calculate_residuals_and_item_fits, nested
  # CalculateOutfit/CalculateInfit and the following Outfit/Infit summaries.
  # This is the diagnostic extended surface; compact finalization remains in
  # calculate_conditional_item_fit_values().
  items <- bundle$model$items
  data <- bundle$data
  n_items <- nrow(items)
  if (is.null(conditional)) {
    conditional <- item_conditional_moments(bundle, fit$item_gamma, include_probabilities = TRUE)
  }

  item_matrix <- data[, items$name, drop = FALSE]
  complete_items <- apply(item_matrix >= 0L, 1L, all)
  item_scores <- rowSums(item_matrix * (item_matrix >= 0L))
  score_weights <- tabulate(
    item_scores[complete_items] + 1L,
    nbins = bundle$model$max_total_score + 1L
  )
  source_score_counts <- fit$counts$score_counts
  score_names <- as.character(seq.int(0L, bundle$model$max_total_score))
  source_score_counts <- stats::setNames(
    as.numeric(source_score_counts[score_names]),
    score_names
  )
  source_score_counts[is.na(source_score_counts)] <- 0
  # Source trace: skbias15.pas::Count_Observed updates Nlowscore/Nhighscore
  # before rejecting rows on exogenous range, while Add_count_to_tables fills
  # the interior ItemMargTables only after exogenous validation. The extended
  # item-restscore gamma tables therefore use complete-item endpoint counts and
  # valid-background interior score counts.
  source_score_counts[[1L]] <- score_weights[[1L]]
  source_score_counts[[length(source_score_counts)]] <- score_weights[[length(score_weights)]]
  observed_rows <- complete_items
  if (nrow(bundle$model$backgrounds) > 0L) {
    observed_rows <- observed_rows & data$status == 1L
  }
  source_tables <- extended_item_restscore_source_tables(
    bundle,
    conditional,
    source_score_counts,
    complete_items,
    observed_rows,
    item_scores,
    score_weights
  )

  score_rows <- list()
  summary_rows <- list()
  for (item_index in seq_len(n_items)) {
    item_score_rows <- list()

    for (score in seq.int(1L, bundle$model$max_total_score - 1L)) {
      n_score <- score_weights[[score + 1L]]
      row <- extended_item_fit_score_row(
        items = items,
        data = data,
        item_index = item_index,
        score = score,
        n_score = n_score,
        observed_rows = observed_rows,
        item_scores = item_scores,
        conditional = conditional
      )
      if (is.null(row)) {
        next
      }
      score_rows[[length(score_rows) + 1L]] <- row
      item_score_rows[[length(item_score_rows) + 1L]] <- row
    }

    item_scores_df <- do.call(rbind_fill, item_score_rows)
    summary_rows[[length(summary_rows) + 1L]] <- extended_item_fit_summary_row(
      items,
      item_index,
      item_scores_df
    )
  }

  list(
    scores = type.convert(do.call(rbind_fill, score_rows), as.is = TRUE),
    summaries = type.convert(do.call(rbind_fill, summary_rows), as.is = TRUE),
    score_n = data.frame(
      score = seq.int(0L, bundle$model$max_total_score),
      n = as.numeric(score_weights),
      stringsAsFactors = FALSE
    ),
    restscore_tables = source_tables$restscore_tables,
    local_restscore_tables = source_tables$local_restscore_tables,
    local_gamma = source_tables$local_gamma
  )
}

#' Compute one extended item-fit score row
#'
#' Source trace: `source/PAS_skunits/skbias15.pas::CalculateOutfit` and
#' `source/PAS_skunits/skbias15.pas::CalculateInfit`.
#' Mathematical step: form observed/expected conditional moments at one total
#' score and preserve the source's positive-residual and score-count weights.
#' @param items Parsed item metadata.
#' @param data Encoded source-shaped analysis data.
#' @param item_index One-based item index.
#' @param score Zero-based total score.
#' @param n_score Complete-item score frequency.
#' @param observed_rows Valid observed-margin record mask.
#' @param item_scores Zero-based total item score per record.
#' @param conditional Conditional item moments by item and score.
#' @return One extended score row, or `NULL` for an unused score.
#' @keywords internal
#' @noRd
extended_item_fit_score_row <- function(items,
                                        data,
                                        item_index,
                                        score,
                                        n_score,
                                        observed_rows,
                                        item_scores,
                                        conditional) {
  item_max <- items$raw_max[[item_index]] - 1L
  item_values <- seq.int(0L, item_max)
  score_rows_mask <- observed_rows & item_scores == score
  observed_count <- tabulate(
    data[[items$name[[item_index]]]][score_rows_mask] + 1L,
    nbins = item_max + 1L
  )
  observed_total <- sum(observed_count)
  if (n_score <= 0L || observed_total <= 0L) {
    return(NULL)
  }
  moments <- conditional[[item_index]][[score + 1L]]
  variance <- moments$variance
  if (variance <= 0) {
    return(NULL)
  }

  probabilities <- moments$probabilities
  observed_frequency <- observed_count / observed_total
  residual <- item_values - moments$mean
  squared_residual <- residual^2
  standardized <- residual / sqrt(variance)
  squared_standardized <- squared_residual / variance

  # The source report rounds this value for display only. Its computational
  # row retains the mixed valid-background frequency times the broader
  # complete-item score count.
  observed_source_value <- observed_frequency * n_score
  positive <- residual > 0
  positive_expected_squared <- sum(squared_residual[positive] * probabilities[positive])
  positive_expected_variance <-
    sum(squared_residual[positive]^2 * probabilities[positive]) -
    positive_expected_squared^2
  positive_observed_squared <- sum(squared_residual[positive] * observed_frequency[positive])
  positive_observed_standardized <- sum(standardized[positive] * observed_frequency[positive])

  outfit_contribution <- sum(squared_standardized * observed_frequency)
  outfit_expected <- sum(squared_standardized * probabilities)
  outfit_variance <- sum(squared_standardized^2 * probabilities) - outfit_expected^2
  outfit_standard_error <- sqrt(max(outfit_variance, 0) / n_score)
  outfit_z <- if (outfit_standard_error > 0) {
    (outfit_contribution - 1) / outfit_standard_error
  } else {
    0
  }

  infit_average <- sum(squared_residual * observed_frequency)
  infit_expected <- sum(squared_residual * probabilities)
  infit_variance <- sum(squared_residual^2 * probabilities) - infit_expected^2
  row <- data.frame(
    item_label = items$label_code[[item_index]],
    item_name = items$name[[item_index]],
    score = score,
    n = n_score,
    observed_mean = sum(item_values * observed_frequency),
    observed_variance = sum(item_values^2 * observed_frequency) -
      sum(item_values * observed_frequency)^2,
    expected_mean = moments$mean,
    expected_variance = variance,
    residual_mean = 0,
    residual_variance = variance,
    n_variance = n_score * variance,
    squared_residual_mean = positive_expected_squared,
    squared_residual_variance = positive_expected_variance,
    squared_residual_observed_average = positive_observed_squared,
    standardized_mean = 0,
    standardized_variance = 1,
    standardized_observed_average = positive_observed_standardized,
    squared_standardized_mean = outfit_expected,
    squared_standardized_variance = outfit_variance,
    squared_standardized_observed_average = safe_ratio(positive_observed_squared, variance),
    outfit_contribution = outfit_contribution,
    outfit_standard_error = outfit_standard_error,
    outfit_z = outfit_z,
    outfit_p = two_sided_source_normal_p(outfit_z),
    infit_average = infit_average,
    infit_expected = infit_expected,
    infit_ratio = safe_ratio(infit_average, infit_expected),
    infit_variance = infit_variance,
    stringsAsFactors = FALSE
  )
  for (score_value in item_values) {
    suffix <- as.character(score_value)
    row[[paste0("observed_", suffix)]] <- observed_source_value[[score_value + 1L]]
    row[[paste0("probability_", suffix)]] <- probabilities[[score_value + 1L]]
    row[[paste0("residual_", suffix)]] <- residual[[score_value + 1L]]
    row[[paste0("squared_residual_", suffix)]] <- squared_residual[[score_value + 1L]]
    row[[paste0("standardized_", suffix)]] <- standardized[[score_value + 1L]]
    row[[paste0("squared_standardized_", suffix)]] <- squared_standardized[[score_value + 1L]]
  }
  row
}

#' Summarize extended score rows for one item
#'
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param items Parsed item metadata.
#' @param item_index One-based item index.
#' @param item_scores_df Extended rows for the item.
#' @return One source-weighted extended item summary row.
#' @keywords internal
#' @noRd
extended_item_fit_summary_row <- function(items, item_index, item_scores_df) {
  summary_keep <- item_scores_df$n > 1
  row <- data.frame(
    item_label = items$label_code[[item_index]],
    item_name = items$name[[item_index]],
    outfit_total_n = sum(item_scores_df$n[summary_keep]),
    outfit_total_observed = safe_ratio(
      sum(item_scores_df$n[summary_keep] * item_scores_df$squared_residual_observed_average[summary_keep]),
      sum(item_scores_df$n[summary_keep])
    ),
    outfit_total_expected = safe_ratio(
      sum(item_scores_df$n[summary_keep] * item_scores_df$squared_residual_mean[summary_keep]),
      sum(item_scores_df$n[summary_keep])
    ),
    outfit_total_value = safe_ratio(
      sum(item_scores_df$n[summary_keep] *
        safe_ratio(item_scores_df$squared_residual_observed_average[summary_keep],
          item_scores_df$squared_residual_mean[summary_keep])),
      sum(item_scores_df$n[summary_keep])
    ),
    infit_observed = sum(item_scores_df$n * item_scores_df$infit_average),
    infit_expected = sum(item_scores_df$n * item_scores_df$infit_expected),
    infit_variance = sum(item_scores_df$n * item_scores_df$infit_variance),
    stringsAsFactors = FALSE
  )
  row$infit_value <- safe_ratio(row$infit_observed, row$infit_expected)
  row
}

#' Internal extended item restscore source tables helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle Source-shaped analysis bundle.
#' @param conditional Internal `conditional` value used by this helper.
#' @param base_score_counts Internal `base_score_counts` value used by this helper.
#' @param complete_items Internal `complete_items` value used by this helper.
#' @param observed_rows Internal `observed_rows` value used by this helper.
#' @param item_scores Internal `item_scores` value used by this helper.
#' @param score_weights Internal `score_weights` value used by this helper.
#' @return The internal `extended_item_restscore_source_tables()` computation result.
#' @keywords internal
#' @noRd
extended_item_restscore_source_tables <- function(
    bundle,
    conditional,
    base_score_counts,
    complete_items,
    observed_rows,
    item_scores,
    score_weights) {
  items <- bundle$model$items
  data <- bundle$data
  n_items <- nrow(items)
  max_score <- bundle$model$max_total_score

  observed <- vector("list", max_score + 1L)
  observed_gamma <- vector("list", max_score + 1L)
  expected <- vector("list", max_score + 1L)
  for (score in seq.int(0L, max_score)) {
    observed[[score + 1L]] <- vector("list", n_items)
    observed_gamma[[score + 1L]] <- vector("list", n_items)
    expected[[score + 1L]] <- vector("list", n_items)
    for (item_index in seq_len(n_items)) {
      item_bins <- items$raw_max[[item_index]]
      observed[[score + 1L]][[item_index]] <- numeric(item_bins)
      observed_gamma[[score + 1L]][[item_index]] <- numeric(item_bins)
      expected[[score + 1L]][[item_index]] <- numeric(item_bins)
    }
  }

  for (score in seq.int(0L, max_score)) {
    source_score_rows <- observed_rows & item_scores == score
    gamma_score_rows <- data$status == 1L & data$score == score
    for (item_index in seq_len(n_items)) {
      item_name <- items$name[[item_index]]
      observed[[score + 1L]][[item_index]] <- tabulate(
        data[[item_name]][source_score_rows] + 1L,
        nbins = items$raw_max[[item_index]]
      )

      item_values <- data[[item_name]][gamma_score_rows]
      item_values <- item_values[item_values >= 0L & item_values < items$raw_max[[item_index]]]
      observed_gamma[[score + 1L]][[item_index]] <- tabulate(
        item_values + 1L,
        nbins = items$raw_max[[item_index]]
      )

      moments <- conditional[[item_index]][[score + 1L]]
      if (score_weights[[score + 1L]] > 0L && length(moments$probabilities) > 0L) {
        expected[[score + 1L]][[item_index]] <-
          score_weights[[score + 1L]] * moments$probabilities
      }
    }
  }

  global <- extended_global_restscore_tables(bundle, conditional, observed_gamma, base_score_counts)
  local <- extended_local_restscore_tables(bundle, observed, expected, base_score_counts)
  list(
    restscore_tables = global$rows,
    local_restscore_tables = local$table_rows,
    local_gamma = list(rows = local$gamma_rows, limits = local$limits)
  )
}

#' Internal restscore matrix helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param item_max Internal `item_max` value used by this helper.
#' @param rest_max Internal `rest_max` value used by this helper.
#' @return The internal `restscore_matrix()` computation result.
#' @keywords internal
#' @noRd
restscore_matrix <- function(item_max, rest_max) {
  matrix(0, nrow = item_max + 1L, ncol = rest_max + 1L)
}

#' Internal extended global restscore tables helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle Source-shaped analysis bundle.
#' @param conditional Internal `conditional` value used by this helper.
#' @param observed_gamma Internal `observed_gamma` value used by this helper.
#' @param base_score_counts Internal `base_score_counts` value used by this helper.
#' @return The internal `extended_global_restscore_tables()` computation result.
#' @keywords internal
#' @noRd
extended_global_restscore_tables <- function(bundle, conditional, observed_gamma, base_score_counts) {
  items <- bundle$model$items
  max_score <- bundle$model$max_total_score
  counts <- as.numeric(base_score_counts[as.character(seq.int(0L, max_score))])
  counts[is.na(counts)] <- 0
  rows <- list()

  for (item_index in seq_len(nrow(items))) {
    item_max <- items$raw_max[[item_index]] - 1L
    rest_max <- max_score - item_max
    observed_table <- restscore_matrix(item_max, rest_max)
    expected_table <- restscore_matrix(item_max, rest_max)

    # Source trace: CalculateItemRestScoreGammaValues seeds only the
    # deterministic all-low/all-high cells from the base complete score counts.
    observed_table[1L, 1L] <- counts[[1L]]
    expected_table[1L, 1L] <- counts[[1L]]
    observed_table[item_max + 1L, rest_max + 1L] <- counts[[max_score + 1L]]
    expected_table[item_max + 1L, rest_max + 1L] <- counts[[max_score + 1L]]

    for (score in seq.int(1L, max_score - 1L)) {
      if (counts[[score + 1L]] <= 0) {
        next
      }
      probabilities <- conditional[[item_index]][[score + 1L]]$probabilities
      if (length(probabilities) == 0L) {
        next
      }
      for (item_score in seq.int(0L, item_max)) {
        rest_score <- score - item_score
        if (rest_score < 0L || rest_score > rest_max) {
          next
        }
        observed_table[item_score + 1L, rest_score + 1L] <-
          observed_table[item_score + 1L, rest_score + 1L] +
          observed_gamma[[score + 1L]][[item_index]][[item_score + 1L]]
        # Source trace: CalculateItemRestScoreGammaValues uses Counts.ScoreCounts
        # as the base source score count for expected cells, not the later
        # local-table ObsTotal renormalization.
        expected_table[item_score + 1L, rest_score + 1L] <-
          expected_table[item_score + 1L, rest_score + 1L] +
          counts[[score + 1L]] * probabilities[[item_score + 1L]]
      }
    }

    rows[[length(rows) + 1L]] <- restscore_table_rows(
      items,
      item_index,
      NA_integer_,
      observed_table,
      expected_table
    )
  }

  list(rows = type.convert(do.call(rbind, rows), as.is = TRUE))
}

#' Internal extended local restscore tables helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle Source-shaped analysis bundle.
#' @param observed Internal `observed` value used by this helper.
#' @param expected Internal `expected` value used by this helper.
#' @param base_score_counts Internal `base_score_counts` value used by this helper.
#' @return The internal `extended_local_restscore_tables()` computation result.
#' @keywords internal
#' @noRd
extended_local_restscore_tables <- function(bundle, observed, expected, base_score_counts) {
  items <- bundle$model$items
  max_score <- bundle$model$max_total_score
  counts <- as.numeric(base_score_counts[as.character(seq.int(0L, max_score))])
  counts[is.na(counts)] <- 0
  max_item_score <- max(items$raw_max - 1L)
  if (max_item_score < 2L) {
    empty_tables <- data.frame(
      item_label = character(),
      item_name = character(),
      local_restscore = integer(),
      restscore_table = character(),
      restscore = integer(),
      item_score = integer(),
      value = numeric(),
      stringsAsFactors = FALSE
    )
    empty_gamma <- data.frame(
      item_label = character(),
      item_name = character(),
      local_restscore = integer(),
      gamma_observed = numeric(),
      gamma_expected = numeric(),
      gamma_sd = numeric(),
      gamma_p = numeric(),
      gamma_risk = integer(),
      observed_n = numeric(),
      expected_n = numeric(),
      stringsAsFactors = FALSE
    )
    empty_limits <- data.frame(
      local_restscore = integer(),
      bh_05 = numeric(),
      bh_01 = numeric(),
      bh_001 = numeric(),
      stringsAsFactors = FALSE
    )
    return(list(table_rows = empty_tables, gamma_rows = empty_gamma, limits = empty_limits))
  }

  table_rows <- list()
  gamma_rows <- list()
  limit_rows <- list()

  for (adjacent_score in seq.int(0L, max_item_score - 1L)) {
    p_values <- numeric(nrow(items))
    adjacent_gamma <- vector("list", nrow(items))

    for (item_index in seq_len(nrow(items))) {
      item_result <- extended_local_restscore_item(
        items = items,
        item_index = item_index,
        adjacent_score = adjacent_score,
        max_score = max_score,
        counts = counts,
        observed = observed,
        expected = expected
      )
      table_rows[[length(table_rows) + 1L]] <- item_result$table_rows
      p_values[[item_index]] <- item_result$gamma_row$gamma_p
      adjacent_gamma[[item_index]] <- item_result$gamma_row
    }

    risks <- item_fdr_risk(p_values)
    bh_05 <- source_bh_critical(p_values, 0.05)
    bh_01 <- source_bh_critical(p_values, 0.01)
    bh_001 <- source_bh_critical(p_values, 0.001)
    limit_rows[[length(limit_rows) + 1L]] <- data.frame(
      local_restscore = adjacent_score,
      bh_05 = bh_05,
      bh_01 = bh_01,
      bh_001 = bh_001,
      stringsAsFactors = FALSE
    )

    for (item_index in seq_len(nrow(items))) {
      adjacent_gamma[[item_index]]$gamma_risk <- risks[[item_index]]
      gamma_rows[[length(gamma_rows) + 1L]] <- adjacent_gamma[[item_index]]
    }
  }

  list(
    table_rows = type.convert(do.call(rbind, table_rows), as.is = TRUE),
    gamma_rows = type.convert(do.call(rbind, gamma_rows), as.is = TRUE),
    limits = type.convert(do.call(rbind, limit_rows), as.is = TRUE)
  )
}

#' Build one item/adjacent-score local restscore diagnostic
#'
#' Source trace: `source/PAS_skunits/skbias14.pas::Calculate_item_restscore_gamma1`.
#' Mathematical step: materialize the adjacent-category item-by-restscore
#' observed and expected tables, rescale expected cells to observed score
#' totals, and derive the fitted-gamma comparison in source loop order.
#' @param items Parsed item metadata.
#' @param item_index One-based item index.
#' @param adjacent_score Lower category of the adjacent pair.
#' @param max_score Maximum total score.
#' @param counts Complete-record score counts including deterministic extremes.
#' @param observed Observed item margins by score and item.
#' @param expected Expected item margins by score and item.
#' @return Table rows and one local-gamma result row.
#' @keywords internal
#' @noRd
extended_local_restscore_item <- function(items,
                                          item_index,
                                          adjacent_score,
                                          max_score,
                                          counts,
                                          observed,
                                          expected) {
  item_max <- items$raw_max[[item_index]] - 1L
  rest_max <- max_score - item_max
  observed_table <- matrix(0, nrow = item_max + 1L, ncol = max_score + 1L)
  expected_table <- matrix(0, nrow = item_max + 1L, ncol = max_score + 1L)

  # Local adjacent tables inherit exactly one deterministic boundary: low only
  # for the 0/1 slice and high only for the top adjacent slice of that item.
  if (adjacent_score == 0L) {
    observed_table[1L, 1L] <- counts[[1L]]
    expected_table[1L, 1L] <- counts[[1L]]
  } else if (adjacent_score == item_max - 1L) {
    observed_table[item_max + 1L, rest_max + 1L] <- counts[[max_score + 1L]]
    expected_table[item_max + 1L, rest_max + 1L] <- counts[[max_score + 1L]]
  }

  observed_n <- 0
  expected_n <- 0
  for (score in seq.int(1L, max_score - 1L)) {
    for (item_score in seq.int(adjacent_score, adjacent_score + 1L)) {
      if (item_score > item_max || item_score > score) {
        next
      }
      rest_score <- score - item_score
      if (rest_score < 0L || rest_score > rest_max) {
        next
      }

      observed_counts <- observed[[score + 1L]][[item_index]]
      expected_counts <- expected[[score + 1L]][[item_index]]
      observed_total <- sum(observed_counts)
      expected_total <- sum(expected_counts)
      if (observed_total > 0) {
        observed_table[item_score + 1L, rest_score + 1L] <-
          observed_table[item_score + 1L, rest_score + 1L] +
          observed_counts[[item_score + 1L]]
        observed_n <- observed_n + observed_table[item_score + 1L, rest_score + 1L]
      }
      if (expected_total > 0) {
        # Calculate_item_restscore_gamma1 uses
        # (Expected[cell] / ExpTotal) * ObsTotal.
        expected_table[item_score + 1L, rest_score + 1L] <-
          expected_table[item_score + 1L, rest_score + 1L] +
          (expected_counts[[item_score + 1L]] / expected_total) * observed_total
        expected_n <- expected_n + expected_table[item_score + 1L, rest_score + 1L]
      }
    }
  }

  table_rows <- restscore_table_rows(
    items,
    item_index,
    adjacent_score,
    observed_table[, seq_len(rest_max + 1L), drop = FALSE],
    expected_table[, seq_len(rest_max + 1L), drop = FALSE]
  )
  gamma_observed <- goodman_kruskal_gamma(observed_table)
  fitted <- fitted_gamma_stats(expected_table)
  gamma_expected <- fitted$gamma
  gamma_sd <- sqrt(fitted$variance)
  gamma_p <- if (gamma_sd > 0) {
    two_sided_source_normal_p((gamma_observed - gamma_expected) / gamma_sd)
  } else {
    1
  }
  list(
    table_rows = table_rows,
    gamma_row = data.frame(
      item_label = items$label_code[[item_index]],
      item_name = items$name[[item_index]],
      local_restscore = adjacent_score,
      gamma_observed = gamma_observed,
      gamma_expected = gamma_expected,
      gamma_sd = gamma_sd,
      gamma_p = gamma_p,
      observed_n = observed_n,
      expected_n = expected_n,
      stringsAsFactors = FALSE
    )
  )
}

#' Internal restscore table rows helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param items Item selection or item metadata.
#' @param item_index One-based item index.
#' @param local_restscore Internal `local_restscore` value used by this helper.
#' @param observed_table Internal `observed_table` value used by this helper.
#' @param expected_table Internal `expected_table` value used by this helper.
#' @return The internal `restscore_table_rows()` computation result.
#' @keywords internal
#' @noRd
restscore_table_rows <- function(items, item_index, local_restscore, observed_table, expected_table) {
  tables <- list(observed = observed_table, expected = expected_table)
  item_rows <- lapply(names(tables), function(table_name) {
    tab <- tables[[table_name]]
    item_scores <- seq.int(0L, nrow(tab) - 1L)
    rest_scores <- seq.int(0L, ncol(tab) - 1L)
    n_cells <- length(item_scores) * length(rest_scores)
    row <- data.frame(
      item_label = rep(items$label_code[[item_index]], n_cells),
      item_name = rep(items$name[[item_index]], n_cells),
      restscore_table = rep(table_name, n_cells),
      restscore = rep(rest_scores, times = length(item_scores)),
      item_score = rep(item_scores, each = length(rest_scores)),
      value = as.vector(t(tab)),
      stringsAsFactors = FALSE
    )
    if (!is.na(local_restscore)) {
      row$local_restscore <- local_restscore
      row <- row[c(
        "item_label",
        "item_name",
        "local_restscore",
        "restscore_table",
        "restscore",
        "item_score",
        "value"
      )]
    }
    row
  })
  do.call(rbind, item_rows)
}
