#' Internal ari values helper
#'
#' Supports the ari values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_Ari`.
#' @param fit Fitted gRm model.
#' @return The internal `ari_values()` computation result.
#' @keywords internal
#' @noRd
ari_values <- function(fit) {
  bundle <- fit$bundle
  observed <- ari_observed_tables(bundle)
  expected <- ari_expected_moments(fit)
  ari_row_table(bundle, observed, expected)
}

#' Internal ari observed tables helper
#'
#' Supports the ari values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_Ari`.
#' @param bundle Source-shaped analysis bundle.
#' @return The internal `ari_observed_tables()` computation result.
#' @keywords internal
#' @noRd
ari_observed_tables <- function(bundle) {
  items <- bundle$model$items
  data <- bundle$data
  n_items <- nrow(items)
  global_item_max <- ari_global_item_max(items)
  max_total_score <- as.integer(bundle$model$max_total_score)
  counts <- array(
    0,
    dim = c(max_total_score + 1L, n_items, global_item_max + 1L),
    dimnames = list(
      score = as.character(seq.int(0L, max_total_score)),
      item = items$name,
      item_score = as.character(seq.int(0L, global_item_max))
    )
  )

  if (max_total_score < 2L || !nrow(data)) {
    return(list(counts = counts, score_min = 1L, score_max = max_total_score - 1L))
  }

  row_weight <- if ("count" %in% names(data)) {
    as.numeric(data$count)
  } else {
    rep(1, nrow(data))
  }
  row_weight[is.na(row_weight)] <- 0

  valid <- data$status == 1L &
    data$score >= 1L &
    data$score <= max_total_score - 1L &
    row_weight > 0

  for (item_index in seq_len(n_items)) {
    values <- data[[items$name[[item_index]]]]
    valid <- valid &
      !is.na(values) &
      values >= 0L &
      values <= items$raw_max[[item_index]] - 1L
  }

  rows <- which(valid)
  # Calculate_Ari.Count_Observed traverses accepted records, then items, and
  # adds the record frequency to the score/item/category cell in that order.
  for (row in rows) {
    score <- as.integer(data$score[[row]])
    weight <- row_weight[[row]]
    for (item_index in seq_len(n_items)) {
      item_score <- as.integer(data[[items$name[[item_index]]]][[row]])
      counts[score + 1L, item_index, item_score + 1L] <-
        counts[score + 1L, item_index, item_score + 1L] + weight
    }
  }

  list(counts = counts, score_min = 1L, score_max = max_total_score - 1L)
}

#' Internal ari expected moments helper
#'
#' Supports the ari values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_Ari`.
#' @param fit Fitted gRm model.
#' @return The internal `ari_expected_moments()` computation result.
#' @keywords internal
#' @noRd
ari_expected_moments <- function(fit) {
  if (is_gllrm_public_fit(fit)) {
    context <- fit$fit$context
    state <- fit$fit
    state$context <- NULL
    return(gllrm_item_conditional_moments(
      context,
      state,
      include_probabilities = TRUE,
      probability_cache = new_gllrm_probability_cache(context, state)
    ))
  }

  item_conditional_moments(
    fit$bundle,
    fit$fit$item_gamma,
    include_probabilities = TRUE
  )
}

#' Internal ari row table helper
#'
#' Supports the ari values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_Ari`.
#' @param bundle Source-shaped analysis bundle.
#' @param observed Internal `observed` value used by this helper.
#' @param expected Internal `expected` value used by this helper.
#' @return The internal `ari_row_table()` computation result.
#' @keywords internal
#' @noRd
ari_row_table <- function(bundle, observed, expected) {
  items <- bundle$model$items
  n_items <- nrow(items)
  global_item_max <- ari_global_item_max(items)
  columns <- ari_table_columns(global_item_max)
  rows <- list()

  if (observed$score_max < observed$score_min) {
    return(ari_empty_table(global_item_max))
  }

  for (item_index in seq_len(n_items)) {
    supported_scores <- seq.int(0L, items$raw_max[[item_index]] - 1L)
    for (score in seq.int(observed$score_min, observed$score_max)) {
      counts <- observed$counts[score + 1L, item_index, ]
      n <- sum(counts)
      if (n <= 0) {
        next
      }

      observed_prob <- as.numeric(counts) / n
      observed_supported <- observed_prob[seq_along(supported_scores)]
      obs_mean <- sum(supported_scores * observed_supported)
      raw_obs_var <- sum(supported_scores^2 * observed_supported) - obs_mean^2
      # Calculate_Ari.SaveARI prints the unusual (n-1)/n multiplier. Retain it
      # exactly; this is not the usual unbiased n/(n-1) variance correction.
      obs_var <- if (n > 1) ((n - 1) / n) * raw_obs_var else 0

      moment <- expected[[item_index]][[score + 1L]]
      expected_prob <- numeric(global_item_max + 1L)
      moment_prob <- as.numeric(moment$probabilities %||% numeric())
      if (length(moment_prob)) {
        expected_prob[seq_along(moment_prob)] <- moment_prob
      }
      exp_mean <- as.numeric(moment$mean %||% 0)
      exp_var <- as.numeric(moment$variance %||% 0)
      # Calculate_Ari.SaveARI standardizes the observed-minus-expected mean by
      # the fitted item variance at this total score, with a zero-variance
      # sentinel z of zero.
      z <- if (isTRUE(exp_var > 0)) {
        sqrt(n) * (obs_mean - exp_mean) / sqrt(exp_var)
      } else {
        0
      }

      values <- c(
        list(
          ItemNo = item_index,
          Item = items$name[[item_index]],
          Score = score,
          n = n
        ),
        stats::setNames(as.list(observed_prob), paste0("Obs", seq.int(0L, global_item_max))),
        list(ObsMean = obs_mean, ObsVar = obs_var),
        stats::setNames(as.list(expected_prob), paste0("Exp", seq.int(0L, global_item_max))),
        list(ExpMean = exp_mean, ExpVar = exp_var, z = z)
      )
      rows[[length(rows) + 1L]] <- as.data.frame(values[columns], stringsAsFactors = FALSE)
    }
  }

  if (!length(rows)) {
    return(ari_empty_table(global_item_max))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Internal ari global item max helper
#'
#' Supports the ari values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_Ari`.
#' @param items Item selection or item metadata.
#' @return The internal `ari_global_item_max()` computation result.
#' @keywords internal
#' @noRd
ari_global_item_max <- function(items) {
  max(as.integer(items$raw_max), 0L) - 1L
}

#' Internal ari table columns helper
#'
#' Supports the ari values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_Ari`.
#' @param global_item_max Internal `global_item_max` value used by this helper.
#' @return The internal `ari_table_columns()` computation result.
#' @keywords internal
#' @noRd
ari_table_columns <- function(global_item_max) {
  c(
    "ItemNo", "Item", "Score", "n",
    paste0("Obs", seq.int(0L, global_item_max)),
    "ObsMean", "ObsVar",
    paste0("Exp", seq.int(0L, global_item_max)),
    "ExpMean", "ExpVar", "z"
  )
}

#' Internal ari empty table helper
#'
#' Supports the ari values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_Ari`.
#' @param global_item_max Internal `global_item_max` value used by this helper.
#' @return The internal `ari_empty_table()` computation result.
#' @keywords internal
#' @noRd
ari_empty_table <- function(global_item_max) {
  columns <- ari_table_columns(global_item_max)
  out <- data.frame(
    ItemNo = integer(),
    Item = character(),
    Score = integer(),
    n = numeric(),
    stringsAsFactors = FALSE
  )
  for (name in setdiff(columns, names(out))) {
    out[[name]] <- numeric()
  }
  out[columns]
}
