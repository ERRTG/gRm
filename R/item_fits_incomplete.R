#' Collect source-style incomplete response records
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle An item-parameters bundle.
#' @return A list with grouped incomplete records, counts, observed scores, and
#'   observed maximum scores.
#' @keywords internal
#' @noRd
empty_source_incomplete_records <- function(bundle) {
  records <- as.data.frame(
    stats::setNames(
      replicate(nrow(bundle$model$items), integer(), simplify = FALSE),
      bundle$model$items$name
    ),
    stringsAsFactors = FALSE
  )
  list(
    records = records,
    count = integer(),
    score = integer(),
    max_score = integer()
  )
}

#' Internal collect source incomplete records helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle Source-shaped analysis bundle.
#' @return The internal `collect_source_incomplete_records()` computation result.
#' @keywords internal
#' @noRd
collect_source_incomplete_records <- function(bundle) {
  # Source trace: SKbias2.pas::collect_incomplete_response_records. DIGRAM keeps
  # incomplete rows only when some but not all item responses are missing, at
  # least two observed items remain, and the observed score is non-extreme for
  # the observed subset.
  items <- bundle$model$items
  backgrounds <- bundle$model$backgrounds
  data <- bundle$data
  records <- list()
  counts <- integer()
  scores <- integer()
  max_scores <- integer()

  for (row_index in seq_len(nrow(data))) {
    if (nrow(backgrounds) > 0L &&
      any(unlist(data[row_index, backgrounds$name, drop = FALSE], use.names = FALSE) < 0L)) {
      next
    }

    item_values <- as.integer(unlist(data[row_index, items$name, drop = FALSE], use.names = FALSE))
    missing <- item_values > (items$raw_max - 1L) | item_values < 0L
    n_missing <- sum(missing)
    if (n_missing == 0L || n_missing == nrow(items)) {
      next
    }

    observed <- !missing
    observed_score <- sum(item_values[observed])
    observed_max <- sum(items$raw_max[observed] - 1L)
    if (sum(observed) < 2L || observed_score == 0L || observed_score == observed_max) {
      next
    }

    source_record <- item_values
    source_record[missing] <- 99L
    key <- paste(source_record, collapse = "\t")
    index <- match(key, names(records))
    if (is.na(index)) {
      records[[key]] <- source_record
      counts[[length(counts) + 1L]] <- 1L
      scores[[length(scores) + 1L]] <- observed_score
      max_scores[[length(max_scores) + 1L]] <- observed_max
    } else {
      counts[[index]] <- counts[[index]] + 1L
    }
  }

  if (length(records) == 0L) {
    record_df <- as.data.frame(matrix(nrow = 0L, ncol = nrow(items)))
  } else {
    record_df <- as.data.frame(do.call(rbind, records))
  }
  names(record_df) <- items$name
  list(records = record_df, count = counts, score = scores, max_score = max_scores)
}

#' Calculate DIGRAM's expected total score for an incomplete record
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle An item-parameters bundle.
#' @param item_gamma Estimated item gamma matrix.
#' @param incomplete Output from `collect_source_incomplete_records()`.
#' @param record_index One-based incomplete record index.
#' @return Observed score plus expected score over missing items.
#' @keywords internal
#' @noRd
incomplete_expected_total_score <- function(bundle, item_gamma, incomplete, record_index) {
  # Source trace: skbias12a.pas::ExpectedIncompleteResponses first estimates the
  # person parameter from observed items, then evaluates the missing-items score
  # generating function at that theta and adds the observed score.
  items <- bundle$model$items
  record <- as.integer(incomplete$records[record_index, items$name, drop = TRUE])
  use_item <- record <= (items$raw_max - 1L)
  missing_item <- !use_item
  observed_score <- incomplete$score[[record_index]]
  observed_max <- incomplete$max_score[[record_index]]

  observed_gamma <- build_source_subset_gamma(bundle, item_gamma, use_item)
  theta <- estimate_person_parameter(observed_score, observed_max, observed_gamma, 1000L)

  missing_gamma <- build_source_subset_gamma(bundle, item_gamma, missing_item)
  missing_max <- sum(items$raw_max[missing_item] - 1L)
  observed_score + true_score_from_gamma(theta, missing_max, missing_gamma)
}

#' Build a source score-generating function for a subset of items
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle An item-parameters bundle.
#' @param item_gamma Estimated item gamma matrix.
#' @param use_item Logical vector selecting included items.
#' @return Numeric score-generating function indexed by score plus one.
#' @keywords internal
#' @noRd
build_source_subset_gamma <- function(bundle, item_gamma, use_item) {
  items <- bundle$model$items
  max_total_score <- sum(items$raw_max[use_item] - 1L)
  gamma_values <- numeric(max_total_score + 1L)
  gamma_values[[1L]] <- 1
  current_max <- 0L

  for (item_index in which(use_item)) {
    next_values <- numeric(max_total_score + 1L)
    next_max <- current_max + items$raw_max[[item_index]] - 1L
    for (score in seq.int(0L, current_max)) {
      current_weight <- gamma_values[[score + 1L]]
      if (current_weight == 0) {
        next
      }
      for (item_score in seq.int(0L, items$raw_max[[item_index]] - 1L)) {
        next_values[[score + item_score + 1L]] <-
          next_values[[score + item_score + 1L]] +
          current_weight * item_gamma[item_index, as.character(item_score)]
      }
    }
    gamma_values <- next_values
    current_max <- next_max
  }

  gamma_values
}

#' Prepare incomplete-record quantities for item-restscore gamma
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle An item-parameters bundle.
#' @param item_gamma Estimated item gamma matrix.
#' @param incomplete Output from `collect_source_incomplete_records()`.
#' @return A list of expected totals and per-item conditional probabilities.
#' @keywords internal
#' @noRd
prepare_incomplete_item_restscore_gamma <- function(bundle, item_gamma, incomplete) {
  # Source trace: skbias14/skbias15 call IncompleteResponseProbabilities inside
  # each item loop, but the value depends only on the incomplete record's
  # observed-item set and score. Caching preserves the algorithm while avoiding
  # repeated identical convolutions for BFI-sized data.
  n_records <- nrow(incomplete$records)
  items <- bundle$model$items
  expected_total <- numeric(n_records)
  probabilities <- vector("list", n_records)
  expected_cache <- new.env(parent = emptyenv())
  probability_cache <- new.env(parent = emptyenv())
  subset_gamma_cache <- new.env(parent = emptyenv())
  excluding_gamma_cache <- new.env(parent = emptyenv())

  cached_subset_gamma <- function(use_item) {
    key <- paste(as.integer(use_item), collapse = "")
    if (!exists(key, subset_gamma_cache, inherits = FALSE)) {
      assign(key, build_source_subset_gamma_fast(items, item_gamma, use_item), envir = subset_gamma_cache)
    }
    get(key, subset_gamma_cache, inherits = FALSE)
  }

  cached_excluding_gamma <- function(use_item, item_index) {
    key <- paste(paste(as.integer(use_item), collapse = ""), item_index, sep = ":")
    if (!exists(key, excluding_gamma_cache, inherits = FALSE)) {
      assign(
        key,
        build_source_subset_gamma_fast(items, item_gamma, replace(use_item, item_index, FALSE)),
        envir = excluding_gamma_cache
      )
    }
    get(key, excluding_gamma_cache, inherits = FALSE)
  }

  for (record_index in seq_len(n_records)) {
    record <- as.integer(incomplete$records[record_index, items$name, drop = TRUE])
    use_item <- as.logical(record <= (items$raw_max - 1L))
    expected_key <- paste(paste(as.integer(use_item), collapse = ""), incomplete$score[[record_index]], sep = ":")
    if (!exists(expected_key, expected_cache, inherits = FALSE)) {
      observed_gamma <- cached_subset_gamma(use_item)
      theta <- estimate_person_parameter(
        incomplete$score[[record_index]],
        incomplete$max_score[[record_index]],
        observed_gamma,
        1000L
      )
      missing_item <- !use_item
      missing_gamma <- cached_subset_gamma(missing_item)
      missing_max <- sum(items$raw_max[missing_item] - 1L)
      assign(
        expected_key,
        incomplete$score[[record_index]] + true_score_from_gamma(theta, missing_max, missing_gamma),
        envir = expected_cache
      )
    }
    expected_total[[record_index]] <- get(expected_key, expected_cache, inherits = FALSE)

    probabilities[[record_index]] <- vector("list", nrow(items))
    for (item_index in which(use_item)) {
      probability_key <- paste(expected_key, item_index, sep = ":")
      if (!exists(probability_key, probability_cache, inherits = FALSE)) {
        item_scores <- seq.int(0L, items$raw_max[[item_index]] - 1L)
        without_item <- cached_excluding_gamma(use_item, item_index)
        weights <- numeric(length(item_scores))
        for (offset in seq_along(item_scores)) {
          item_score <- item_scores[[offset]]
          if (incomplete$score[[record_index]] >= item_score) {
            weights[[offset]] <- item_gamma[item_index, as.character(item_score)] *
              without_item[[incomplete$score[[record_index]] - item_score + 1L]]
          }
        }
        denominator <- sum(weights)
        probs <- if (denominator > 0) weights / denominator else numeric(length(item_scores))
        assign(probability_key, probs, envir = probability_cache)
      }
      probabilities[[record_index]][[item_index]] <- get(probability_key, probability_cache, inherits = FALSE)
    }
  }

  list(expected_total = expected_total, probabilities = probabilities)
}

#' Internal item conditional moments helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle Source-shaped analysis bundle.
#' @param item_gamma Internal `item_gamma` value used by this helper.
#' @param include_probabilities Internal `include_probabilities` value used by this helper.
#' @return The internal `item_conditional_moments()` computation result.
#' @keywords internal
#' @noRd
item_conditional_moments <- function(bundle, item_gamma, include_probabilities = FALSE) {
  items <- bundle$model$items
  max_total_score <- bundle$model$max_total_score
  result <- vector("list", nrow(items))

  for (item_index in seq_len(nrow(items))) {
    without_item <- build_gamma_excluding_item(bundle, item_gamma, item_index)
    item_result <- vector("list", max_total_score + 1L)
    for (score in seq.int(0L, max_total_score)) {
      item_scores <- seq.int(0L, items$raw_max[[item_index]] - 1L)
      weights <- numeric(length(item_scores))
      for (offset in seq_along(item_scores)) {
        item_score <- item_scores[[offset]]
        if (score >= item_score) {
          weights[[offset]] <- item_gamma[item_index, as.character(item_score)] *
            without_item[[score - item_score + 1L]]
        }
      }
      denominator <- sum(weights)
      if (denominator > 0) {
        probabilities <- weights / denominator
        mean_value <- sum(item_scores * probabilities)
        centered <- item_scores - mean_value
        variance <- sum(centered^2 * probabilities)
        fourth <- sum(centered^4 * probabilities)
      } else {
        probabilities <- numeric(length(item_scores))
        mean_value <- 0
        variance <- 0
        fourth <- 0
      }
      item_result[[score + 1L]] <- list(
        mean = mean_value,
        variance = variance,
        fourth = fourth,
        probabilities = if (include_probabilities) probabilities else numeric()
      )
    }
    result[[item_index]] <- item_result
  }

  result
}

#' Internal item conditional moment for subset helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle Source-shaped analysis bundle.
#' @param item_gamma Internal `item_gamma` value used by this helper.
#' @param item_index One-based item index.
#' @param score Zero-based total or item score.
#' @param use_item Internal `use_item` value used by this helper.
#' @return The internal `item_conditional_moment_for_subset()` computation result.
#' @keywords internal
#' @noRd
item_conditional_moment_for_subset <- function(bundle, item_gamma, item_index, score, use_item) {
  items <- bundle$model$items
  item_scores <- seq.int(0L, items$raw_max[[item_index]] - 1L)
  without_item <- build_gamma_excluding_item_subset(bundle, item_gamma, item_index, use_item)
  weights <- numeric(length(item_scores))
  for (offset in seq_along(item_scores)) {
    item_score <- item_scores[[offset]]
    if (score >= item_score) {
      weights[[offset]] <- item_gamma[item_index, as.character(item_score)] *
        without_item[[score - item_score + 1L]]
    }
  }
  denominator <- sum(weights)
  if (denominator > 0) {
    probabilities <- weights / denominator
    mean_value <- sum(item_scores * probabilities)
    centered <- item_scores - mean_value
    variance <- sum(centered^2 * probabilities)
    fourth <- sum(centered^4 * probabilities)
  } else {
    probabilities <- numeric(length(item_scores))
    mean_value <- 0
    variance <- 0
    fourth <- 0
  }
  list(mean = mean_value, variance = variance, fourth = fourth, probabilities = probabilities)
}

#' Internal build gamma excluding item subset helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle Source-shaped analysis bundle.
#' @param item_gamma Internal `item_gamma` value used by this helper.
#' @param excluded_item Internal `excluded_item` value used by this helper.
#' @param use_item Internal `use_item` value used by this helper.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
build_gamma_excluding_item_subset <- function(bundle, item_gamma, excluded_item, use_item) {
  items <- bundle$model$items
  # Source trace: this is the same score-polynomial convolution as
  # Gamma_calculation1/IncompleteresponseProbabilities, but bounded to the
  # selected observed subset. The previous full-test bound was algorithmically
  # equivalent for values read by DIGRAM, but wasteful for incomplete BFI rows.
  included <- use_item
  included[[excluded_item]] <- FALSE
  max_total_score <- sum(items$raw_max[included] - 1L)
  gamma_values <- numeric(max_total_score + 1L)
  gamma_values[[1L]] <- 1
  current_max <- 0L

  for (item_index in seq_len(nrow(items))) {
    if (!use_item[[item_index]] || item_index == excluded_item) {
      next
    }
    next_values <- numeric(max_total_score + 1L)
    next_max <- current_max + items$raw_max[[item_index]] - 1L
    for (score in seq.int(0L, current_max)) {
      current_weight <- gamma_values[[score + 1L]]
      if (current_weight == 0) {
        next
      }
      for (item_score in seq.int(0L, items$raw_max[[item_index]] - 1L)) {
        next_values[[score + item_score + 1L]] <-
          next_values[[score + item_score + 1L]] +
          current_weight * item_gamma[item_index, as.character(item_score)]
      }
    }
    gamma_values[seq_len(next_max + 1L)] <- next_values[seq_len(next_max + 1L)]
    current_max <- next_max
  }

  gamma_values
}

#' Fast source score-generating function for a subset of items
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param items Item metadata data frame.
#' @param item_gamma Estimated item gamma matrix.
#' @param use_item Logical vector selecting included items.
#' @return Numeric score-generating function indexed by score plus one.
#' @keywords internal
#' @noRd
build_source_subset_gamma_fast <- function(items, item_gamma, use_item) {
  # Source trace: same polynomial convolution as Gamma_calculation1, expressed
  # with R's vectorized open convolution. This changes only execution strategy,
  # not the generated gamma polynomial.
  gamma_values <- 1
  for (item_index in which(use_item)) {
    item_values <- as.numeric(item_gamma[item_index, as.character(seq.int(0L, items$raw_max[[item_index]] - 1L))])
    gamma_values <- convolve(gamma_values, rev(item_values), type = "open")
  }
  gamma_values
}
