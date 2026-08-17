#' Internal gllrm global homogeneity item mean rows helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param context Prepared GLLRM computation context.
#' @param full_state Internal `full_state` value used by this helper.
#' @param group Internal `group` value used by this helper.
#' @param probability_cache Internal `probability_cache` value used by this helper.
#' @return The internal `gllrm_global_homogeneity_item_mean_rows()` computation result.
#' @keywords internal
#' @noRd
gllrm_global_homogeneity_item_mean_rows <- function(context,
                                                     full_state,
                                                     group,
                                                     probability_cache = NULL) {
  items <- context$items
  counts <- context$counts$item_counts
  expected_tables <- gllrm_global_homogeneity_expected_item_margin_tables(
    context,
    full_state,
    probability_cache = probability_cache
  )
  expected_variance <- gllrm_global_homogeneity_item_variances(
    context,
    expected_tables = expected_tables
  )
  rows <- vector("list", nrow(items))

  for (item_index in seq_len(nrow(items))) {
    item_max <- items$raw_max[[item_index]] - 1L
    item_scores <- seq.int(0L, item_max)
    score_cols <- item_scores + 1L
    observed_counts <- as.numeric(counts[item_index, as.character(item_scores)])
    n <- sum(observed_counts)
    observed_mean <- if (n > 0) {
      sum(item_scores * observed_counts) / n
    } else {
      0
    }
    expected_counts <- colSums(expected_tables[[item_index]])[score_cols]
    expected_n <- sum(expected_counts)
    expected_mean <- if (expected_n > 0) {
      sum(item_scores * expected_counts) / expected_n
    } else {
      0
    }
    item_expected_variance <- expected_variance[[item_index]]
    # The available Pascal source backs the GLLRM counts, expected
    # means, score-group likelihoods, and CLR, but it does not reproduce the
    # runtime denominator used for this compact residual column. Historical
    # harness traces in this repository end at the same boundary, so the
    # residual and marker are deliberately exposed as unmodeled values.
    residual <- NA_real_
    marker <- NA_character_
    rows[[item_index]] <- data.frame(
      group = as.integer(group),
      item_label = items$label_code[[item_index]],
      item_name = items$name[[item_index]],
      n = as.integer(n),
      observed_mean = observed_mean,
      expected_mean = expected_mean,
      expected_variance = item_expected_variance,
      visible_source_residual = residual,
      residual = residual,
      marker = marker,
      residual_runtime_source_backed = FALSE,
      marker_runtime_source_backed = FALSE,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

#' Internal gllrm global homogeneity expected item margin tables helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param context Prepared GLLRM computation context.
#' @param full_state Internal `full_state` value used by this helper.
#' @param probability_cache Internal `probability_cache` value used by this helper.
#' @return The internal `gllrm_global_homogeneity_expected_item_margin_tables()` computation result.
#' @keywords internal
#' @noRd
gllrm_global_homogeneity_expected_item_margin_tables <- function(context, full_state, probability_cache = NULL) {
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  probability_cache <- probability_cache %||%
    new_gllrm_probability_cache(context, full_state, components = components)
  tables <- lapply(seq_len(context$n_items), function(item_index) {
    item_max <- context$item_raw_max[[item_index]] - 1L
    matrix(
      0,
      nrow = context$max_total_score + 1L,
      ncol = item_max + 1L,
      dimnames = list(
        as.character(seq.int(0L, context$max_total_score)),
        as.character(seq.int(0L, item_max))
      )
    )
  })

  for (group_index in seq_len(nrow(context$score_exo_groups))) {
    group <- context$score_exo_groups[group_index, , drop = FALSE]
    total_score <- group$score[[1L]]
    background_values <- gllrm_group_background_values(context, group)
    probabilities <- gllrm_cached_item_probabilities(
      probability_cache,
      total_score = total_score,
      background_values = background_values
    )
    for (item_index in seq_len(context$n_items)) {
      p <- probabilities[[item_index]]
      tables[[item_index]][total_score + 1L, seq_along(p)] <-
        tables[[item_index]][total_score + 1L, seq_along(p)] +
        group$count[[1L]] * p
    }
  }

  tables
}

#' Internal gllrm global homogeneity item variances helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param context Prepared GLLRM computation context.
#' @param expected_tables Internal `expected_tables` value used by this helper.
#' @return The internal `gllrm_global_homogeneity_item_variances()` computation result.
#' @keywords internal
#' @noRd
gllrm_global_homogeneity_item_variances <- function(context, expected_tables) {
  out <- numeric(context$n_items)
  score_item_n <- gllrm_global_homogeneity_score_item_n(context)
  scfra <- if (min(context$score[context$valid_rows]) == 0L) {
    1L
  } else {
    min(context$score[context$valid_rows])
  }
  sctil <- if (max(context$score[context$valid_rows]) == context$max_total_score) {
    context$max_total_score - 1L
  } else {
    max(context$score[context$valid_rows])
  }

  for (item_index in seq_len(context$n_items)) {
    item_max <- context$item_raw_max[[item_index]] - 1L
    expected_table <- expected_tables[[item_index]]
    total_summary <- global_homogeneity_summarize_tal(
      colSums(expected_table),
      item_max
    )
    if (total_summary$n <= 0 || scfra > sctil) {
      next
    }
    for (score in seq.int(scfra, sctil)) {
      observed_score_n <- score_item_n[score + 1L, item_index]
      if (observed_score_n <= 0) {
        next
      }
      score_summary <- global_homogeneity_summarize_tal(
        expected_table[score + 1L, ],
        item_max
      )
      out[[item_index]] <- out[[item_index]] +
        (observed_score_n / total_summary$n) * score_summary$variance
    }
  }

  out
}

#' Internal gllrm global homogeneity score item n helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `gllrm_global_homogeneity_score_item_n()` computation result.
#' @keywords internal
#' @noRd
gllrm_global_homogeneity_score_item_n <- function(context) {
  out <- matrix(
    0,
    nrow = context$max_total_score + 1L,
    ncol = context$n_items,
    dimnames = list(
      as.character(seq.int(0L, context$max_total_score)),
      context$items$name
    )
  )
  for (item_index in seq_len(context$n_items)) {
    out[, item_index] <- tabulate(
      context$score[context$valid_rows] + 1L,
      nbins = context$max_total_score + 1L
    )
  }
  out
}
