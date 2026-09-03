#' Internal gllrm component restscore tables helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param component_items Internal `component_items` value used by this helper.
#' @return The internal `gllrm_component_restscore_tables()` computation result.
#' @keywords internal
#' @noRd
gllrm_component_restscore_tables <- function(context, state, component_items) {
  component_max <- sum(context$item_raw_max[component_items] - 1L)
  rest_max <- context$max_total_score - component_max
  source_score_min <- 1L
  source_score_max <- context$max_total_score - 1L
  observed <- matrix(0, nrow = component_max + 1L, ncol = rest_max + 1L)
  expected <- matrix(0, nrow = component_max + 1L, ncol = rest_max + 1L)

  # Source trace: skbias15 Item_fit_analysis calls Calculate_item_fits with
  # LeastScore=1 and LargestScore=highest_possible_score-1. Count_Observed
  # updates Nlowscore/Nhighscore from scores outside that interval before
  # rejecting rows on exogenous values, and
  # Calculate_comp_restscore_gamma seeds both observed and expected endpoint
  # cells from those counts, then loops scores 1..largest_possible_score-1.
  # Interior component-restscore cells come from the valid ScoreTab/component
  # margin pass.
  complete_item_rows <- which(rowSums(context$item_matrix < 0L) == 0L)
  if (length(complete_item_rows) > 0L) {
    complete_scores <- rowSums(context$item_matrix[complete_item_rows, , drop = FALSE])
    low_count <- sum(complete_scores < source_score_min)
    high_count <- sum(complete_scores > source_score_max)
    observed[1L, 1L] <- low_count
    observed[component_max + 1L, rest_max + 1L] <- high_count
    expected[1L, 1L] <- low_count
    expected[component_max + 1L, rest_max + 1L] <- high_count
  }

  interior_rows <- context$valid_rows[
    context$score[context$valid_rows] >= source_score_min &
      context$score[context$valid_rows] <= source_score_max
  ]
  if (length(interior_rows) > 0L) {
    component_scores <- rowSums(context$item_matrix[interior_rows, component_items, drop = FALSE])
    rest_scores <- context$score[interior_rows] - component_scores
    keep <- component_scores >= 0L & component_scores <= component_max &
      rest_scores >= 0L & rest_scores <= rest_max
    if (any(keep)) {
      observed[] <- observed + tabulate(
        component_scores[keep] + rest_scores[keep] * (component_max + 1L) + 1L,
        nbins = length(observed)
      )
    }
  }

  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  component_index <- match(
    gllrm_component_key(component_items),
    vapply(components, gllrm_component_key, character(1L))
  )
  if (is.na(component_index)) {
    return(list(observed = observed, expected = expected))
  }

  background_cache <- new.env(parent = emptyenv())
  for (group_index in seq_len(nrow(context$score_exo_groups))) {
    group <- context$score_exo_groups[group_index, , drop = FALSE]
    total_score <- group$score[[1L]]
    if (total_score < source_score_min || total_score > source_score_max) {
      next
    }
    background_values <- gllrm_group_background_values(context, group)
    cache_key <- paste(background_values, collapse = "\r")
    if (!exists(cache_key, envir = background_cache, inherits = FALSE)) {
      component_gamma <- lapply(components, function(items) {
        gllrm_component_gamma(context, state, items, background_values)
      })
      assign(
        cache_key,
        list(
          component_gamma = component_gamma,
          convolutions = gllrm_component_convolutions(
            lapply(component_gamma, `[[`, "gamma"),
            context$max_total_score
          )
        ),
        envir = background_cache
      )
    }
    cached <- get(cache_key, envir = background_cache, inherits = FALSE)
    full_gamma <- cached$convolutions$full
    denominator <- full_gamma[[total_score + 1L]]
    if (denominator <= 0) {
      next
    }
    component_info <- cached$component_gamma[[component_index]]
    rest_gamma <- cached$convolutions$rest[[component_index]]
    key <- gllrm_component_key(component_items)
    config_matrix <- context$component_config_matrices[[key]]
    config_scores <- context$component_config_scores[[key]]
    if (is.null(config_matrix)) {
      configs <- gllrm_component_configurations(context, component_items)
      config_matrix <- as.matrix(configs)
      storage.mode(config_matrix) <- "integer"
      config_scores <- rowSums(config_matrix)
    }
    for (config_index in seq_len(nrow(config_matrix))) {
      component_score <- config_scores[[config_index]]
      rest_score <- total_score - component_score
      if (rest_score < 0L || rest_score > rest_max) {
        next
      }
      config_weight <- component_info$config_weights[[config_index]]
      if (config_weight <= 0) {
        next
      }
      expected[component_score + 1L, rest_score + 1L] <-
        expected[component_score + 1L, rest_score + 1L] +
        group$count[[1L]] * (config_weight / component_info$scale) *
          rest_gamma[[rest_score + 1L]] / denominator
    }
  }

  list(observed = observed, expected = expected)
}

#' Internal component restscore table rows helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param component Internal `component` value used by this helper.
#' @param representative_item_label Internal `representative_item_label` value used by this helper.
#' @param observed_table Internal `observed_table` value used by this helper.
#' @param expected_table Internal `expected_table` value used by this helper.
#' @return The internal `component_restscore_table_rows()` computation result.
#' @keywords internal
#' @noRd
component_restscore_table_rows <- function(component,
                                           representative_item_label,
                                           observed_table,
                                           expected_table) {
  tables <- list(observed = observed_table, expected = expected_table)
  out <- lapply(names(tables), function(table_name) {
    tab <- tables[[table_name]]
    component_scores <- seq.int(0L, nrow(tab) - 1L)
    rest_scores <- seq.int(0L, ncol(tab) - 1L)
    n_cells <- length(component_scores) * length(rest_scores)
    data.frame(
      component = rep(component, n_cells),
      representative_item_label = rep(representative_item_label, n_cells),
      restscore_table = rep(table_name, n_cells),
      restscore = rep(rest_scores, each = length(component_scores)),
      component_score = rep(component_scores, times = length(rest_scores)),
      value = as.vector(tab),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

#' Internal gllrm item conditional moments helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param include_probabilities Internal `include_probabilities` value used by this helper.
#' @param probability_cache Internal `probability_cache` value used by this helper.
#' @return The internal `gllrm_item_conditional_moments()` computation result.
#' @keywords internal
#' @noRd
gllrm_item_conditional_moments <- function(context,
                                                  state,
                                                  include_probabilities = FALSE,
                                                  probability_cache = NULL) {
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  probability_cache <- probability_cache %||%
    new_gllrm_probability_cache(context, state, components = components)
  probability_by_item <- lapply(seq_len(context$n_items), function(item_index) {
    lapply(
      seq.int(0L, context$max_total_score),
      function(score) numeric(length(context$item_score_values[[item_index]]))
    )
  })

  for (score in seq.int(0L, context$max_total_score)) {
    group_rows <- context$score_exo_groups[
      context$score_exo_groups$score == score,
      ,
      drop = FALSE
    ]
    total_count <- sum(group_rows$count)
    if (total_count <= 0L) {
      next
    }

    for (group_index in seq_len(nrow(group_rows))) {
      group <- group_rows[group_index, , drop = FALSE]
      background_values <- gllrm_group_background_values(context, group)
      group_probabilities <- gllrm_cached_item_probabilities(
        probability_cache,
        total_score = score,
        background_values = background_values
      )
      for (item_index in seq_len(context$n_items)) {
        probability_by_item[[item_index]][[score + 1L]] <-
          probability_by_item[[item_index]][[score + 1L]] +
          group$count[[1L]] * group_probabilities[[item_index]]
      }
    }

    for (item_index in seq_len(context$n_items)) {
      probabilities <- probability_by_item[[item_index]][[score + 1L]] / total_count
      probability_sum <- sum(probabilities)
      if (probability_sum > 0) {
        probabilities <- probabilities / probability_sum
      }
      probability_by_item[[item_index]][[score + 1L]] <- probabilities
    }
  }

  result <- vector("list", context$n_items)
  for (item_index in seq_len(context$n_items)) {
    item_result <- vector("list", context$max_total_score + 1L)
    item_scores <- context$item_score_values[[item_index]]
    for (score in seq.int(0L, context$max_total_score)) {
      probabilities <- probability_by_item[[item_index]][[score + 1L]]
      mean_value <- sum(item_scores * probabilities)
      centered <- item_scores - mean_value
      variance <- sum(centered^2 * probabilities)
      fourth <- sum(centered^4 * probabilities)
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

#' Internal gllrm group item probabilities helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param total_score Internal `total_score` value used by this helper.
#' @param background_values Internal `background_values` value used by this helper.
#' @param components Internal `components` value used by this helper.
#' @return The internal `gllrm_group_item_probabilities()` computation result.
#' @keywords internal
#' @noRd
gllrm_group_item_probabilities <- function(context,
                                                  state,
                                                  total_score,
                                                  background_values,
                                                  components = NULL) {
  components <- components %||% context$ld_components_items %||% gllrm_ld_components(context)$items
  out <- lapply(context$item_raw_max, numeric)
  if (total_score < 0L || total_score > context$max_total_score) {
    return(out)
  }

  component_gamma <- lapply(components, function(component_items) {
    one <- gllrm_component_gamma(context, state, component_items, background_values)
    one$gamma * one$scale
  })
  convolutions <- gllrm_component_convolutions(component_gamma, context$max_total_score)
  full_gamma <- convolutions$full
  denominator <- full_gamma[[total_score + 1L]]
  if (denominator <= 0) {
    return(out)
  }

  for (component_index in seq_along(components)) {
    component_items <- components[[component_index]]
    rest_gamma <- convolutions$rest[[component_index]]

    key <- gllrm_component_key(component_items)
    config_matrix <- context$component_config_matrices[[key]]
    config_scores <- context$component_config_scores[[key]]
    if (is.null(config_matrix)) {
      configs <- gllrm_component_configurations(context, component_items)
      config_matrix <- as.matrix(configs)
      storage.mode(config_matrix) <- "integer"
      config_scores <- rowSums(config_matrix)
    }

    weights <- gllrm_component_config_weights_fast(
      context,
      state,
      component_items,
      config_matrix,
      background_values,
      key = key
    )
    valid <- config_scores <= total_score & weights > 0
    if (!any(valid)) {
      next
    }
    valid_index <- which(valid)
    rest_values <- rest_gamma[total_score - config_scores[valid_index] + 1L]
    probabilities <- weights[valid_index] * rest_values / denominator
    positive <- probabilities > 0
    if (!any(positive)) {
      next
    }
    valid_index <- valid_index[positive]
    probabilities <- probabilities[positive]

    for (local_index in seq_along(component_items)) {
      item_index <- component_items[[local_index]]
      item_scores <- config_matrix[valid_index, local_index] + 1L
      by_score <- rowsum(probabilities, item_scores, reorder = FALSE)
      score_index <- as.integer(rownames(by_score))
      out[[item_index]][score_index] <-
        out[[item_index]][score_index] + as.numeric(by_score[, 1L])
    }
  }

  out
}
