#' @noRd
item_fits_values <- function(project, max_step = 5000L, max_delta = 0.0001, include_extended = TRUE) {
  if ((inherits(project, "gRm_fit") || inherits(project, "gRm_gllrm_fit")) && inherits(project$values, "gRm_active_gllrm_values")) {
    return(active_gllrm_item_fits_values(project, include_extended = include_extended))
  }
  bundle <- build_item_parameters_bundle(project)
  fit <- fit_rasch_base(bundle, max_step = max_step, max_delta = max_delta)
  conditional <- item_conditional_moments(bundle, fit$item_gamma, include_probabilities = TRUE)

  item_fit <- calculate_conditional_item_fit_values(bundle, fit, conditional = conditional)
  gamma_fit <- calculate_item_restscore_gamma_values(bundle, fit, conditional = conditional)
  extended <- if (include_extended) {
    calculate_extended_item_fit_values(bundle, fit, conditional = conditional)
  } else {
    NULL
  }

  rows <- data.frame(
    item_label = bundle$model$items$label,
    item_name = bundle$model$items$name,
    outfit = item_fit$outfit,
    outfit_sd = item_fit$outfit_sd,
    p_outfit = item_fit$p_outfit,
    outfit_fdr = item_fdr_risk(item_fit$p_outfit),
    infit = item_fit$infit,
    infit_sd = item_fit$infit_sd,
    p_infit = item_fit$p_infit,
    infit_fdr = item_fdr_risk(item_fit$p_infit),
    observed_gamma = gamma_fit$observed_gamma,
    expected_gamma = gamma_fit$expected_gamma,
    gamma_sd = gamma_fit$gamma_sd,
    p_gamma = gamma_fit$p_gamma,
    gamma_fdr = item_fdr_risk(gamma_fit$p_gamma),
    stringsAsFactors = FALSE
  )
  rows$direction <- item_fit_direction(rows)

  side_file <- data.frame(
    item = rows$item_label,
    outfit = rows$outfit,
    p_outfit = rows$p_outfit,
    infit = rows$infit,
    p_infit = rows$p_infit,
    ObsGamma = rows$observed_gamma,
    ExpGamma = rows$expected_gamma,
    p_gamma = rows$p_gamma,
    stringsAsFactors = FALSE
  )

  all_p <- c(rows$p_infit, rows$p_outfit, rows$p_gamma)
  result <- list(
    items = rows,
    side_file = side_file,
    bh_limits = c(
      fdr_5 = source_bh_critical(all_p, 0.05),
      fdr_1 = source_bh_critical(all_p, 0.01)
    ),
    fit = fit,
    extended = extended
  )
  class(result) <- c("gRm_item_fits_values", class(result))
  result
}

active_gllrm_item_fits_values <- function(fit, include_extended = TRUE) {
  context <- fit$fit$context
  state <- fit$fit
  state$context <- NULL
  bundle <- context$bundle
  fit_like <- list(
    model = "active_gllrm",
    item_gamma = state$item_gamma,
    counts = context$counts,
    context = context,
    state = state
  )
  conditional <- active_gllrm_item_conditional_moments(
    context,
    state,
    include_probabilities = TRUE,
    probability_cache = new_active_gllrm_probability_cache(context, state)
  )

  item_fit <- calculate_conditional_item_fit_values(bundle, fit_like, conditional = conditional)
  gamma_fit <- calculate_item_restscore_gamma_values(bundle, fit_like, conditional = conditional)
  component_gamma <- active_gllrm_component_restscore_values(context, state)
  extended <- if (include_extended) {
    out <- calculate_extended_item_fit_values(bundle, fit_like, conditional = conditional)
    out$component_restscore_tables <- component_gamma$table_rows
    out
  } else {
    NULL
  }

  rows <- data.frame(
    item_label = bundle$model$items$label,
    item_name = bundle$model$items$name,
    outfit = item_fit$outfit,
    outfit_sd = item_fit$outfit_sd,
    p_outfit = item_fit$p_outfit,
    outfit_fdr = item_fdr_risk(item_fit$p_outfit),
    infit = item_fit$infit,
    infit_sd = item_fit$infit_sd,
    p_infit = item_fit$p_infit,
    infit_fdr = item_fdr_risk(item_fit$p_infit),
    observed_gamma = gamma_fit$observed_gamma,
    expected_gamma = gamma_fit$expected_gamma,
    gamma_sd = gamma_fit$gamma_sd,
    p_gamma = gamma_fit$p_gamma,
    gamma_fdr = item_fdr_risk(gamma_fit$p_gamma),
    stringsAsFactors = FALSE
  )
  rows$direction <- item_fit_direction(rows)

  side_file <- data.frame(
    item = rows$item_label,
    outfit = rows$outfit,
    p_outfit = rows$p_outfit,
    infit = rows$infit,
    p_infit = rows$p_infit,
    ObsGamma = rows$observed_gamma,
    ExpGamma = rows$expected_gamma,
    p_gamma = rows$p_gamma,
    stringsAsFactors = FALSE
  )

  all_p <- c(rows$p_infit, rows$p_outfit, rows$p_gamma, component_gamma$rows$p_gamma)
  result <- list(
    items = rows,
    component_gamma = component_gamma$rows,
    side_file = side_file,
    bh_limits = c(
      fdr_5 = source_bh_critical(all_p, 0.05),
      fdr_1 = source_bh_critical(all_p, 0.01)
    ),
    fit = fit_like,
    extended = extended,
    gllrm_copy_exists = TRUE
  )
  class(result) <- c("gRm_item_fits_values", class(result))
  result
}

active_gllrm_component_restscore_values <- function(context, state) {
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  active_components <- components[lengths(components) > 1L]
  if (length(active_components) == 0L) {
    empty_rows <- data.frame(
      component = character(),
      representative_item = integer(),
      representative_item_label = character(),
      observed_gamma = numeric(),
      expected_gamma = numeric(),
      gamma_sd = numeric(),
      p_gamma = numeric(),
      stringsAsFactors = FALSE
    )
    empty_tables <- data.frame(
      component = character(),
      representative_item_label = character(),
      restscore_table = character(),
      restscore = integer(),
      component_score = integer(),
      value = numeric(),
      stringsAsFactors = FALSE
    )
    return(list(rows = empty_rows, table_rows = empty_tables))
  }

  rows <- list()
  table_rows <- list()
  for (component_items in active_components) {
    tables <- active_gllrm_component_restscore_tables(context, state, component_items)
    observed_gamma <- goodman_kruskal_gamma(tables$observed)
    fitted <- fitted_gamma_stats(tables$expected)
    expected_gamma <- fitted$gamma
    gamma_sd <- sqrt(fitted$variance)
    p_gamma <- if (gamma_sd > 0) {
      two_sided_source_normal_p((observed_gamma - expected_gamma) / gamma_sd)
    } else {
      1
    }
    component <- paste(context$items$label[component_items], collapse = "")
    representative_item <- component_items[[length(component_items)]]
    representative_item_label <- context$items$label[[representative_item]]
    rows[[length(rows) + 1L]] <- data.frame(
      component = component,
      representative_item = representative_item,
      representative_item_label = representative_item_label,
      observed_gamma = observed_gamma,
      expected_gamma = expected_gamma,
      gamma_sd = gamma_sd,
      p_gamma = p_gamma,
      stringsAsFactors = FALSE
    )
    table_rows[[length(table_rows) + 1L]] <- component_restscore_table_rows(
      component = component,
      representative_item_label = representative_item_label,
      observed_table = tables$observed,
      expected_table = tables$expected
    )
  }

  list(
    rows = type.convert(do.call(rbind, rows), as.is = TRUE),
    table_rows = type.convert(do.call(rbind, table_rows), as.is = TRUE)
  )
}

active_gllrm_component_restscore_tables <- function(context, state, component_items) {
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

active_gllrm_item_conditional_moments <- function(context,
                                                  state,
                                                  include_probabilities = FALSE,
                                                  probability_cache = NULL) {
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  probability_cache <- probability_cache %||%
    new_active_gllrm_probability_cache(context, state, components = components)
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
      group_probabilities <- active_gllrm_cached_item_probabilities(
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

active_gllrm_group_item_probabilities <- function(context,
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

#' @keywords internal
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
    item_max <- items$raw_max[[item_index]] - 1L
    item_values <- seq.int(0L, item_max)
    item_score_rows <- list()

    for (score in seq.int(1L, bundle$model$max_total_score - 1L)) {
      n_score <- score_weights[[score + 1L]]
      score_rows_mask <- observed_rows & item_scores == score
      observed_count <- tabulate(
        data[[items$name[[item_index]]]][score_rows_mask] + 1L,
        nbins = item_max + 1L
      )
      observed_total <- sum(observed_count)
      if (n_score <= 0L || observed_total <= 0L) {
        next
      }
      moments <- conditional[[item_index]][[score + 1L]]
      variance <- moments$variance
      if (variance <= 0) {
        next
      }

      probabilities <- moments$probabilities
      observed_frequency <- observed_count / observed_total
      residual <- item_values - moments$mean
      squared_residual <- residual^2
      standardized <- residual / sqrt(variance)
      squared_standardized <- squared_residual / variance

      # Source trace: the extended "Observed freqs" line prints this source
      # value with fixed zero decimals, but the source TSV/value layer stores
      # observed_frequency * ScoreDistribution without report-side rounding.
      # The mixed valid-background frequency and complete-row score weight is
      # mathematically strange but source-faithful.
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
      # Source trace: skbias15.pas::CalculateOutfit computes
      # z := (outfit-1)/StdError.
      outfit_z <- if (outfit_standard_error > 0) {
        (outfit_contribution - 1) / outfit_standard_error
      } else {
        0
      }

      infit_average <- sum(squared_residual * observed_frequency)
      infit_expected <- sum(squared_residual * probabilities)
      infit_variance <- sum(squared_residual^2 * probabilities) - infit_expected^2
      infit_ratio <- safe_ratio(infit_average, infit_expected)

      row <- data.frame(
        item_label = items$label[[item_index]],
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
        infit_ratio = infit_ratio,
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
      score_rows[[length(score_rows) + 1L]] <- row
      item_score_rows[[length(item_score_rows) + 1L]] <- row
    }

    item_scores_df <- do.call(rbind_fill, item_score_rows)
    summary_keep <- item_scores_df$n > 1
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      item_label = items$label[[item_index]],
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
    summary_rows[[length(summary_rows)]]$infit_value <- safe_ratio(
      summary_rows[[length(summary_rows)]]$infit_observed,
      summary_rows[[length(summary_rows)]]$infit_expected
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

#' @keywords internal
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

#' @keywords internal
restscore_matrix <- function(item_max, rest_max) {
  matrix(0, nrow = item_max + 1L, ncol = rest_max + 1L)
}

#' @keywords internal
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

#' @keywords internal
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
      item_max <- items$raw_max[[item_index]] - 1L
      rest_max <- max_score - item_max
      observed_table <- matrix(0, nrow = item_max + 1L, ncol = max_score + 1L)
      expected_table <- matrix(0, nrow = item_max + 1L, ncol = max_score + 1L)

      # Source trace: local adjacent tables inherit exactly one deterministic
      # boundary: low only for the 0/1 slice, high only for the top adjacent
      # slice of that item.
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
            # Source trace: Calculate_item_restscore_gamma1 rescales the
            # complete-score expected margin to the valid observed score total:
            # (Expected[cell] / ExpTotal) * ObsTotal.
            expected_table[item_score + 1L, rest_score + 1L] <-
              expected_table[item_score + 1L, rest_score + 1L] +
              (expected_counts[[item_score + 1L]] / expected_total) * observed_total
            expected_n <- expected_n + expected_table[item_score + 1L, rest_score + 1L]
          }
        }
      }

      table_rows[[length(table_rows) + 1L]] <- restscore_table_rows(
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
      p_values[[item_index]] <- gamma_p
      adjacent_gamma[[item_index]] <- data.frame(
        item_label = items$label[[item_index]],
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

#' @keywords internal
restscore_table_rows <- function(items, item_index, local_restscore, observed_table, expected_table) {
  tables <- list(observed = observed_table, expected = expected_table)
  item_rows <- lapply(names(tables), function(table_name) {
    tab <- tables[[table_name]]
    item_scores <- seq.int(0L, nrow(tab) - 1L)
    rest_scores <- seq.int(0L, ncol(tab) - 1L)
    n_cells <- length(item_scores) * length(rest_scores)
    row <- data.frame(
      item_label = rep(items$label[[item_index]], n_cells),
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

#' @keywords internal
calculate_conditional_item_fit_values <- function(bundle, fit, conditional = NULL) {
  # Source trace: skbias15.pas::Calculate_item_fits / CalculateInAndOutfits
  # prints ItemOutfit, sqrt(itemOutfitvariance), pOutfit, ItemInfit,
  # sqrt(iteminfitvariance), and pInfit after source score-window aggregation.
  items <- bundle$model$items
  data <- bundle$data
  n_items <- nrow(items)

  outfit_sum <- outfit_mean <- outfit_var <- numeric(n_items)
  infit_sum <- infit_mean <- infit_var <- infit_weight <- numeric(n_items)
  n_used <- integer(n_items)

  if (is.null(conditional)) {
    conditional <- item_conditional_moments(bundle, fit$item_gamma, include_probabilities = TRUE)
  }
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

  # Source trace: skbias12a.pas::IncompleteItemfits is called from
  # skbias14.pas when NincompleteRecs > 0. It contributes only observed items
  # in incomplete records. Conditional probabilities are recomputed for the
  # subset of observed items and the partial observed score.
  # Source trace: skbias15.pas only calls IncompleteItemFits when the DIGRAM
  # runtime has populated NincompleteRecs. The supplied example ItemFits extended
  # report has no "persons with incomplete responses were included" line, so the
  # R report does not synthesize incomplete item-fit records from bundle rows.
  incomplete_rows <- integer()
  for (row_index in incomplete_rows) {
    use_item <- vapply(items$name, function(name) data[[name]][[row_index]] >= 0L, logical(1))
    partial_score <- sum(as.integer(data[row_index, items$name[use_item], drop = FALSE]))
    for (item_index in which(use_item)) {
      item_score <- data[[items$name[[item_index]]]][[row_index]]
      moments <- item_conditional_moment_for_subset(bundle, fit$item_gamma, item_index, partial_score, use_item)
      variance <- moments$variance
      if (variance <= 0) {
        next
      }
      centered <- seq.int(0L, items$raw_max[[item_index]] - 1L) - moments$mean
      outfit_range <- centered^2 / variance
      infit_range <- centered^2

      outfit_sum[[item_index]] <- outfit_sum[[item_index]] + outfit_range[[item_score + 1L]]
      outfit_mean[[item_index]] <- outfit_mean[[item_index]] + sum(outfit_range * moments$probabilities)
      outfit_var[[item_index]] <- outfit_var[[item_index]] +
        sum(outfit_range^2 * moments$probabilities) - sum(outfit_range * moments$probabilities)^2

      infit_sum[[item_index]] <- infit_sum[[item_index]] + infit_range[[item_score + 1L]]
      infit_mean[[item_index]] <- infit_mean[[item_index]] + sum(infit_range * moments$probabilities)
      infit_var[[item_index]] <- infit_var[[item_index]] +
        sum(infit_range^2 * moments$probabilities) - sum(infit_range * moments$probabilities)^2
      infit_weight[[item_index]] <- infit_weight[[item_index]] + variance
      n_used[[item_index]] <- n_used[[item_index]] + 1L
    }
  }

  outfit <- safe_ratio(outfit_sum, n_used)
  outfit_expected <- safe_ratio(outfit_mean, n_used)
  outfit_variance <- safe_ratio(outfit_var, n_used^2)
  outfit_z <- safe_z(outfit, outfit_expected, outfit_variance)

  infit <- safe_ratio(infit_sum, infit_weight)
  infit_expected <- safe_ratio(infit_mean, infit_weight)
  infit_variance <- safe_ratio(infit_var, infit_weight^2)
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

#' @keywords internal
calculate_item_restscore_gamma_values <- function(bundle, fit, conditional = NULL) {
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

  for (item_index in seq_len(n_items)) {
    item_max <- items$raw_max[[item_index]] - 1L
    rest_max <- max_score - item_max
    observed <- matrix(0, nrow = item_max + 1L, ncol = rest_max + 1L)
    expected <- observed

    # Source trace: skbias14.pas initializes the item-restscore table with
    # Nlowscore and Nhighscore at the two deterministic extreme cells.
    observed[1L, 1L] <- low_count
    expected[1L, 1L] <- low_count
    observed[item_max + 1L, rest_max + 1L] <- high_count_source
    expected[item_max + 1L, rest_max + 1L] <- high_count_source

    for (total_score in seq.int(1L, max_score - 1L)) {
      n_score <- score_weights[[total_score + 1L]]
      if (n_score <= 0L) {
        next
      }
      moments <- conditional[[item_index]][[total_score + 1L]]
      probabilities <- moments$probabilities
      if (length(probabilities) == 0L) {
        next
      }

      score_rows <- observed_rows & item_scores == total_score
      observed_count <- tabulate(
        data[[items$name[[item_index]]]][score_rows] + 1L,
        nbins = item_max + 1L
      )
      observed_total <- sum(observed_count)
      observed_frequency <- if (observed_total > 0L) observed_count / observed_total else observed_count
      for (candidate_score in seq.int(0L, item_max)) {
        candidate_rest <- total_score - candidate_score
        if (candidate_rest >= 0L && candidate_rest <= rest_max) {
          # Source trace: skbias15.pas::CalculateMeans replaces
          # ItemMargTables(score,item,iscore) with relative frequencies and
          # stores n in item_max+1 before Calculate_item_restscore_gamma
          # multiplies both observed and expected margins by that n.
          observed[candidate_score + 1L, candidate_rest + 1L] <-
            observed[candidate_score + 1L, candidate_rest + 1L] +
            observed_frequency[[candidate_score + 1L]] * n_score
          expected[candidate_score + 1L, candidate_rest + 1L] <-
            expected[candidate_score + 1L, candidate_rest + 1L] +
            probabilities[[candidate_score + 1L]] * n_score
        }
      }
    }

    # Source trace: skbias14/skbias15 only include IncompleteRecs when the DIGRAM
    # runtime has populated NincompleteRecs; in that case the report prints
    # "records will be included during calculation of item-restscore gamma". The
    # BFI and example ItemFits runtime reports do not contain that line, so compact
    # report reproduction must not synthesize incomplete rows from raw data here.

    observed_gamma[[item_index]] <- goodman_kruskal_gamma(observed)
    fitted <- fitted_gamma_stats(expected)
    expected_gamma[[item_index]] <- fitted$gamma
    gamma_sd[[item_index]] <- sqrt(fitted$variance)
    if (gamma_sd[[item_index]] > 0) {
      z <- (observed_gamma[[item_index]] - expected_gamma[[item_index]]) / gamma_sd[[item_index]]
      p_gamma[[item_index]] <- two_sided_source_normal_p(z)
    } else {
      p_gamma[[item_index]] <- 1
    }
  }

  list(
    observed_gamma = observed_gamma,
    expected_gamma = expected_gamma,
    gamma_sd = gamma_sd,
    p_gamma = p_gamma
  )
}

#' Collect source-style incomplete response records
#'
#' @param bundle An item-parameters bundle.
#' @return A list with grouped incomplete records, counts, observed scores, and
#'   observed maximum scores.
#' @keywords internal
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
#' @param bundle An item-parameters bundle.
#' @param item_gamma Estimated item gamma matrix.
#' @param incomplete Output from [collect_source_incomplete_records()].
#' @param record_index One-based incomplete record index.
#' @return Observed score plus expected score over missing items.
#' @keywords internal
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
#' @param bundle An item-parameters bundle.
#' @param item_gamma Estimated item gamma matrix.
#' @param use_item Logical vector selecting included items.
#' @return Numeric score-generating function indexed by score plus one.
#' @keywords internal
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
#' @param bundle An item-parameters bundle.
#' @param item_gamma Estimated item gamma matrix.
#' @param incomplete Output from [collect_source_incomplete_records()].
#' @return A list of expected totals and per-item conditional probabilities.
#' @keywords internal
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

#' @keywords internal
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

#' @keywords internal
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

#' @keywords internal
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
#' @param items Item metadata data frame.
#' @param item_gamma Estimated item gamma matrix.
#' @param use_item Logical vector selecting included items.
#' @return Numeric score-generating function indexed by score plus one.
#' @keywords internal
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

#' @keywords internal
goodman_kruskal_gamma <- function(tab) {
  cd <- gamma_concordance_discordance(tab)
  denominator <- cd$concordant + cd$discordant
  if (denominator <= 0) {
    return(0)
  }
  (cd$concordant - cd$discordant) / denominator
}

#' @keywords internal
gamma_concordance_discordance <- function(tab) {
  concordant <- 0
  discordant <- 0
  for (row in seq_len(nrow(tab))) {
    for (col in seq_len(ncol(tab))) {
      count <- tab[row, col]
      if (count == 0) {
        next
      }
      if (row < nrow(tab) && col < ncol(tab)) {
        concordant <- concordant + count * sum(tab[(row + 1L):nrow(tab), (col + 1L):ncol(tab), drop = FALSE])
      }
      if (row < nrow(tab) && col > 1L) {
        discordant <- discordant + count * sum(tab[(row + 1L):nrow(tab), seq_len(col - 1L), drop = FALSE])
      }
    }
  }
  list(concordant = concordant, discordant = discordant)
}

#' @keywords internal
fitted_gamma_stats <- function(expected) {
  # Source trace: skbias15.pas::CalculateFittedGAMMA. The source computes
  # gamma = PMQ / PPQ and S1 = 16 / PPQ^4 * sum(tab[i,j] *
  # (Q * AIJ[i,j] - P * DIJ[i,j])^2).
  cd <- gamma_cell_tables(expected)
  ppq <- cd$p + cd$q
  pmq <- cd$p - cd$q
  if (ppq <= 0) {
    return(list(gamma = 0, variance = 0))
  }
  factor <- 16 / (ppq^4)
  variance <- 0
  for (row in seq_len(nrow(expected))) {
    for (col in seq_len(ncol(expected))) {
      m <- cd$q * cd$aij[row, col] - cd$p * cd$dij[row, col]
      variance <- variance + expected[row, col] * m * m
    }
  }
  list(gamma = pmq / ppq, variance = factor * variance)
}

#' @keywords internal
gamma_cell_tables <- function(tab) {
  aij <- matrix(0, nrow = nrow(tab), ncol = ncol(tab))
  dij <- matrix(0, nrow = nrow(tab), ncol = ncol(tab))
  p <- 0
  q <- 0
  for (row in seq_len(nrow(tab))) {
    for (col in seq_len(ncol(tab))) {
      for (other_row in seq_len(nrow(tab))) {
        for (other_col in seq_len(ncol(tab))) {
          concordant <- (row > other_row && col > other_col) ||
            (row < other_row && col < other_col)
          discordant <- (row < other_row && col > other_col) ||
            (row > other_row && col < other_col)
          if (concordant) {
            aij[row, col] <- aij[row, col] + tab[other_row, other_col]
          } else if (discordant) {
            dij[row, col] <- dij[row, col] + tab[other_row, other_col]
          }
        }
      }
      p <- p + tab[row, col] * aij[row, col]
      q <- q + tab[row, col] * dij[row, col]
    }
  }
  list(aij = aij, dij = dij, p = p, q = q)
}

#' @keywords internal
item_fdr_risk <- function(p_values) {
  risks <- integer(length(p_values))
  for (alpha_index in seq_along(c(0.05, 0.01, 0.001))) {
    alpha <- c(0.05, 0.01, 0.001)[[alpha_index]]
    critical <- source_bh_critical(p_values, alpha)
    risks[p_values <= critical] <- alpha_index
  }
  risks
}

#' @keywords internal
item_fit_direction <- function(rows) {
  any_flag <- rows$outfit_fdr > 0L | rows$infit_fdr > 0L | rows$gamma_fdr > 0L
  direction <- rep("", nrow(rows))
  direction[any_flag & rows$observed_gamma < rows$expected_gamma] <- "low"
  direction[any_flag & rows$observed_gamma > rows$expected_gamma] <- "high"
  direction
}

#' @keywords internal
safe_ratio <- function(numerator, denominator) {
  result <- rep(0, length(numerator))
  ok <- denominator > 0
  result[ok] <- numerator[ok] / denominator[ok]
  result
}

#' @keywords internal
safe_z <- function(observed, expected, variance) {
  z <- rep(0, length(observed))
  ok <- variance > 0
  z[ok] <- (observed[ok] - expected[ok]) / sqrt(variance[ok])
  z
}

#' @keywords internal
inf_replace <- function(x, value) {
  x[!is.finite(x)] <- value
  x
}

#' @keywords internal
two_sided_source_normal_p <- function(z) {
  vapply(z, function(value) pmin(1, 2 * source_tail_norm(abs(value), TRUE)), numeric(1))
}

#' Source fixed-field integer rounding
#'
#' Pascal's fixed-width `:0` numeric formatting rounds half values away from
#' zero for the non-negative item-fit frequencies printed by `skbias15.pas`.
#'
#' @param value Numeric vector.
#' @return Integer-like numeric vector rounded as DIGRAM prints it.
#' @keywords internal
source_print_round <- function(value) {
  sign(value) * floor(abs(value) + 0.5)
}
