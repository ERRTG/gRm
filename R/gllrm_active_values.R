gllrm_active_values <- function(active_fit, spec) {
  context <- active_fit$context
  state <- active_fit$state
  item_values <- gllrm_active_item_values(context, state)
  item_values$log_likelihood <- state$log_likelihood
  item_values$active_context <- context
  item_values$ld_parameters <- state$ld_parameters
  item_values$dif_parameters <- state$dif_parameters
  item_values$expected_items <- state$expected_items
  item_values$expected_ld <- state$expected_ld
  item_values$expected_dif <- state$expected_dif
  item_values$update_items <- state$update_items
  item_values$update_ld <- state$update_ld
  item_values$update_dif <- state$update_dif
  class(item_values) <- c("gRm_active_gllrm_values", class(item_values))
  item_values
}

gllrm_active_item_values <- function(context, state) {
  items <- context$items
  max_category_count <- ncol(state$item_gamma)
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
    gamma_values <- as.numeric(state$item_gamma[item_index, seq_len(items$raw_max[[item_index]])])
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

  list(
    n_step = state$n_step,
    delta = state$report_delta,
    log_likelihood = state$log_likelihood,
    likelihood_n = context$counts$n_valid,
    input_stats = item_parameters_input_stats(context$bundle),
    item_names = item_names,
    item_labels = items$label_code,
    background_labels = context$backgrounds$label_code,
    background_max = context$backgrounds$raw_max,
    n_parameters = calculate_active_gllrm_n_parameters(context),
    observed_score_range = calculate_observed_score_range(context$counts$item_counts),
    item_gamma = state$item_gamma,
    thresholds = thresholds,
    locations = locations,
    ice = ice,
    mice = mice,
    ice_item_effect = item_effect,
    mice_item_effect = mice_item_effect,
    ld_parameter_tables = gllrm_active_ld_parameter_values(context, state),
    dif_parameter_tables = gllrm_active_dif_parameter_values(context, state),
    item_statistics = data.frame(
      item = item_names,
      location = unname(locations),
      midpoint = unname(midpoints),
      target = unname(targets),
      info_at_target = unname(info_at_target),
      info_per_step = unname(info_per_step),
      stringsAsFactors = FALSE,
      row.names = item_names
    ),
    ice_fields = NULL
  )
}

gllrm_active_detail_tables <- function(values) {
  context <- values$active_context
  list(
    ld_parameters = gllrm_ld_parameter_table(context, values$ld_parameters),
    dif_parameters = gllrm_dif_parameter_table(context, values$dif_parameters),
    expected_items = gllrm_item_margin_table(context, values$expected_items, "expected"),
    expected_ld = gllrm_ld_margin_table(context, values$expected_ld, "expected"),
    expected_dif = gllrm_dif_margin_table(context, values$expected_dif, "expected"),
    update_items = gllrm_item_margin_table(context, values$update_items, "update"),
    update_ld = gllrm_ld_margin_table(context, values$update_ld, "update"),
    update_dif = gllrm_dif_margin_table(context, values$update_dif, "update")
  )
}

gllrm_ld_parameter_table <- function(context, parameters) {
  rows <- list()
  for (ld_index in seq_along(parameters)) {
    spec <- context$ld_specs[[ld_index]]
    gamma <- parameters[[ld_index]]
    for (i in seq_len(nrow(gamma))) {
      for (j in seq_len(ncol(gamma))) {
        rows[[length(rows) + 1L]] <- data.frame(
          term = spec$term,
          item1 = spec$item1_name,
          item2 = spec$item2_name,
          score1 = i - 1L,
          score2 = j - 1L,
          gamma = gamma[i, j],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

gllrm_dif_parameter_table <- function(context, parameters) {
  rows <- list()
  for (dif_index in seq_along(parameters)) {
    spec <- context$dif_specs[[dif_index]]
    gamma <- parameters[[dif_index]]
    for (i in seq_len(nrow(gamma))) {
      for (j in seq_len(ncol(gamma))) {
        rows[[length(rows) + 1L]] <- data.frame(
          term = spec$term,
          item = spec$item_name,
          exogenous = spec$background_name,
          score = i - 1L,
          value = j,
          gamma = gamma[i, j],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

gllrm_item_margin_table <- function(context, mat, value_name) {
  rows <- list()
  for (item in seq_len(context$n_items)) {
    for (score in context$item_score_values[[item]]) {
      rows[[length(rows) + 1L]] <- data.frame(
        item = context$items$name[[item]],
        score = score,
        value = mat[item, score + 1L],
        stringsAsFactors = FALSE
      )
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  names(out)[names(out) == "value"] <- value_name
  out
}

gllrm_ld_margin_table <- function(context, margins, value_name) {
  rows <- list()
  for (ld_index in seq_along(margins)) {
    spec <- context$ld_specs[[ld_index]]
    margin <- margins[[ld_index]]
    for (i in seq_len(nrow(margin))) {
      for (j in seq_len(ncol(margin))) {
        rows[[length(rows) + 1L]] <- data.frame(
          term = spec$term,
          item1 = spec$item1_name,
          item2 = spec$item2_name,
          score1 = i - 1L,
          score2 = j - 1L,
          value = margin[i, j],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  names(out)[names(out) == "value"] <- value_name
  out
}

gllrm_dif_margin_table <- function(context, margins, value_name) {
  rows <- list()
  for (dif_index in seq_along(margins)) {
    spec <- context$dif_specs[[dif_index]]
    margin <- margins[[dif_index]]
    for (i in seq_len(nrow(margin))) {
      for (j in seq_len(ncol(margin))) {
        rows[[length(rows) + 1L]] <- data.frame(
          term = spec$term,
          item = spec$item_name,
          exogenous = spec$background_name,
          score = i - 1L,
          value = j,
          margin = margin[i, j],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  names(out)[names(out) == "margin"] <- value_name
  out
}

source_gamma_from_table <- function(table) {
  p <- 0
  q <- 0
  x_dim <- nrow(table)
  y_dim <- ncol(table)
  for (x in seq_len(x_dim)) {
    for (y in seq_len(y_dim)) {
      aij <- 0
      dij <- 0
      for (k in seq_len(x_dim)) {
        for (l in seq_len(y_dim)) {
          concordant <- (x > k && y > l) || (x < k && y < l)
          discordant <- (x < k && y > l) || (x > k && y < l)
          if (concordant) {
            aij <- aij + table[k, l]
          } else if (discordant) {
            dij <- dij + table[k, l]
          }
        }
      }
      p <- p + table[x, y] * aij
      q <- q + table[x, y] * dij
    }
  }
  ppq <- p + q
  if (ppq > 0) (p - q) / ppq else 0
}

standardize_parameter_table_source <- function(table, row_margins, col_margins) {
  standardized <- table
  for (step in seq_len(30L)) {
    for (row in seq_len(nrow(standardized))) {
      total <- sum(standardized[row, ])
      if (total > 0) {
        standardized[row, ] <- standardized[row, ] * row_margins[[row]] / total
      }
    }
    for (col in seq_len(ncol(standardized))) {
      total <- sum(standardized[, col])
      if (total > 0) {
        standardized[, col] <- standardized[, col] * col_margins[[col]] / total
      }
    }
  }
  list(
    table = standardized,
    gamma = source_gamma_from_table(standardized),
    odds_ratio = source_odds_ratio_from_gamma(source_gamma_from_table(standardized))
  )
}

source_odds_ratio_from_gamma <- function(gamma) {
  if (gamma == 1) {
    999.99
  } else {
    (1 + gamma) / (1 - gamma)
  }
}

gllrm_active_ld_parameter_values <- function(context, state) {
  lapply(seq_along(context$ld_specs), function(ld_index) {
    spec <- context$ld_specs[[ld_index]]
    table <- state$ld_parameters[[ld_index]]
    row_margins <- context$counts$item_counts[spec$item1, seq_len(nrow(table))]
    col_margins <- context$counts$item_counts[spec$item2, seq_len(ncol(table))]
    standardized <- standardize_parameter_table_source(table, row_margins, col_margins)
    list(
      spec = spec,
      item1_code = context$items$label_code[[spec$item1]],
      item2_code = context$items$label_code[[spec$item2]],
      raw = table,
      standardized = standardized$table,
      gamma = standardized$gamma,
      odds_ratio = standardized$odds_ratio
    )
  })
}

gllrm_active_dif_parameter_values <- function(context, state) {
  lapply(seq_along(context$dif_specs), function(dif_index) {
    spec <- context$dif_specs[[dif_index]]
    table <- state$dif_parameters[[dif_index]]
    row_margins <- context$counts$item_counts[spec$item, seq_len(nrow(table))]
    col_margins <- colSums(context$observed_dif[[dif_index]])
    standardized <- standardize_parameter_table_source(table, row_margins, col_margins)
    list(
      spec = spec,
      item_code = context$items$label_code[[spec$item]],
      background_code = context$backgrounds$label_code[[spec$background]],
      raw = table,
      standardized = standardized$table,
      gamma = standardized$gamma,
      odds_ratio = standardized$odds_ratio
    )
  })
}

calculate_active_gllrm_n_parameters <- function(context) {
  base <- calculate_source_n_parameters(context$counts$item_counts)
  ld_df <- sum(vapply(context$observed_ld, source_observed_margin_n_parameters, integer(1L)))
  dif_df <- sum(vapply(seq_along(context$dif_specs), function(dif_index) {
    source_observed_margin_n_parameters(context$observed_dif[[dif_index]])
  }, integer(1L)))
  base + ld_df + dif_df
}

source_observed_margin_n_parameters <- function(observed) {
  row_nonzero <- rowSums(observed > 0)
  col_nonzero <- colSums(observed > 0)
  if (!any(row_nonzero > 0) || !any(col_nonzero > 0)) {
    return(0L)
  }
  row_ref <- which.max(row_nonzero)
  col_ref <- which.max(col_nonzero)
  sum(observed[-row_ref, -col_ref, drop = FALSE] > 0)
}
