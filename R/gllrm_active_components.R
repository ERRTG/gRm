gllrm_ld_components <- function(context) {
  parent <- seq_len(context$n_items)
  find_root <- function(item) {
    while (parent[[item]] != item) {
      parent[[item]] <<- parent[[parent[[item]]]]
      item <- parent[[item]]
    }
    item
  }
  if (length(context$ld_specs) > 0L) {
    for (spec in context$ld_specs) {
      parent[[find_root(spec$item2)]] <- find_root(spec$item1)
    }
  }

  root_to_component <- integer(context$n_items)
  components <- list()
  component_of <- integer(context$n_items)
  for (item in seq_len(context$n_items)) {
    root <- find_root(item)
    if (root_to_component[[root]] == 0L) {
      root_to_component[[root]] <- length(components) + 1L
      components[[length(components) + 1L]] <- integer()
    }
    component <- root_to_component[[root]]
    component_of[[item]] <- component
    components[[component]] <- c(components[[component]], item)
  }
  list(items = components, component_of = component_of)
}

initialize_gllrm_active_state <- function(context) {
  item_gamma <- initial_gllrm_active_item_gamma(context)
  ld_parameters <- lapply(context$ld_specs, function(spec) {
    rows <- context$item_raw_max[[spec$item1]]
    cols <- context$item_raw_max[[spec$item2]]
    matrix(
      0,
      nrow = rows,
      ncol = cols,
      dimnames = list(as.character(seq.int(0L, rows - 1L)), as.character(seq.int(0L, cols - 1L)))
    )
  })
  for (ld_index in seq_along(ld_parameters)) {
    ld_parameters[[ld_index]][context$observed_ld[[ld_index]] > 0] <- 1
  }
  dif_parameters <- lapply(context$dif_specs, function(spec) {
    rows <- context$item_raw_max[[spec$item]]
    cols <- context$background_raw_max[[spec$background]]
    matrix(
      0,
      nrow = rows,
      ncol = cols,
      dimnames = list(as.character(seq.int(0L, rows - 1L)), as.character(seq_len(cols)))
    )
  })
  for (dif_index in seq_along(dif_parameters)) {
    dif_parameters[[dif_index]][context$observed_dif[[dif_index]] > 0] <- 1
  }
  list(
    item_gamma = item_gamma,
    ld_parameters = ld_parameters,
    dif_parameters = dif_parameters,
    expected_items = item_gamma * 0,
    expected_ld = lapply(ld_parameters, function(x) x * 0),
    expected_dif = lapply(dif_parameters, function(x) x * 0),
    update_items = item_gamma * 0 + 1,
    update_ld = lapply(ld_parameters, function(x) x * 0),
    update_dif = lapply(dif_parameters, function(x) x * 0),
    delta = 0,
    report_delta = 0,
    n_step = 0L,
    converged = FALSE,
    log_likelihood = NA_real_
  )
}

initial_gllrm_active_item_gamma <- function(context) {
  # Source trace: source/GLLRM.txt::Estimate_GLLRM initializes item
  # parameters only for item scores observed in the active estimation range.
  gamma <- initial_item_gamma(context$bundle)
  gamma[,] <- 0
  observed <- context$counts$item_counts
  for (item_index in seq_len(context$n_items)) {
    cols <- seq_len(context$item_raw_max[[item_index]])
    gamma[item_index, cols] <- as.numeric(observed[item_index, cols] > 0)
  }
  gamma
}

gllrm_active_dif_map <- function(context) {
  out <- list()
  if (length(context$dif_specs) == 0L) {
    return(out)
  }
  for (dif_index in seq_along(context$dif_specs)) {
    spec <- context$dif_specs[[dif_index]]
    out[[paste(spec$item, spec$background, sep = ":")]] <- dif_index
  }
  out
}

gllrm_active_dif_by_item <- function(context) {
  out <- vector("list", context$n_items)
  for (item in seq_len(context$n_items)) {
    out[[item]] <- data.frame(
      background = integer(),
      dif_index = integer()
    )
  }
  for (dif_index in seq_along(context$dif_specs)) {
    spec <- context$dif_specs[[dif_index]]
    out[[spec$item]] <- rbind(
      out[[spec$item]],
      data.frame(
        background = spec$background,
        dif_index = dif_index
      )
    )
  }
  out
}

gllrm_component_key <- function(component_items) {
  paste(as.integer(component_items), collapse = ":")
}

gllrm_component_configurations <- function(context, component_items) {
  key <- gllrm_component_key(component_items)
  if (!is.null(context$component_configurations) &&
      !is.null(context$component_configurations[[key]])) {
    return(context$component_configurations[[key]])
  }
  gllrm_build_component_configurations(context, component_items)
}

gllrm_build_component_configurations <- function(context, component_items) {
  grids <- lapply(component_items, function(item) context$item_score_values[[item]])
  if (length(grids) == 1L) {
    out <- data.frame(value = grids[[1L]])
    names(out) <- as.character(component_items[[1L]])
    return(out)
  }
  out <- expand.grid(grids, KEEP.OUT.ATTRS = FALSE)
  names(out) <- as.character(component_items)
  out
}

gllrm_component_config_weight <- function(context, state, component_items, item_values, background_values) {
  weight <- 1
  dif_map <- context$dif_map %||% gllrm_active_dif_map(context)
  for (item in component_items) {
    item_score <- item_values[[as.character(item)]]
    weight <- weight * state$item_gamma[item, item_score + 1L]
    if (context$n_backgrounds > 0L) {
      for (background in seq_len(context$n_backgrounds)) {
        dif_index <- dif_map[[paste(item, background, sep = ":")]]
        if (!is.null(dif_index)) {
          value <- background_values[[background]]
          weight <- weight * state$dif_parameters[[dif_index]][item_score + 1L, value]
        }
      }
    }
  }
  if (length(context$ld_specs) > 0L) {
    key <- gllrm_component_key(component_items)
    component_lookup <- context$component_lookup[[key]]
    if (is.null(component_lookup)) {
      component_lookup <- logical(context$n_items)
      component_lookup[component_items] <- TRUE
    }
    ld_indices <- context$component_ld_indices[[key]] %||% seq_along(context$ld_specs)
    for (ld_index in ld_indices) {
      spec <- context$ld_specs[[ld_index]]
      if (isTRUE(component_lookup[[spec$item1]]) && isTRUE(component_lookup[[spec$item2]])) {
        weight <- weight * state$ld_parameters[[ld_index]][
          item_values[[as.character(spec$item1)]] + 1L,
          item_values[[as.character(spec$item2)]] + 1L
        ]
      }
    }
  }
  weight
}

gllrm_component_config_weight_fast <- function(context,
                                               state,
                                               component_items,
                                               item_values,
                                               background_values,
                                               key = gllrm_component_key(component_items)) {
  weight <- 1
  for (local_index in seq_along(component_items)) {
    item <- component_items[[local_index]]
    item_score <- item_values[[local_index]]
    weight <- weight * state$item_gamma[item, item_score + 1L]
    dif_rows <- context$dif_by_item_matrices[[item]] %||% context$dif_by_item[[item]]
    if (nrow(dif_rows) > 0L) {
      for (dif_row in seq_len(nrow(dif_rows))) {
        background <- dif_rows[dif_row, 1L]
        dif_index <- dif_rows[dif_row, 2L]
        value <- background_values[[background]]
        weight <- weight * state$dif_parameters[[dif_index]][item_score + 1L, value]
      }
    }
  }

  ld_local <- context$component_ld_local_matrices[[key]] %||% context$component_ld_local_indices[[key]]
  if (!is.null(ld_local) && nrow(ld_local) > 0L) {
    for (ld_row in seq_len(nrow(ld_local))) {
      weight <- weight * state$ld_parameters[[ld_local[ld_row, 1L]]][
        item_values[[ld_local[ld_row, 2L]]] + 1L,
        item_values[[ld_local[ld_row, 3L]]] + 1L
      ]
    }
  }
  weight
}

gllrm_component_config_weights_fast <- function(context,
                                                state,
                                                component_items,
                                                config_matrix,
                                                background_values,
                                                key = gllrm_component_key(component_items)) {
  weights <- rep(1, nrow(config_matrix))
  for (local_index in seq_along(component_items)) {
    item <- component_items[[local_index]]
    item_scores <- config_matrix[, local_index] + 1L
    weights <- weights * state$item_gamma[item, item_scores]
    dif_rows <- context$dif_by_item_matrices[[item]] %||% context$dif_by_item[[item]]
    if (nrow(dif_rows) > 0L) {
      for (dif_row in seq_len(nrow(dif_rows))) {
        background <- dif_rows[dif_row, 1L]
        dif_index <- dif_rows[dif_row, 2L]
        value <- background_values[[background]]
        weights <- weights * state$dif_parameters[[dif_index]][item_scores, value]
      }
    }
  }

  ld_local <- context$component_ld_local_matrices[[key]] %||% context$component_ld_local_indices[[key]]
  if (!is.null(ld_local) && nrow(ld_local) > 0L) {
    for (ld_row in seq_len(nrow(ld_local))) {
      item1_scores <- config_matrix[, ld_local[ld_row, 2L]] + 1L
      item2_scores <- config_matrix[, ld_local[ld_row, 3L]] + 1L
      weights <- weights * state$ld_parameters[[ld_local[ld_row, 1L]]][
        cbind(item1_scores, item2_scores)
      ]
    }
  }
  weights
}

gllrm_component_gamma <- function(context, state, component_items, background_values) {
  gamma <- numeric(context$max_total_score + 1L)
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
  positive <- weights > 0
  if (any(positive)) {
    gamma_sum <- rowsum(weights[positive], config_scores[positive] + 1L, reorder = FALSE)
    gamma[as.integer(rownames(gamma_sum))] <- as.numeric(gamma_sum[, 1L])
  }
  out <- gllrm_normalize_component_gamma(gamma)
  out$config_weights <- weights
  out
}

gllrm_normalize_component_gamma <- function(gamma) {
  scale <- max(gamma)
  if (scale <= 0) {
    scale <- 1
  }
  list(gamma = gamma / scale, scale = scale)
}
