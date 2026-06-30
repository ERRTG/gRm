# M2/M3 expected-table helpers.
#
# Source trace: source/PAS_skunits/SKbias2.pas::calculate_expected_IJ_table,
# calculate_expected_IX_table, calculate_expected_IJK_table,
# calculate_EIJX_table, and calculate_EIXZ_table accumulate expected counts
# from current-GLLRM score-conditioned focus probabilities, not independence
# tables. The conditional focus probabilities correspond to Calculate_Cprob1,
# Calculate_Cprob2, and Calculate_Cprob3 for one-, two-, and three-item focus
# sets. The item-exogenous-exogenous and item-exogenous-score helpers retain an
# explicit compatibility note because the preserved calculate_EIXZ_table branch
# has validation-sensitive Cprob1 loop structure in the Pascal source.

new_m2_m3_focus_probability_cache <- function(context, state) {
  list(
    context = context,
    state = state,
    components = context$ld_components_items %||% gllrm_ld_components(context)$items,
    focus = new.env(parent = emptyenv())
  )
}

m2_m3_focus_cache_key <- function(focus_items, total_score, background_values, context) {
  dif_backgrounds <- context$dif_background_indices %||% seq_along(background_values)
  background_key <- if (length(dif_backgrounds) == 0L) {
    "."
  } else {
    paste(background_values[dif_backgrounds], collapse = "\r")
  }
  paste(c(paste(as.integer(focus_items), collapse = ":"), as.integer(total_score), background_key), collapse = "\r")
}

m2_m3_focus_probabilities <- function(cache, focus_items, total_score, background_values) {
  context <- cache$context
  focus_items <- as.integer(focus_items)
  key <- m2_m3_focus_cache_key(focus_items, total_score, background_values, context)
  if (!exists(key, envir = cache$focus, inherits = FALSE)) {
    assign(
      key,
      m2_m3_calculate_focus_probabilities(
        context = context,
        state = cache$state,
        components = cache$components,
        focus_items = focus_items,
        total_score = total_score,
        background_values = background_values
      ),
      envir = cache$focus
    )
  }
  get(key, envir = cache$focus, inherits = FALSE)
}

m2_m3_calculate_focus_probabilities <- function(context,
                                                state,
                                                components,
                                                focus_items,
                                                total_score,
                                                background_values) {
  if (length(focus_items) < 1L || length(focus_items) > 3L) {
    stop("M2/M3 focus probabilities support focus size 1, 2, or 3.", call. = FALSE)
  }
  if (anyNA(focus_items) || any(focus_items < 1L | focus_items > context$n_items)) {
    stop("M2/M3 focus item index is outside the item range.", call. = FALSE)
  }
  if (anyDuplicated(focus_items)) {
    stop("M2/M3 focus items must be unique.", call. = FALSE)
  }

  dims <- as.integer(context$item_raw_max[focus_items])
  out <- array(0, dim = dims)
  if (total_score < 0L || total_score > context$max_total_score) {
    return(out)
  }

  component_gammas <- lapply(components, function(component_items) {
    gamma <- gllrm_component_gamma(context, state, component_items, background_values)
    gamma$gamma * gamma$scale
  })
  full_gamma <- m2_m3_convolve_component_gammas(component_gammas, context$max_total_score)
  denominator <- full_gamma[[total_score + 1L]]
  if (denominator <= 0) {
    return(out)
  }

  focus_component_indices <- unique(context$ld_component_of[focus_items])
  rest_component_indices <- setdiff(seq_along(components), focus_component_indices)
  rest_gamma <- m2_m3_convolve_component_gammas(component_gammas[rest_component_indices], context$max_total_score)
  focus_blocks <- m2_m3_focus_component_blocks(
    context = context,
    state = state,
    components = components,
    focus_component_indices = focus_component_indices,
    background_values = background_values
  )
  if (any(vapply(focus_blocks, nrow, integer(1L)) == 0L)) {
    return(out)
  }

  row_grid <- expand.grid(lapply(focus_blocks, function(block) seq_len(nrow(block))), KEEP.OUT.ATTRS = FALSE)
  for (grid_row in seq_len(nrow(row_grid))) {
    block_score <- 0L
    block_weight <- 1
    focus_values <- integer(length(focus_items))
    for (block_index in seq_along(focus_blocks)) {
      block <- focus_blocks[[block_index]]
      block_row <- block[row_grid[grid_row, block_index], , drop = FALSE]
      block_score <- block_score + block_row$score[[1L]]
      block_weight <- block_weight * block_row$weight[[1L]]

      component_items <- components[[focus_component_indices[[block_index]]]]
      key <- gllrm_component_key(component_items)
      config <- context$component_config_matrices[[key]][block_row$row[[1L]], ]
      matched <- match(focus_items, component_items)
      present <- !is.na(matched)
      focus_values[present] <- config[matched[present]]
    }
    if (block_score > total_score) {
      next
    }
    rest_value <- rest_gamma[[total_score - block_score + 1L]]
    if (rest_value <= 0) {
      next
    }
    out[matrix(focus_values + 1L, nrow = 1L)] <-
      out[matrix(focus_values + 1L, nrow = 1L)] + block_weight * rest_value / denominator
  }

  out
}

m2_m3_convolve_component_gammas <- function(gammas, max_total_score) {
  out <- c(1, numeric(max_total_score))
  for (gamma in gammas) {
    out <- convolve_score_vectors(out, gamma, max_total_score)
  }
  out
}

m2_m3_focus_component_blocks <- function(context,
                                         state,
                                         components,
                                         focus_component_indices,
                                         background_values) {
  lapply(focus_component_indices, function(component_index) {
    component_items <- components[[component_index]]
    key <- gllrm_component_key(component_items)
    config_matrix <- context$component_config_matrices[[key]]
    weights <- gllrm_component_config_weights_fast(
      context,
      state,
      component_items,
      config_matrix,
      background_values,
      key = key
    )
    block <- data.frame(
      row = seq_len(nrow(config_matrix)),
      score = context$component_config_scores[[key]],
      weight = weights
    )
    block[block$weight > 0, , drop = FALSE]
  })
}

m2_m3_expected_table <- function(context, state, spec, score_group_lookup, probability_cache) {
  cache <- probability_cache %||% new_m2_m3_focus_probability_cache(context, state)
  switch(
    spec$kind,
    item_item = m2_m3_expected_item_item(context, spec, NULL, cache),
    item_exogenous = m2_m3_expected_item_exogenous(context, spec, NULL, cache),
    item_score_group = m2_m3_expected_item_score_group(context, spec, score_group_lookup, cache),
    item_item_item = m2_m3_expected_item_item_item(context, spec, NULL, cache),
    item_item_exogenous = m2_m3_expected_item_item_exogenous(context, spec, NULL, cache),
    item_item_score_group = m2_m3_expected_item_item_score_group(context, spec, score_group_lookup, cache),
    item_exogenous_exogenous = m2_m3_expected_item_exogenous_exogenous(context, spec, NULL, cache),
    item_exogenous_score_group = m2_m3_expected_item_exogenous_score_group(context, spec, score_group_lookup, cache),
    stop("Unknown M2/M3 expected-table margin kind.", call. = FALSE)
  )
}

m2_m3_expected_item_item <- function(context, spec, score_groups, cache) {
  m2_m3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = c(spec$item1, spec$item2)
  )
}

m2_m3_expected_item_exogenous <- function(context, spec, score_groups, cache) {
  m2_m3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = spec$item,
    extra_dims = context$background_raw_max[[spec$exogenous]],
    extra_value = function(group, background_values) background_values[[spec$exogenous]]
  )
}

m2_m3_expected_item_score_group <- function(context, spec, score_group_lookup, cache) {
  m2_m3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = spec$item,
    extra_dims = m2_m3_score_group_count(score_group_lookup),
    extra_value = function(group, background_values) {
      global_homogeneity_lookup_score(score_group_lookup, group$score[[1L]])
    }
  )
}

m2_m3_expected_item_item_item <- function(context, spec, score_groups, cache) {
  m2_m3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = c(spec$item1, spec$item2, spec$item3)
  )
}

m2_m3_expected_item_item_exogenous <- function(context, spec, score_groups, cache) {
  m2_m3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = c(spec$item1, spec$item2),
    extra_dims = context$background_raw_max[[spec$exogenous]],
    extra_value = function(group, background_values) background_values[[spec$exogenous]]
  )
}

m2_m3_expected_item_item_score_group <- function(context, spec, score_group_lookup, cache) {
  m2_m3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = c(spec$item1, spec$item2),
    extra_dims = m2_m3_score_group_count(score_group_lookup),
    extra_value = function(group, background_values) {
      global_homogeneity_lookup_score(score_group_lookup, group$score[[1L]])
    }
  )
}

m2_m3_expected_item_exogenous_exogenous <- function(context, spec, score_groups, cache) {
  # Compatibility-sensitive EIXZ branch: source/PAS_skunits/SKbias2.pas::
  # calculate_EIXZ_table uses Calculate_Cprob1 for the single focus item.
  m2_m3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = spec$item,
    extra_dims = c(context$background_raw_max[[spec$exogenous1]], context$background_raw_max[[spec$exogenous2]]),
    extra_value = function(group, background_values) {
      c(background_values[[spec$exogenous1]], background_values[[spec$exogenous2]])
    }
  )
}

m2_m3_expected_item_exogenous_score_group <- function(context, spec, score_group_lookup, cache) {
  # Compatibility-sensitive EIXZ branch: source/PAS_skunits/SKbias2.pas::
  # calculate_EIXZ_table also covers the score-group pseudo-variable case.
  m2_m3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = spec$item,
    extra_dims = c(context$background_raw_max[[spec$exogenous]], m2_m3_score_group_count(score_group_lookup)),
    extra_value = function(group, background_values) {
      c(
        background_values[[spec$exogenous]],
        global_homogeneity_lookup_score(score_group_lookup, group$score[[1L]])
      )
    }
  )
}

m2_m3_expected_focus_table <- function(context,
                                       cache,
                                       focus_items,
                                       extra_dims = integer(),
                                       extra_value = NULL) {
  dims <- c(as.integer(context$item_raw_max[focus_items]), as.integer(extra_dims))
  out <- array(0, dim = dims)
  groups <- m2_m3_source_score_background_groups(context)
  for (group_index in seq_len(nrow(groups))) {
    group <- groups[group_index, , drop = FALSE]
    background_values <- gllrm_group_background_values(context, group)
    suffix <- if (is.null(extra_value)) {
      integer()
    } else {
      as.integer(extra_value(group, background_values))
    }
    if (length(suffix) > 0L && any(is.na(suffix))) {
      next
    }
    probabilities <- m2_m3_focus_probabilities(
      cache,
      focus_items = focus_items,
      total_score = group$score[[1L]],
      background_values = background_values
    )
    out <- m2_m3_add_expected_probabilities(
      out = out,
      probabilities = probabilities,
      suffix = suffix,
      count = group$count[[1L]]
    )
  }
  out
}

m2_m3_source_score_background_groups <- function(context) {
  # Source trace: source/PAS_skunits/SKbias2.pas::calculate_expected_IX_table,
  # calculate_EIJX_table, and calculate_EIXZ_table loop over source records,
  # attach the observed exogenous/score-group value, and only use DIF
  # background values for the conditional probabilities. The fitting
  # score_exo_groups object is a smaller sufficient statistic and is not enough
  # for CM2/CM3 expected margins that mention non-DIF exogenous variables.
  rows <- context$valid_rows
  if (!length(rows)) {
    out <- data.frame(score = integer(), count = integer())
    for (name in context$backgrounds$name) {
      out[[name]] <- integer()
    }
    return(out)
  }

  groups <- list()
  index <- new.env(parent = emptyenv())
  for (row in rows) {
    values <- if (context$n_backgrounds > 0L) {
      as.integer(context$background_matrix[row, ])
    } else {
      integer()
    }
    if (length(values) && any(is.na(values) | values < 1L)) {
      next
    }
    key <- paste(c(context$score[[row]], values), collapse = "\r")
    pos <- index[[key]]
    if (is.null(pos)) {
      pos <- length(groups) + 1L
      index[[key]] <- pos
      groups[[pos]] <- c(score = context$score[[row]], count = 0L, values)
    }
    groups[[pos]][["count"]] <- groups[[pos]][["count"]] + 1L
  }
  if (!length(groups)) {
    out <- data.frame(score = integer(), count = integer())
    for (name in context$backgrounds$name) {
      out[[name]] <- integer()
    }
    return(out)
  }

  mat <- do.call(rbind, groups)
  out <- data.frame(score = as.integer(mat[, "score"]), count = as.integer(mat[, "count"]))
  if (context$n_backgrounds > 0L) {
    for (background_index in seq_len(context$n_backgrounds)) {
      out[[context$backgrounds$name[[background_index]]]] <-
        as.integer(mat[, 2L + background_index])
    }
  }
  out
}

m2_m3_add_expected_probabilities <- function(out, probabilities, suffix, count) {
  if (length(suffix) == 0L) {
    return(out + count * probabilities)
  }

  values <- as.numeric(probabilities)
  positive <- values != 0
  if (!any(positive)) {
    return(out)
  }
  prefix_indices <- arrayInd(seq_along(probabilities), dim(probabilities))
  for (cell_index in which(positive)) {
    target <- c(prefix_indices[cell_index, ], suffix)
    out[matrix(target, nrow = 1L)] <- out[matrix(target, nrow = 1L)] + count * values[[cell_index]]
  }
  out
}
