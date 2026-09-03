# CM2/CM3 expected-table helpers.
#
# Source trace: source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table,
# calculate_expected_IX_table, calculate_expected_IJK_table,
# calculate_EIJX_table, and calculate_EIXZ_table accumulate expected counts
# from current-GLLRM score-conditioned focus probabilities, not independence
# tables. The conditional focus probabilities correspond to Calculate_Cprob1,
# Calculate_Cprob2, and Calculate_Cprob3 for one-, two-, and three-item focus
# sets. The item-exogenous-exogenous and item-exogenous-score helpers retain an
# explicit compatibility note because the preserved calculate_EIXZ_table branch
# has validation-sensitive Cprob1 loop structure in the Pascal source.

#' Internal new cm2 cm3 focus probability cache helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
new_cm2_cm3_focus_probability_cache <- function(context, state) {
  list(
    context = context,
    state = state,
    components = context$ld_components_items %||% gllrm_ld_components(context)$items,
    focus = new.env(parent = emptyenv()),
    background = new.env(parent = emptyenv())
  )
}

#' Cache background-conditioned component state for CM2/CM3 probabilities
#'
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_Cprob1`,
#' `calculate_Cprob2`, and `calculate_Cprob3`. For one fitted state and one DIF
#' background cell, the component gammas and configuration weights are
#' invariant across calls for different scores and focus sets. Caching those
#' values changes neither their source traversal nor floating-point order.
#' @param cache CM2/CM3 focus probability cache.
#' @param background_values One-based exogenous category values.
#' @return Background-conditioned component gammas, full gamma, configuration
#'   blocks, and a cache for focus-complement convolutions.
#' @keywords internal
#' @noRd
cm2_cm3_focus_background_state <- function(cache, background_values) {
  context <- cache$context
  dif_backgrounds <- context$dif_background_indices %||% seq_along(background_values)
  key <- if (length(dif_backgrounds)) {
    paste(background_values[dif_backgrounds], collapse = "\r")
  } else {
    "."
  }
  if (!exists(key, envir = cache$background, inherits = FALSE)) {
    # SKbias2.Calculate_Cprob1/2/3 form component-score gammas from the
    # current item, DIF, and LD weights before combining focus and complement
    # components. Keep that component order unchanged here.
    component_gammas <- lapply(cache$components, function(component_items) {
      gamma <- gllrm_component_gamma(
        context,
        cache$state,
        component_items,
        background_values
      )
      gamma$gamma * gamma$scale
    })
    blocks <- lapply(seq_along(cache$components), function(component_index) {
      component_items <- cache$components[[component_index]]
      component_key <- context$component_keys[[component_index]]
      configurations <- context$component_config_matrices[[component_key]]
      weights <- gllrm_component_config_weights_fast(
        context,
        cache$state,
        component_items,
        configurations,
        background_values,
        key = component_key
      )
      positive <- which(weights > 0)
      list(
        items = component_items,
        row = positive,
        score = context$component_config_scores[[component_key]][positive],
        weight = weights[positive],
        configurations = configurations[positive, , drop = FALSE]
      )
    })
    assign(
      key,
      list(
        component_gammas = component_gammas,
        full_gamma = cm2_cm3_convolve_component_gammas(
          component_gammas,
          context$max_total_score
        ),
        blocks = blocks,
        rest_gamma = new.env(parent = emptyenv()),
        focus_plan = new.env(parent = emptyenv())
      ),
      envir = cache$background
    )
  }
  get(key, envir = cache$background, inherits = FALSE)
}

#' Precompute one source-ordered CM2/CM3 focus traversal plan
#'
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_Cprob1`,
#' `calculate_Cprob2`, and `calculate_Cprob3`.
#'
#' The configuration scores, focus cells, and component-weight products are
#' invariant across total-score calls for one fitted state and DIF background
#' cell. The plan retains `expand.grid`'s original row order and multiplies
#' component weights in component order, exactly as the scalar source-shaped
#' loop did. Only integer addressing work is moved out of repeated score calls.
#' @param background_state Cached background-conditioned component state.
#' @param context Prepared GLLRM context.
#' @param focus_items Ordered item indices in the focus set.
#' @return Source-ordered component indices, score and weight vectors, and
#'   one-based linear focus-array cells.
#' @keywords internal
#' @noRd
cm2_cm3_focus_traversal_plan <- function(background_state, context, focus_items) {
  key <- paste(as.integer(focus_items), collapse = ":")
  if (!exists(key, envir = background_state$focus_plan, inherits = FALSE)) {
    focus_component_indices <- unique(context$ld_component_of[focus_items])
    focus_blocks <- lapply(focus_component_indices, function(component_index) {
      block <- background_state$blocks[[component_index]]
      matched <- match(focus_items, block$items)
      present <- which(!is.na(matched))
      list(
        row = block$row,
        score = block$score,
        weight = block$weight,
        focus_positions = present,
        focus_values = block$configurations[, matched[present], drop = FALSE]
      )
    })
    if (any(vapply(focus_blocks, function(block) length(block$row), integer(1L)) == 0L)) {
      plan <- NULL
    } else {
      row_grid <- as.matrix(expand.grid(
        lapply(focus_blocks, function(block) seq_along(block$row)),
        KEEP.OUT.ATTRS = FALSE
      ))
      n_combinations <- nrow(row_grid)
      block_score <- integer(n_combinations)
      block_weight <- rep(1, n_combinations)
      focus_values <- matrix(0L, nrow = n_combinations, ncol = length(focus_items))
      for (block_index in seq_along(focus_blocks)) {
        # Calculate_Cprob1/2/3 multiply the participating component weights in
        # component order. Precompute that identical product without changing
        # the multiplication sequence.
        block <- focus_blocks[[block_index]]
        selected <- row_grid[, block_index]
        block_score <- block_score + block$score[selected]
        block_weight <- block_weight * block$weight[selected]
        focus_values[, block$focus_positions] <-
          block$focus_values[selected, , drop = FALSE]
      }
      dims <- as.integer(context$item_raw_max[focus_items])
      focus_cell <- rep(1L, n_combinations)
      multiplier <- 1L
      for (focus_index in seq_along(focus_items)) {
        focus_cell <- focus_cell + focus_values[, focus_index] * multiplier
        multiplier <- multiplier * dims[[focus_index]]
      }
      plan <- list(
        focus_component_indices = focus_component_indices,
        block_score = block_score,
        block_weight = block_weight,
        focus_cell = focus_cell
      )
    }
    assign(key, list(value = plan), envir = background_state$focus_plan)
  }
  get(key, envir = background_state$focus_plan, inherits = FALSE)$value
}

#' Internal cm2 cm3 focus cache key helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param focus_items Internal `focus_items` value used by this helper.
#' @param total_score Internal `total_score` value used by this helper.
#' @param background_values Internal `background_values` value used by this helper.
#' @param context Prepared GLLRM computation context.
#' @return The internal `cm2_cm3_focus_cache_key()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_focus_cache_key <- function(focus_items, total_score, background_values, context) {
  dif_backgrounds <- context$dif_background_indices %||% seq_along(background_values)
  background_key <- if (length(dif_backgrounds) == 0L) {
    "."
  } else {
    paste(background_values[dif_backgrounds], collapse = "\r")
  }
  paste(c(paste(as.integer(focus_items), collapse = ":"), as.integer(total_score), background_key), collapse = "\r")
}

#' Internal cm2 cm3 focus probabilities helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param cache Computation cache.
#' @param focus_items Internal `focus_items` value used by this helper.
#' @param total_score Internal `total_score` value used by this helper.
#' @param background_values Internal `background_values` value used by this helper.
#' @return The internal `cm2_cm3_focus_probabilities()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_focus_probabilities <- function(cache, focus_items, total_score, background_values) {
  context <- cache$context
  focus_items <- as.integer(focus_items)
  key <- cm2_cm3_focus_cache_key(focus_items, total_score, background_values, context)
  if (!exists(key, envir = cache$focus, inherits = FALSE)) {
    background_state <- cm2_cm3_focus_background_state(cache, background_values)
    assign(
      key,
      cm2_cm3_calculate_focus_probabilities(
        context = context,
        state = cache$state,
        components = cache$components,
        focus_items = focus_items,
        total_score = total_score,
        background_values = background_values,
        background_state = background_state
      ),
      envir = cache$focus
    )
  }
  get(key, envir = cache$focus, inherits = FALSE)
}

#' Internal cm2 cm3 calculate focus probabilities helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_Cprob1`,
#' `calculate_Cprob2`, and `calculate_Cprob3`, as called by the source expected-
#' table routines.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param components Internal `components` value used by this helper.
#' @param focus_items Internal `focus_items` value used by this helper.
#' @param total_score Internal `total_score` value used by this helper.
#' @param background_values Internal `background_values` value used by this helper.
#' @param background_state Optional cached component state for the supplied
#'   background values.
#' @return The internal `cm2_cm3_calculate_focus_probabilities()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_calculate_focus_probabilities <- function(context,
                                                state,
                                                components,
                                                focus_items,
                                                total_score,
                                                background_values,
                                                background_state = NULL) {
  if (length(focus_items) < 1L || length(focus_items) > 3L) {
    stop("CM2/CM3 focus probabilities support focus size 1, 2, or 3.", call. = FALSE)
  }
  if (anyNA(focus_items) || any(focus_items < 1L | focus_items > context$n_items)) {
    stop("CM2/CM3 focus item index is outside the item range.", call. = FALSE)
  }
  if (anyDuplicated(focus_items)) {
    stop("CM2/CM3 focus items must be unique.", call. = FALSE)
  }

  dims <- as.integer(context$item_raw_max[focus_items])
  out <- array(0, dim = dims)
  if (total_score < 0L || total_score > context$max_total_score) {
    return(out)
  }

  if (is.null(background_state)) {
    temporary_cache <- list(
      context = context,
      state = state,
      components = components,
      background = new.env(parent = emptyenv())
    )
    background_state <- cm2_cm3_focus_background_state(
      temporary_cache,
      background_values
    )
  }
  denominator <- background_state$full_gamma[[total_score + 1L]]
  if (denominator <= 0) {
    return(out)
  }

  plan <- cm2_cm3_focus_traversal_plan(background_state, context, focus_items)
  if (is.null(plan)) {
    return(out)
  }
  focus_component_indices <- plan$focus_component_indices
  rest_component_indices <- setdiff(seq_along(components), focus_component_indices)
  rest_key <- if (length(rest_component_indices)) {
    paste(rest_component_indices, collapse = ":")
  } else {
    "."
  }
  if (!exists(rest_key, envir = background_state$rest_gamma, inherits = FALSE)) {
    assign(
      rest_key,
      cm2_cm3_convolve_component_gammas(
        background_state$component_gammas[rest_component_indices],
        context$max_total_score
      ),
      envir = background_state$rest_gamma
    )
  }
  rest_gamma <- get(rest_key, envir = background_state$rest_gamma, inherits = FALSE)
  # SKbias2.Calculate_Cprob1/2/3 add one configuration contribution at a time.
  # Retain that traversal and accumulation order; vectorized grouped sums can
  # change the last bits and would break exact bootstrap parity.
  for (grid_row in seq_along(plan$block_score)) {
    block_score <- plan$block_score[[grid_row]]
    if (block_score > total_score) {
      next
    }
    rest_value <- rest_gamma[[total_score - block_score + 1L]]
    if (rest_value <= 0) {
      next
    }
    focus_cell <- plan$focus_cell[[grid_row]]
    out[[focus_cell]] <- out[[focus_cell]] +
      plan$block_weight[[grid_row]] * rest_value / denominator
  }

  out
}

#' Internal cm2 cm3 convolve component gammas helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param gammas Internal `gammas` value used by this helper.
#' @param max_total_score Internal `max_total_score` value used by this helper.
#' @return The internal `cm2_cm3_convolve_component_gammas()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_convolve_component_gammas <- function(gammas, max_total_score) {
  out <- c(1, numeric(max_total_score))
  for (gamma in gammas) {
    out <- convolve_score_vectors(out, gamma, max_total_score)
  }
  out
}

#' Internal cm2 cm3 expected table helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param spec GLLRM model specification.
#' @param score_group_lookup Internal `score_group_lookup` value used by this helper.
#' @param probability_cache Internal `probability_cache` value used by this helper.
#' @return The internal `cm2_cm3_expected_table()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_expected_table <- function(context, state, spec, score_group_lookup, probability_cache) {
  cache <- probability_cache %||% new_cm2_cm3_focus_probability_cache(context, state)
  switch(
    spec$kind,
    item_item = cm2_cm3_expected_item_item(context, spec, NULL, cache),
    item_exogenous = cm2_cm3_expected_item_exogenous(context, spec, NULL, cache),
    item_score_group = cm2_cm3_expected_item_score_group(context, spec, score_group_lookup, cache),
    item_item_item = cm2_cm3_expected_item_item_item(context, spec, NULL, cache),
    item_item_exogenous = cm2_cm3_expected_item_item_exogenous(context, spec, NULL, cache),
    item_item_score_group = cm2_cm3_expected_item_item_score_group(context, spec, score_group_lookup, cache),
    item_exogenous_exogenous = cm2_cm3_expected_item_exogenous_exogenous(context, spec, NULL, cache),
    item_exogenous_score_group = cm2_cm3_expected_item_exogenous_score_group(context, spec, score_group_lookup, cache),
    stop("Unknown CM2/CM3 expected-table margin kind.", call. = FALSE)
  )
}

#' Internal cm2 cm3 expected item item helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @param score_groups Internal `score_groups` value used by this helper.
#' @param cache Computation cache.
#' @return The internal `cm2_cm3_expected_item_item()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_expected_item_item <- function(context, spec, score_groups, cache) {
  cm2_cm3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = c(spec$item1, spec$item2)
  )
}

#' Internal cm2 cm3 expected item exogenous helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @param score_groups Internal `score_groups` value used by this helper.
#' @param cache Computation cache.
#' @return The internal `cm2_cm3_expected_item_exogenous()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_expected_item_exogenous <- function(context, spec, score_groups, cache) {
  cm2_cm3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = spec$item,
    extra_dims = context$background_raw_max[[spec$exogenous]],
    extra_value = function(group, background_values) background_values[[spec$exogenous]]
  )
}

#' Internal cm2 cm3 expected item score group helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @param score_group_lookup Internal `score_group_lookup` value used by this helper.
#' @param cache Computation cache.
#' @return The internal `cm2_cm3_expected_item_score_group()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_expected_item_score_group <- function(context, spec, score_group_lookup, cache) {
  cm2_cm3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = spec$item,
    extra_dims = cm2_cm3_score_group_count(score_group_lookup),
    extra_value = function(group, background_values) {
      global_homogeneity_lookup_score(score_group_lookup, group$score[[1L]])
    }
  )
}

#' Internal cm2 cm3 expected item item item helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @param score_groups Internal `score_groups` value used by this helper.
#' @param cache Computation cache.
#' @return The internal `cm2_cm3_expected_item_item_item()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_expected_item_item_item <- function(context, spec, score_groups, cache) {
  cm2_cm3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = c(spec$item1, spec$item2, spec$item3)
  )
}

#' Internal cm2 cm3 expected item item exogenous helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @param score_groups Internal `score_groups` value used by this helper.
#' @param cache Computation cache.
#' @return The internal `cm2_cm3_expected_item_item_exogenous()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_expected_item_item_exogenous <- function(context, spec, score_groups, cache) {
  cm2_cm3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = c(spec$item1, spec$item2),
    extra_dims = context$background_raw_max[[spec$exogenous]],
    extra_value = function(group, background_values) background_values[[spec$exogenous]]
  )
}

#' Internal cm2 cm3 expected item item score group helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @param score_group_lookup Internal `score_group_lookup` value used by this helper.
#' @param cache Computation cache.
#' @return The internal `cm2_cm3_expected_item_item_score_group()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_expected_item_item_score_group <- function(context, spec, score_group_lookup, cache) {
  cm2_cm3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = c(spec$item1, spec$item2),
    extra_dims = cm2_cm3_score_group_count(score_group_lookup),
    extra_value = function(group, background_values) {
      global_homogeneity_lookup_score(score_group_lookup, group$score[[1L]])
    }
  )
}

#' Internal cm2 cm3 expected item exogenous exogenous helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @param score_groups Internal `score_groups` value used by this helper.
#' @param cache Computation cache.
#' @return The internal `cm2_cm3_expected_item_exogenous_exogenous()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_expected_item_exogenous_exogenous <- function(context, spec, score_groups, cache) {
  # Compatibility-sensitive EIXZ branch: source/digram_source_20260817/skunits/SKbias2.pas::
  # calculate_EIXZ_table uses Calculate_Cprob1 for the single focus item.
  cm2_cm3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = spec$item,
    extra_dims = c(context$background_raw_max[[spec$exogenous1]], context$background_raw_max[[spec$exogenous2]]),
    extra_value = function(group, background_values) {
      c(background_values[[spec$exogenous1]], background_values[[spec$exogenous2]])
    }
  )
}

#' Internal cm2 cm3 expected item exogenous score group helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @param score_group_lookup Internal `score_group_lookup` value used by this helper.
#' @param cache Computation cache.
#' @return The internal `cm2_cm3_expected_item_exogenous_score_group()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_expected_item_exogenous_score_group <- function(context, spec, score_group_lookup, cache) {
  # Compatibility-sensitive EIXZ branch: source/digram_source_20260817/skunits/SKbias2.pas::
  # calculate_EIXZ_table also covers the score-group pseudo-variable case.
  cm2_cm3_expected_focus_table(
    context = context,
    cache = cache,
    focus_items = spec$item,
    extra_dims = c(context$background_raw_max[[spec$exogenous]], cm2_cm3_score_group_count(score_group_lookup)),
    extra_value = function(group, background_values) {
      c(
        background_values[[spec$exogenous]],
        global_homogeneity_lookup_score(score_group_lookup, group$score[[1L]])
      )
    }
  )
}

#' Internal cm2 cm3 expected focus table helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param context Prepared GLLRM computation context.
#' @param cache Computation cache.
#' @param focus_items Internal `focus_items` value used by this helper.
#' @param extra_dims Internal `extra_dims` value used by this helper.
#' @param extra_value Internal `extra_value` value used by this helper.
#' @return The internal `cm2_cm3_expected_focus_table()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_expected_focus_table <- function(context,
                                       cache,
                                       focus_items,
                                       extra_dims = integer(),
                                       extra_value = NULL) {
  dims <- c(as.integer(context$item_raw_max[focus_items]), as.integer(extra_dims))
  out <- array(0, dim = dims)
  groups <- cm2_cm3_source_score_background_groups(context)
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
    probabilities <- cm2_cm3_focus_probabilities(
      cache,
      focus_items = focus_items,
      total_score = group$score[[1L]],
      background_values = background_values
    )
    out <- cm2_cm3_add_expected_probabilities(
      out = out,
      probabilities = probabilities,
      suffix = suffix,
      count = group$count[[1L]]
    )
  }
  out
}

#' Build complete-record groups for CM2/CM3 expected margins
#'
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param context Fitted GLLRM context.
#' @param rows Optional source-row policy. The default mirrors the expected
#'   CM2/CM3 routines, which require complete items and all exogenous values but
#'   do not apply the estimation score window.
#' @return A data frame of total score, record count, and exogenous values in
#'   first-occurrence source order.
#' @keywords internal
cm2_cm3_source_score_background_groups <- function(
    context,
    rows = source_complete_item_exogenous_rows(context)) {
  # Source trace: source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IX_table,
  # calculate_EIJX_table, and calculate_EIXZ_table loop over source records,
  # attach the observed exogenous/score-group value, and only use DIF
  # background values for the conditional probabilities. The fitting
  # score_exo_groups object is a smaller sufficient statistic and is not enough
  # for CM2/CM3 expected margins that mention non-DIF exogenous variables.
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
    # Source trace: calculate_expected_IJ_table/calculate_expected_IX_table in
    # source/digram_source_20260817/skunits/SKbias2.pas recompute the total from Get_Items after
    # record acceptance. Do not reuse the CML status/score sentinel here.
    row_score <- sum(as.integer(context$item_matrix[row, ]))
    values <- if (context$n_backgrounds > 0L) {
      as.integer(context$background_matrix[row, ])
    } else {
      integer()
    }
    if (length(values) && any(is.na(values) | values < 1L)) {
      next
    }
    key <- paste(c(row_score, values), collapse = "\r")
    pos <- index[[key]]
    if (is.null(pos)) {
      pos <- length(groups) + 1L
      index[[key]] <- pos
      groups[[pos]] <- c(score = row_score, count = 0L, values)
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

#' Internal cm2 cm3 add expected probabilities helper
#'
#' Supports the cm2 cm3 expected implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table`.
#' @param out Internal `out` value used by this helper.
#' @param probabilities Internal `probabilities` value used by this helper.
#' @param suffix Internal `suffix` value used by this helper.
#' @param count Internal `count` value used by this helper.
#' @return The internal `cm2_cm3_add_expected_probabilities()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_add_expected_probabilities <- function(out, probabilities, suffix, count) {
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
