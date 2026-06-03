calculate_gllrm_joint_expected_margins <- function(context, state) {
  state$expected_items <- state$item_gamma
  state$expected_items[,] <- 0
  state$expected_ld <- lapply(state$ld_parameters, function(x) x * 0)
  state$expected_dif <- lapply(state$dif_parameters, function(x) x * 0)
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  background_cache <- new.env(parent = emptyenv())

  for (group_index in seq_len(nrow(context$score_exo_groups))) {
    group <- context$score_exo_groups[group_index, , drop = FALSE]
    background_values <- gllrm_group_background_values(context, group)
    cache_key <- gllrm_background_cache_key(context, background_values)
    if (!exists(cache_key, envir = background_cache, inherits = FALSE)) {
      component_gamma <- lapply(components, function(component_items) {
        gllrm_component_gamma(context, state, component_items, background_values)
      })
      gamma_vectors <- lapply(component_gamma, `[[`, "gamma")
      assign(
        cache_key,
        list(
          component_gamma = component_gamma,
          convolutions = gllrm_component_convolutions(gamma_vectors, context$max_total_score)
        ),
        envir = background_cache
      )
    }
    cached <- get(cache_key, envir = background_cache, inherits = FALSE)
    component_gamma <- cached$component_gamma
    convolutions <- cached$convolutions
    full_gamma <- convolutions$full
    total_score <- group$score[[1L]]
    denominator <- full_gamma[[total_score + 1L]]
    if (denominator <= 0) {
      next
    }

    for (component_index in seq_along(components)) {
      rest_gamma <- convolutions$rest[[component_index]]
      state <- gllrm_accumulate_component_expected(
        context = context,
        state = state,
        component_items = components[[component_index]],
        background_values = background_values,
        total_score = total_score,
        group_count = group$count[[1L]],
        rest_gamma = rest_gamma,
        component_scale = component_gamma[[component_index]]$scale,
        component_weights = component_gamma[[component_index]]$config_weights,
        denominator = denominator
      )
    }
  }
  state
}

gllrm_component_convolutions <- function(gamma_vectors, max_total_score) {
  n <- length(gamma_vectors)
  unit <- c(1, numeric(max_total_score))
  prefix <- vector("list", n + 1L)
  suffix <- vector("list", n + 1L)
  prefix[[1L]] <- unit
  if (n > 0L) {
    for (i in seq_len(n)) {
      prefix[[i + 1L]] <- convolve_score_vectors(prefix[[i]], gamma_vectors[[i]], max_total_score)
    }
  }
  suffix[[n + 1L]] <- unit
  if (n > 0L) {
    for (i in rev(seq_len(n))) {
      suffix[[i]] <- convolve_score_vectors(gamma_vectors[[i]], suffix[[i + 1L]], max_total_score)
    }
  }
  rest <- vector("list", n)
  if (n > 0L) {
    for (i in seq_len(n)) {
      rest[[i]] <- convolve_score_vectors(prefix[[i]], suffix[[i + 1L]], max_total_score)
    }
  }
  list(full = prefix[[n + 1L]], rest = rest)
}

gllrm_group_background_values <- function(context, group) {
  if (context$n_backgrounds == 0L) {
    return(integer())
  }
  values <- as.integer(group[1L, context$backgrounds$name, drop = TRUE])
  names(values) <- context$backgrounds$name
  values
}

gllrm_background_cache_key <- function(context, background_values) {
  active <- context$active_background_indices %||% seq_along(background_values)
  if (length(active) == 0L) {
    return(".")
  }
  paste(background_values[active], collapse = "\r")
}

gllrm_accumulate_component_expected <- function(context,
                                               state,
                                               component_items,
                                               background_values,
                                               total_score,
                                               group_count,
                                               rest_gamma,
                                               component_scale,
                                               component_weights = NULL,
                                               denominator) {
  key <- gllrm_component_key(component_items)
  config_matrix <- context$component_config_matrices[[key]]
  config_scores <- context$component_config_scores[[key]]
  if (is.null(config_matrix)) {
    configs <- gllrm_component_configurations(context, component_items)
    config_matrix <- as.matrix(configs)
    storage.mode(config_matrix) <- "integer"
    config_scores <- rowSums(config_matrix)
  }
  ld_local <- context$component_ld_local_matrices[[key]] %||%
    context$component_ld_local_indices[[key]]

  if (is.null(component_weights)) {
    component_weights <- gllrm_component_config_weights_fast(
      context,
      state,
      component_items,
      config_matrix,
      background_values,
      key = key
    )
  }
  valid <- config_scores <= total_score & component_weights > 0
  if (!any(valid)) {
    return(state)
  }
  valid_index <- which(valid)
  rest_values <- rest_gamma[total_score - config_scores[valid_index] + 1L]
  expected_values <- group_count * (component_weights[valid_index] / component_scale) *
    rest_values / denominator
  positive <- expected_values > 0
  if (!any(positive)) {
    return(state)
  }
  valid_index <- valid_index[positive]
  expected_values <- expected_values[positive]

  for (local_index in seq_along(component_items)) {
    item <- component_items[[local_index]]
    item_scores <- config_matrix[valid_index, local_index] + 1L
    expected_by_score <- rowsum(expected_values, item_scores, reorder = FALSE)
    state$expected_items[item, as.integer(rownames(expected_by_score))] <-
      state$expected_items[item, as.integer(rownames(expected_by_score))] +
      as.numeric(expected_by_score[, 1L])
    dif_rows <- context$dif_by_item_matrices[[item]] %||% context$dif_by_item[[item]]
    if (nrow(dif_rows) > 0L) {
      for (dif_row in seq_len(nrow(dif_rows))) {
        background <- dif_rows[dif_row, 1L]
        dif_index <- dif_rows[dif_row, 2L]
        background_value <- background_values[[background]]
        state$expected_dif[[dif_index]][as.integer(rownames(expected_by_score)), background_value] <-
          state$expected_dif[[dif_index]][as.integer(rownames(expected_by_score)), background_value] +
          as.numeric(expected_by_score[, 1L])
      }
    }
  }
  if (!is.null(ld_local) && nrow(ld_local) > 0L) {
    for (ld_row in seq_len(nrow(ld_local))) {
      ld_index <- ld_local[ld_row, 1L]
      item1_scores <- config_matrix[valid_index, ld_local[ld_row, 2L]] + 1L
      item2_scores <- config_matrix[valid_index, ld_local[ld_row, 3L]] + 1L
      pair_index <- item1_scores + (item2_scores - 1L) * nrow(state$expected_ld[[ld_index]])
      expected_by_pair <- rowsum(expected_values, pair_index, reorder = FALSE)
      state$expected_ld[[ld_index]][as.integer(rownames(expected_by_pair))] <-
        state$expected_ld[[ld_index]][as.integer(rownames(expected_by_pair))] +
        as.numeric(expected_by_pair[, 1L])
    }
  }
  state
}

update_gllrm_active_parameters <- function(context, state, apply_update, track_delta) {
  ratio_state <- calculate_rasch_update_ratios(
    context$bundle,
    context$counts,
    state$expected_items,
    state$item_gamma,
    apply_update = apply_update
  )
  state$update_items <- ratio_state$update
  if (isTRUE(track_delta)) {
    state$delta <- max(state$delta, ratio_state$delta)
  }
  if (isTRUE(apply_update)) {
    state$item_gamma <- ratio_state$item_gamma
  }

  for (ld_index in seq_along(context$ld_specs)) {
    observed <- context$observed_ld[[ld_index]]
    expected <- state$expected_ld[[ld_index]]
    update <- expected
    update[,] <- 0
    active <- observed > 0 & expected > 0
    update[active] <- observed[active] / expected[active]
    if (isTRUE(track_delta) && length(expected) > 0L) {
      state$delta <- max(state$delta, max(abs(expected - observed)))
    }
    if (isTRUE(apply_update)) {
      state$ld_parameters[[ld_index]] <- state$ld_parameters[[ld_index]] * update
    }
    state$update_ld[[ld_index]] <- update
  }

  for (dif_index in seq_along(context$dif_specs)) {
    observed <- context$observed_dif[[dif_index]]
    expected <- state$expected_dif[[dif_index]]
    update <- expected
    update[,] <- 0
    active <- expected > 0
    update[active] <- observed[active] / expected[active]
    if (isTRUE(track_delta) && length(expected) > 0L) {
      state$delta <- max(state$delta, max(abs(expected - observed)))
    }
    if (isTRUE(apply_update)) {
      state$dif_parameters[[dif_index]] <- state$dif_parameters[[dif_index]] * update
    }
    state$update_dif[[dif_index]] <- update
  }
  state
}

update_gllrm_active_parameters_once <- function(context, state) {
  ratio_state <- calculate_rasch_update_ratios(
    context$bundle,
    context$counts,
    state$expected_items,
    state$item_gamma,
    apply_update = TRUE
  )
  state$update_items <- ratio_state$update
  state$delta <- max(state$delta, ratio_state$delta)
  state$item_gamma <- ratio_state$item_gamma

  for (ld_index in seq_along(context$ld_specs)) {
    observed <- context$observed_ld[[ld_index]]
    expected <- state$expected_ld[[ld_index]]
    update <- expected
    update[,] <- 0
    active <- observed > 0 & expected > 0
    update[active] <- observed[active] / expected[active]
    if (length(expected) > 0L) {
      state$delta <- max(state$delta, max(abs(expected - observed)))
    }
    state$ld_parameters[[ld_index]] <- state$ld_parameters[[ld_index]] * update
    state$update_ld[[ld_index]] <- update
  }

  for (dif_index in seq_along(context$dif_specs)) {
    observed <- context$observed_dif[[dif_index]]
    expected <- state$expected_dif[[dif_index]]
    update <- expected
    update[,] <- 0
    active <- expected > 0
    update[active] <- observed[active] / expected[active]
    if (length(expected) > 0L) {
      state$delta <- max(state$delta, max(abs(expected - observed)))
    }
    state$dif_parameters[[dif_index]] <- state$dif_parameters[[dif_index]] * update
    state$update_dif[[dif_index]] <- update
  }

  state$report_delta <- state$delta
  state
}

adjust_gllrm_dependency_parameters <- function(context,
                                               state,
                                               absorb_ld_item_factors = FALSE,
                                               reset_ld_reference = isTRUE(absorb_ld_item_factors),
                                               preserve_current_ties = TRUE) {
  item_log_factors <- state$item_gamma
  item_log_factors[,] <- 0
  item_zero_factors <- state$item_gamma
  item_zero_factors[,] <- FALSE
  i_ref <- 1L
  j_ref <- 1L

  for (ld_index in seq_along(context$ld_specs)) {
    spec <- context$ld_specs[[ld_index]]
    if (isTRUE(reset_ld_reference)) {
      i_ref <- 1L
      j_ref <- 1L
    }
    adjustment <- adjust_ld_gamma_source_reference_details(
      context$observed_ld[[ld_index]],
      state$ld_parameters[[ld_index]],
      i_ref = i_ref,
      j_ref = j_ref,
      preserve_current_ties = preserve_current_ties
    )
    i_ref <- adjustment$i_ref
    j_ref <- adjustment$j_ref
    state$ld_parameters[[ld_index]] <- adjustment$adjusted
    if (isTRUE(absorb_ld_item_factors) && isTRUE(adjustment$adjusted_with_item_factors)) {
      item_factors <- add_ld_item_log_factors(
        item_log_factors,
        item = spec$item1,
        factors = adjustment$i_first,
        zero_factors = item_zero_factors
      )
      item_log_factors <- item_factors$log_factors
      item_zero_factors <- item_factors$zero_factors
      item_factors <- add_ld_item_log_factors(
        item_log_factors,
        item = spec$item2,
        factors = adjustment$j_first,
        zero_factors = item_zero_factors
      )
      item_log_factors <- item_factors$log_factors
      item_zero_factors <- item_factors$zero_factors
    }
  }
  for (dif_index in seq_along(context$dif_specs)) {
    state$dif_parameters[[dif_index]] <- adjust_gllrm_dif_reference(
      context$observed_dif[[dif_index]],
      state$dif_parameters[[dif_index]]
    )
  }
  state$item_gamma <- adjust_item_gammas_source_scale_log(context$bundle, state$item_gamma, item_log_factors, item_zero_factors)
  state
}

gllrm_output_parameter_state <- function(context, state) {
  if (length(context$ld_specs) == 0L) {
    return(state)
  }
  score_reference <- context$item_score_reference %||% 0L
  large_ld_component <- max(lengths(context$ld_components_items %||% list()), 0L) >= 5L
  many_ld_terms <- length(context$ld_specs) > 3L
  if (score_reference == 0L && !large_ld_component && !many_ld_terms) {
    return(state)
  }
  preserve_output_ties <- score_reference == 0L && !large_ld_component && many_ld_terms
  # Source trace: the report-facing active GLLRM parameter gauge is obtained by
  # readjusting LD tables to carried reference rows/columns and absorbing the
  # extracted row/column factors into item gammas before the final item-gamma
  # source scaling. This is the output analogue of Adjust_IJparameters0. DIGRAM
  # applies this presentation gauge for non-zero item-score references, large
  # LD components, and multi-term LD structures such as the example EmoReg model.
  adjust_gllrm_dependency_parameters(
    context,
    state,
    absorb_ld_item_factors = TRUE,
    reset_ld_reference = score_reference != 0L,
    preserve_current_ties = preserve_output_ties
  )
}

add_ld_item_log_factors <- function(log_factors, item, factors, zero_factors) {
  for (score in seq_along(factors)) {
    if (factors[[score]] > 0) {
      log_factors[item, score] <- log_factors[item, score] + log(factors[[score]])
    } else {
      zero_factors[item, score] <- TRUE
    }
  }
  list(log_factors = log_factors, zero_factors = zero_factors)
}

adjust_item_gammas_source_scale_log <- function(bundle, item_gamma, log_factors = NULL, zero_factors = NULL) {
  items <- bundle$model$items
  log_values <- item_gamma
  log_values[,] <- -Inf
  if (is.null(log_factors)) {
    log_factors <- item_gamma
    log_factors[,] <- 0
  }
  if (is.null(zero_factors)) {
    zero_factors <- item_gamma
    zero_factors[,] <- FALSE
  }

  positive <- item_gamma > 0 & !zero_factors
  log_values[positive] <- log(item_gamma[positive]) + log_factors[positive]

  ifra <- integer(nrow(items))
  itil <- integer(nrow(items))
  valid_item <- logical(nrow(items))
  for (item_index in seq_len(nrow(items))) {
    scores <- seq.int(0L, items$raw_max[[item_index]] - 1L)
    item_logs <- as.numeric(log_values[item_index, as.character(scores)])
    positive_scores <- scores[is.finite(item_logs)]
    if (length(positive_scores) == 0L) {
      ifra[[item_index]] <- items$raw_max[[item_index]]
      itil[[item_index]] <- 0L
    } else {
      ifra[[item_index]] <- min(positive_scores)
      itil[[item_index]] <- max(positive_scores)
      valid_item[[item_index]] <- ifra[[item_index]] < itil[[item_index]]
    }
  }

  top_log_sum <- 0
  s_max <- 0L
  s_min <- 0L
  normalized_logs <- log_values
  for (item_index in seq_len(nrow(items))) {
    if (!isTRUE(valid_item[[item_index]])) {
      next
    }
    fra <- ifra[[item_index]]
    til <- itil[[item_index]]
    s_max <- s_max + til
    s_min <- s_min + fra
    alpha <- normalized_logs[item_index, as.character(fra)]
    for (score in seq.int(fra, til)) {
      normalized_logs[item_index, as.character(score)] <-
        normalized_logs[item_index, as.character(score)] - alpha
    }
    top_log_sum <- top_log_sum + normalized_logs[item_index, as.character(til)]
  }

  alpha <- 0
  if ((s_max - s_min) > 0) {
    alpha <- -top_log_sum / (s_max - s_min)
  }

  out <- item_gamma
  out[,] <- 0
  for (item_index in seq_len(nrow(items))) {
    if (!isTRUE(valid_item[[item_index]])) {
      next
    }
    fra <- ifra[[item_index]]
    til <- itil[[item_index]]
    for (score in seq.int(fra, til)) {
      out[item_index, as.character(score)] <-
        exp(normalized_logs[item_index, as.character(score)] + (score - fra) * alpha)
    }
  }
  for (item_index in seq_len(nrow(items))) {
    if (!isTRUE(valid_item[[item_index]]) && ifra[[item_index]] <= itil[[item_index]]) {
      out[item_index, as.character(itil[[item_index]])] <- 1
    }
  }
  out
}

adjust_gllrm_dif_reference <- function(observed, gamma) {
  rows <- nrow(gamma)
  cols <- ncol(gamma)
  i_cells <- rowSums(observed > 0)
  j_cells <- colSums(observed > 0)
  i_ref <- 1L
  if (i_cells[[i_ref]] < cols) {
    i_ref <- which.max(i_cells)
  }
  j_ref <- 1L
  if (j_cells[[j_ref]] < rows) {
    j_ref <- which.max(j_cells)
  }
  ix_ref <- gamma[i_ref, j_ref]
  i_pos <- sum(i_cells > 0)
  j_pos <- sum(j_cells > 0)
  i_first <- gamma[, j_ref]
  for (i in seq_len(rows)) {
    if (i_first[[i]] == 0) {
      hit <- which(observed[i, ] > 0)
      if (length(hit) > 0L) {
        i_first[[i]] <- gamma[i, hit[[1L]]]
      }
    }
  }
  j_first <- gamma[i_ref, ]
  for (j in seq_len(cols)) {
    if (j_first[[j]] == 0) {
      hit <- which(observed[, j] > 0)
      if (length(hit) > 0L) {
        j_first[[j]] <- gamma[hit[[1L]], j]
      }
    }
  }
  adjusted <- gamma
  if (i_pos > 1L && j_pos > 1L && ix_ref > 0) {
    old <- gamma
    for (i in seq_len(rows)) {
      for (j in seq_len(cols)) {
        if (observed[i, j] > 0) {
          denominator <- j_first[[j]] * i_first[[i]]
          adjusted[i, j] <- if (denominator > 0) old[i, j] * ix_ref / denominator else 0
        } else {
          adjusted[i, j] <- 0
        }
      }
    }
  } else {
    adjusted[,] <- 0
    adjusted[observed > 0] <- 1
  }
  adjusted
}

active_gllrm_loglike <- function(context, state) {
  loglike <- 0
  score_gamma_cache <- new.env(parent = emptyenv())
  for (row in context$valid_rows) {
    score <- context$score[[row]]
    background_values <- if (context$n_backgrounds > 0L) context$background_matrix[row, ] else integer()
    cache_key <- gllrm_background_cache_key(context, background_values)
    if (!exists(cache_key, envir = score_gamma_cache, inherits = FALSE)) {
      assign(cache_key, active_gllrm_score_gamma(context, state, background_values), envir = score_gamma_cache)
    }
    score_gamma <- get(cache_key, envir = score_gamma_cache, inherits = FALSE)
    product_gamma <- 1
    for (item in seq_len(context$n_items)) {
      item_score <- context$item_matrix[row, item]
      product_gamma <- product_gamma * state$item_gamma[item, item_score + 1L]
      dif_rows <- context$dif_by_item_matrices[[item]] %||% context$dif_by_item[[item]]
      if (nrow(dif_rows) > 0L) {
        for (dif_row in seq_len(nrow(dif_rows))) {
          background <- dif_rows[dif_row, 1L]
          dif_index <- dif_rows[dif_row, 2L]
          product_gamma <- product_gamma *
            state$dif_parameters[[dif_index]][item_score + 1L, background_values[[background]]]
      }
    }
    }
    for (ld_index in seq_along(context$ld_specs)) {
      spec <- context$ld_specs[[ld_index]]
      product_gamma <- product_gamma *
        state$ld_parameters[[ld_index]][
          context$item_matrix[row, spec$item1] + 1L,
          context$item_matrix[row, spec$item2] + 1L
        ]
    }
    probability <- product_gamma / score_gamma[[score + 1L]]
    if (probability > 0) {
      loglike <- loglike - log(probability)
    }
  }
  loglike
}

active_gllrm_score_gamma <- function(context, state, background_values) {
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  full <- c(1, numeric(context$max_total_score))
  for (component_items in components) {
    component <- gllrm_component_gamma(context, state, component_items, background_values)
    full <- convolve_score_vectors(full, component$gamma * component$scale, context$max_total_score)
  }
  full
}

fit_gllrm_active <- function(spec,
                             max_step = 5000L,
                             max_delta = 0.0001,
                             max_joint_configs = 200000L,
                             bundle = NULL) {
  bundle <- bundle %||% build_item_parameters_bundle(spec$project)
  context <- build_gllrm_active_context(spec, bundle, max_joint_configs = max_joint_configs)
  state <- initialize_gllrm_active_state(context)
  if (context$counts$n_valid == 0L || context$n_items == 0L) {
    state$converged <- TRUE
    return(list(context = context, state = state, bundle = bundle))
  }

  previous_delta <- Inf
  delta_history <- numeric(max_step)
  finish_count <- 0L
  while (state$n_step < max_step) {
    state$n_step <- state$n_step + 1L
    previous_delta <- state$delta
    state <- calculate_gllrm_joint_expected_margins(context, state)
    state$delta <- 0
    state <- update_gllrm_active_parameters_once(context, state)
    # Source trace: the runtime report path reselects the densest observed
    # row/column reference for active LD tables during IPF. This is distinct
    # from the final report-facing gauge adjustment in
    # gllrm_output_parameter_state().
    state <- adjust_gllrm_dependency_parameters(
      context,
      state,
      absorb_ld_item_factors = FALSE,
      preserve_current_ties = FALSE
    )
    if (previous_delta <= state$delta) {
      finish_count <- finish_count + 1L
    } else {
      finish_count <- 0L
    }
    delta_history[[state$n_step]] <- state$delta
    if (state$delta < max_delta) {
      state$converged <- TRUE
      break
    }
    repeated_delta <- state$delta == previous_delta
    repeated_two_back <- state$n_step > 5L && state$delta == delta_history[[state$n_step - 2L]]
    source_periodic_stop <- state$n_step %% 1000L == 0L ||
      (state$n_step %% 50L == 0L && state$delta > 10)
    if (source_periodic_stop || repeated_delta || repeated_two_back || finish_count > 10L) {
      break
    }
  }
  if (!state$converged || state$report_delta == 0) {
    state$report_delta <- state$delta
  }
  state <- calculate_gllrm_joint_expected_margins(context, state)
  state$delta <- 0
  state <- update_gllrm_active_parameters(context, state, apply_update = FALSE, track_delta = TRUE)
  state$log_likelihood <- active_gllrm_loglike(context, state)
  list(context = context, state = state, bundle = bundle)
}

gllrm_spec_ld_term_strings <- function(spec) {
  if (is.null(spec$ld) || nrow(spec$ld) == 0L) {
    return(character())
  }
  paste(spec$ld$item1, spec$ld$item2, sep = ":")
}

gllrm_spec_dif_term_strings <- function(spec) {
  if (is.null(spec$dif) || nrow(spec$dif) == 0L) {
    return(character())
  }
  paste(spec$dif$item, spec$dif$exogenous, sep = ":")
}

fit_gllrm_with_added_ld <- function(object, item1, item2, max_step, max_delta) {
  context <- object$fit$context
  ld_terms <- source_order_ld_terms(
    context$items,
    c(
      gllrm_spec_ld_term_strings(object$spec),
      paste(context$items$name[[item1]], context$items$name[[item2]], sep = ":")
    )
  )
  spec <- gllrm(
    object$analysis,
    ld = gRm_terms_formula(ld_terms),
    dif = gRm_terms_formula(gllrm_spec_dif_term_strings(object$spec))
  )
  fit(spec, max_step = max_step, max_delta = max_delta)
}

fit_gllrm_with_added_dif <- function(object, item, background, max_step, max_delta) {
  context <- object$fit$context
  dif_terms <- source_order_dif_terms(
    context$items,
    context$backgrounds,
    c(
      gllrm_spec_dif_term_strings(object$spec),
      paste(context$items$name[[item]], context$backgrounds$name[[background]], sep = ":")
    )
  )
  spec <- gllrm(
    object$analysis,
    ld = gRm_terms_formula(gllrm_spec_ld_term_strings(object$spec)),
    dif = gRm_terms_formula(dif_terms)
  )
  fit(spec, max_step = max_step, max_delta = max_delta)
}

source_order_ld_terms <- function(items, terms) {
  if (length(terms) == 0L) {
    return(character())
  }
  parts <- strsplit(terms, ":", fixed = TRUE)
  item1 <- vapply(parts, `[[`, character(1L), 1L)
  item2 <- vapply(parts, `[[`, character(1L), 2L)
  index1 <- match(item1, items$name)
  index2 <- match(item2, items$name)
  lo <- pmin(index1, index2)
  hi <- pmax(index1, index2)
  unique(paste(items$name[lo], items$name[hi], sep = ":")[order(lo, hi)])
}

source_order_dif_terms <- function(items, backgrounds, terms) {
  if (length(terms) == 0L) {
    return(character())
  }
  parts <- strsplit(terms, ":", fixed = TRUE)
  item <- vapply(parts, `[[`, character(1L), 1L)
  background <- vapply(parts, `[[`, character(1L), 2L)
  item_index <- match(item, items$name)
  background_index <- match(background, backgrounds$name)
  unique(paste(items$name[item_index], backgrounds$name[background_index], sep = ":")[order(item_index, background_index)])
}
