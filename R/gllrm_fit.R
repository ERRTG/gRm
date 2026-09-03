#' Source trace: source/digram_source_20260817/skunits/skbias22.pas::CalculateBiasedGammaValues2 and
#' source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM compute expected item, IJ,
#' and IX margins from the current GLLRM parameters. The R code performs
#' the same score-conditioned summation using explicit context objects instead
#' of Pascal global arrays.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @return The internal `calculate_gllrm_joint_expected_margins_r()` computation result.
#' @keywords internal
#' @noRd
calculate_gllrm_joint_expected_margins_r <- function(context, state) {
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

#' The R implementation above is the readable source-shaped reference for
#' source/digram_source_20260817/skunits/skbias22.pas::CalculateBiasedGammaValues2. The dispatcher uses a
#' native backend for speed when available; the backend must remain a mechanical
#' translation of this R reference and is tested against it directly.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @return The internal `calculate_gllrm_joint_expected_margins()` computation result.
#' @keywords internal
#' @noRd
calculate_gllrm_joint_expected_margins <- function(context, state) {
  native <- calculate_gllrm_joint_expected_margins_cpp(context, state)
  if (!is.null(native)) {
    return(native)
  }
  calculate_gllrm_joint_expected_margins_r(context, state)
}

#' Internal gllrm expected native available helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @return The internal `gllrm_expected_native_available()` computation result.
#' @keywords internal
#' @noRd
gllrm_expected_native_available <- function() {
  is.loaded("gRm_gllrm_expected_margins", PACKAGE = "gRm")
}

#' Internal calculate gllrm joint expected margins cpp helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @return The internal `calculate_gllrm_joint_expected_margins_cpp()` computation result.
#' @keywords internal
#' @noRd
calculate_gllrm_joint_expected_margins_cpp <- function(context, state) {
  if (!gllrm_expected_native_available()) {
    return(NULL)
  }
  native_input <- gllrm_expected_native_input(context)
  out <- .Call(
    "gRm_gllrm_expected_margins",
    native_input,
    state$item_gamma,
    state$ld_parameters,
    state$dif_parameters,
    PACKAGE = "gRm"
  )
  state$expected_items <- out$expected_items
  state$expected_ld <- out$expected_ld
  state$expected_dif <- out$expected_dif
  state
}

#' Internal gllrm expected native input helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `gllrm_expected_native_input()` computation result.
#' @keywords internal
#' @noRd
gllrm_expected_native_input <- function(context) {
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  keys <- vapply(components, gllrm_component_key, character(1L))
  group_backgrounds <- if (context$n_backgrounds > 0L) {
    out <- as.matrix(context$score_exo_groups[, context$backgrounds$name, drop = FALSE])
    storage.mode(out) <- "integer"
    out
  } else {
    matrix(integer(), nrow = nrow(context$score_exo_groups), ncol = 0L)
  }

  list(
    components = unname(components),
    config_matrices = unname(context$component_config_matrices[keys]),
    config_scores = lapply(unname(context$component_config_scores[keys]), as.integer),
    ld_local_matrices = unname(context$component_ld_local_matrices[keys]),
    dif_by_item = lapply(context$dif_by_item_matrices, function(x) {
      storage.mode(x) <- "integer"
      x
    }),
    item_raw_max = as.integer(context$item_raw_max),
    max_total_score = as.integer(context$max_total_score),
    group_scores = as.integer(context$score_exo_groups$score),
    group_counts = as.numeric(context$score_exo_groups$count),
    group_backgrounds = group_backgrounds,
    dif_backgrounds = as.integer(context$dif_background_indices %||% integer())
  )
}

#' Source trace: source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation combines
#' component score polynomials into whole-model score probabilities. This helper
#' keeps the same convolution calculation in R list form.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param gamma_vectors Internal `gamma_vectors` value used by this helper.
#' @param max_total_score Internal `max_total_score` value used by this helper.
#' @return The internal `gllrm_component_convolutions()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal gllrm group background values helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param context Prepared GLLRM computation context.
#' @param group Internal `group` value used by this helper.
#' @return The internal `gllrm_group_background_values()` computation result.
#' @keywords internal
#' @noRd
gllrm_group_background_values <- function(context, group) {
  if (context$n_backgrounds == 0L) {
    return(integer())
  }
  values <- as.integer(group[1L, context$backgrounds$name, drop = TRUE])
  names(values) <- context$backgrounds$name
  values
}

#' Internal gllrm background cache key helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param context Prepared GLLRM computation context.
#' @param background_values Internal `background_values` value used by this helper.
#' @return The internal `gllrm_background_cache_key()` computation result.
#' @keywords internal
#' @noRd
gllrm_background_cache_key <- function(context, background_values) {
  dif_backgrounds <- context$dif_background_indices %||% seq_along(background_values)
  if (length(dif_backgrounds) == 0L) {
    return(".")
  }
  paste(background_values[dif_backgrounds], collapse = "\r")
}

#' Source trace: source/digram_source_20260817/skunits/skbias22.pas::CalculateBiasedGammaValues2 accumulates
#' fitted item, IJ, and IX cells while traversing component configurations. This
#' R helper contributes one component's fitted margins to the whole-model
#' expected margins.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param component_items Internal `component_items` value used by this helper.
#' @param background_values Internal `background_values` value used by this helper.
#' @param total_score Internal `total_score` value used by this helper.
#' @param group_count Internal `group_count` value used by this helper.
#' @param rest_gamma Internal `rest_gamma` value used by this helper.
#' @param component_scale Internal `component_scale` value used by this helper.
#' @param component_weights Internal `component_weights` value used by this helper.
#' @param denominator Internal `denominator` value used by this helper.
#' @return The internal `gllrm_accumulate_component_expected()` computation result.
#' @keywords internal
#' @noRd
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

#' Source trace: source/digram_source_20260817/skunits/skbias22.pas::Find_new_IJparameters and
#' source/digram_source_20260817/skunits/skbias22.pas::Find_new_IXparameters update included LD and DIF
#' parameters by observed/fitted margin ratios, while the item score parameters
#' receive the analogous source IPF update. This core function performs those
#' multiplicative updates for one iteration.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param apply_update Internal `apply_update` value used by this helper.
#' @param track_delta Internal `track_delta` value used by this helper.
#' @return The internal `update_gllrm_parameters_core()` computation result.
#' @keywords internal
#' @noRd
update_gllrm_parameters_core <- function(context, state, apply_update, track_delta) {
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
    estimable <- observed > 0 & expected > 0
    update[estimable] <- observed[estimable] / expected[estimable]
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
    estimable <- expected > 0
    update[estimable] <- observed[estimable] / expected[estimable]
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

#' Internal update gllrm parameters helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param apply_update Internal `apply_update` value used by this helper.
#' @param track_delta Internal `track_delta` value used by this helper.
#' @return The internal `update_gllrm_parameters()` computation result.
#' @keywords internal
#' @noRd
update_gllrm_parameters <- function(context, state, apply_update, track_delta) {
  update_gllrm_parameters_core(
    context = context,
    state = state,
    apply_update = apply_update,
    track_delta = track_delta
  )
}

#' Internal update gllrm parameters once helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @return The internal `update_gllrm_parameters_once()` computation result.
#' @keywords internal
#' @noRd
update_gllrm_parameters_once <- function(context, state) {
  state <- update_gllrm_parameters_core(
    context = context,
    state = state,
    apply_update = TRUE,
    track_delta = TRUE
  )
  state$report_delta <- state$delta
  state
}

#' Source trace: source/digram_source_20260817/skunits/skbias22.pas::Adjust_IJparameters and
#' source/digram_source_20260817/skunits/skbias22.pas::Adjust_IJparameters0 move local-dependence effects
#' into the source reporting gauge after the multiplicative updates. This R
#' helper preserves that gauge convention before parameters are exposed.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param absorb_ld_item_factors Internal `absorb_ld_item_factors` value used by this helper.
#' @param reset_ld_reference Internal `reset_ld_reference` value used by this helper.
#' @param preserve_current_ties Internal `preserve_current_ties` value used by this helper.
#' @param initial_i_ref Internal `initial_i_ref` value used by this helper.
#' @param initial_j_ref Internal `initial_j_ref` value used by this helper.
#' @return The internal `adjust_gllrm_dependency_parameters()` computation result.
#' @keywords internal
#' @noRd
adjust_gllrm_dependency_parameters <- function(context,
                                               state,
                                               absorb_ld_item_factors = FALSE,
                                               reset_ld_reference = isTRUE(absorb_ld_item_factors),
                                               preserve_current_ties = TRUE,
                                               initial_i_ref = NULL,
                                               initial_j_ref = NULL) {
  item_log_factors <- state$item_gamma
  item_log_factors[,] <- 0
  item_zero_factors <- state$item_gamma
  item_zero_factors[,] <- FALSE
  source_ref <- as.integer(context$item_score_reference %||% 0L) + 1L
  i_ref <- as.integer(initial_i_ref %||% source_ref)
  j_ref <- as.integer(initial_j_ref %||% source_ref)

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

#' Internal gllrm output parameter state helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @return The internal `gllrm_output_parameter_state()` computation result.
#' @keywords internal
#' @noRd
gllrm_output_parameter_state <- function(context, state) {
  if (length(context$ld_specs) == 0L) {
    return(state)
  }
  score_reference <- context$item_score_reference %||% 0L
  if (score_reference == 0L) {
    return(state)
  }
  # Source trace: source/digram_source_20260817/skunits/skbias22.pas::GLLRM_estim applies the
  # final Adjust_IJparameters0/Adjust_IXparameters/Adjust_gamma0 block only when
  # (IscoreRef <> 0) and (Nij > 0). Component size and term count are not source
  # conditions for this report-facing gauge.
  adjust_gllrm_dependency_parameters(
    context,
    state,
    absorb_ld_item_factors = TRUE,
    reset_ld_reference = TRUE,
    preserve_current_ties = FALSE
  )
}

#' Internal add ld item log factors helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param log_factors Internal `log_factors` value used by this helper.
#' @param item One-based item index.
#' @param factors Internal `factors` value used by this helper.
#' @param zero_factors Internal `zero_factors` value used by this helper.
#' @return The internal `add_ld_item_log_factors()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal adjust item gammas source scale log helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param bundle Source-shaped analysis bundle.
#' @param item_gamma Internal `item_gamma` value used by this helper.
#' @param log_factors Internal `log_factors` value used by this helper.
#' @param zero_factors Internal `zero_factors` value used by this helper.
#' @return The internal `adjust_item_gammas_source_scale_log()` computation result.
#' @keywords internal
#' @noRd
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

#' Source trace: source/digram_source_20260817/skunits/skbias22.pas::Adjust_IXparameters adjusts included DIF
#' parameters against the selected reference background category. The R helper
#' applies the same reporting gauge to the fitted IX matrix.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param observed Internal `observed` value used by this helper.
#' @param gamma Internal `gamma` value used by this helper.
#' @return The internal `adjust_gllrm_dif_reference()` computation result.
#' @keywords internal
#' @noRd
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

#' Source trace: source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM and
#' source/digram_source_20260817/skunits/skbias22.pas::GLLRM_estim evaluate the fitted current GLLRM
#' likelihood from the same score-conditioned probabilities used for expected
#' margins. The R helper returns the DIGRAM/source-style negative conditional
#' log likelihood. The public logLik() method converts this stored value to R's
#' usual log-likelihood sign convention.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @return The internal `gllrm_loglike()` computation result.
#' @keywords internal
#' @noRd
gllrm_loglike <- function(context, state) {
  negative_loglike <- 0
  score_gamma_cache <- new.env(parent = emptyenv())
  for (row in context$valid_rows) {
    score <- context$score[[row]]
    background_values <- if (context$n_backgrounds > 0L) context$background_matrix[row, ] else integer()
    cache_key <- gllrm_background_cache_key(context, background_values)
    if (!exists(cache_key, envir = score_gamma_cache, inherits = FALSE)) {
      assign(cache_key, gllrm_score_gamma(context, state, background_values), envir = score_gamma_cache)
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
      negative_loglike <- negative_loglike - log(probability)
    }
  }
  negative_loglike
}

#' Internal gllrm score gamma helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param background_values Internal `background_values` value used by this helper.
#' @return The internal `gllrm_score_gamma()` computation result.
#' @keywords internal
#' @noRd
gllrm_score_gamma <- function(context, state, background_values) {
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  full <- c(1, numeric(context$max_total_score))
  for (component_items in components) {
    component <- gllrm_component_gamma(context, state, component_items, background_values)
    full <- convolve_score_vectors(full, component$gamma * component$scale, context$max_total_score)
  }
  full
}

#' Source trace: source/digram_source_20260817/skunits/skbias22.pas::GLLRM_estim and
#' source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM run the iterative GLLRM
#' estimation loop: initialize parameters, compute expected margins,
#' update by observed/fitted ratios, adjust the reporting gauge, and stop at the
#' source tolerance or iteration limit. This R function is the package's direct
#' implementation of that loop.
#' @param spec GLLRM model specification.
#' @param max_step Maximum fitting iteration.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @param max_joint_configs Internal `max_joint_configs` value used by this helper.
#' @param bundle Source-shaped analysis bundle.
#' @param initial_parameters Optional saved observed-data parameter state used
#'   only by source-shaped bootstrap refits.
#' @return The internal `fit_gllrm()` computation result.
#' @keywords internal
#' @noRd
fit_gllrm <- function(spec,
                             max_step = 5000L,
                             max_delta = 0.0001,
                             max_joint_configs = 200000L,
                             bundle = NULL,
                             initial_parameters = NULL) {
  bundle <- bundle %||% build_item_parameters_bundle(spec$project)
  context <- build_gllrm_context(spec, bundle, max_joint_configs = max_joint_configs)
  # Ordinary fits take InitializeParameters' unit start. CM3 bootstrap refits
  # instead take skbias12b.UseItemParameters after SKbias8.MakeIRTcopy saved the
  # observed fit immediately before simulation.
  state <- if (is.null(initial_parameters)) {
    initialize_gllrm_state(context)
  } else {
    initialize_gllrm_state_from_parameters(context, initial_parameters)
  }
  if (context$counts$n_valid == 0L || context$n_items == 0L) {
    state$converged <- TRUE
    return(list(context = context, state = state, bundle = bundle))
  }

  control <- source_gllrm_control_state(
    n_valid = context$counts$n_valid,
    max_step = max_step,
    max_delta = max_delta
  )
  state$delta <- control$delta
  repeat {
    state$n_step <- control$n_step + 1L
    state <- calculate_gllrm_joint_expected_margins(context, state)
    state$delta <- 0
    state <- update_gllrm_parameters_once(context, state)
    # Source trace: the runtime report path reselects the densest observed
    # row/column reference for included LD tables during IPF. This is distinct
    # from the final report-facing gauge adjustment in
    # gllrm_output_parameter_state().
    state <- adjust_gllrm_dependency_parameters(
      context,
      state,
      absorb_ld_item_factors = FALSE,
      preserve_current_ties = TRUE,
      initial_i_ref = as.integer(context$item_score_reference %||% 0L) + 1L,
      initial_j_ref = as.integer(context$item_score_reference %||% 0L) + 1L
    )
    observed <- source_gllrm_observe_delta(control, state$delta)
    control <- observed$control
    state$n_step <- control$n_step
    if (isTRUE(observed$decision$stop)) {
      break
    }
  }

  # The source report uses the discrepancy that triggered the stop. A final
  # expected-margin pass below is diagnostic and must not overwrite it.
  state$report_delta <- control$delta
  state$stop_reason <- control$stop_reason
  state$previous_delta <- control$previous_delta
  state$initial_delta <- control$initial_delta
  state$min_delta <- control$min_delta
  state$min_delta_step <- control$min_delta_step
  state$finish_count <- control$n_finish
  state$recurring_delta_values <- control$recurring
  state$delta_history <- unname(control$delta_history[as.character(seq_len(control$n_step))])
  state$convergence_before_final_acceptance <- control$convergence
  state$converged <- source_gllrm_final_convergence(control)

  state <- calculate_gllrm_joint_expected_margins(context, state)
  state$delta <- 0
  state <- update_gllrm_parameters(context, state, apply_update = FALSE, track_delta = TRUE)
  state$log_likelihood <- gllrm_loglike(context, state)
  list(context = context, state = state, bundle = bundle)
}

#' Internal fit gllrm with added ld full refit helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param object Object dispatched to this helper.
#' @param item1 Internal `item1` value used by this helper.
#' @param item2 Internal `item2` value used by this helper.
#' @param max_step Maximum fitting iteration.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @return The internal `fit_gllrm_with_added_ld_full_refit()` computation result.
#' @keywords internal
#' @noRd
fit_gllrm_with_added_ld_full_refit <- function(object, item1, item2, max_step, max_delta) {
  model <- object$model %||% object$spec
  context <- object$fit$context
  candidate <- canonical_ld_term(
    vars = c(context$items$name[[item1]], context$items$name[[item2]]),
    items = context$items$name,
    label = paste(context$items$name[[item1]], context$items$name[[item2]], sep = ":"),
    source = "user"
  )
  ld_terms <- source_order_ld_table(
    context$items,
    rbind_fill(model$ld %||% empty_ld_terms(), candidate)
  )
  dif_terms <- source_order_dif_table(
    context$items,
    context$backgrounds,
    model$dif %||% empty_dif_terms()
  )
  spec <- new_gRm_model_from_canonical_terms(
    model,
    ld = ld_terms,
    dif = dif_terms,
    call = match.call()
  )
  fit(spec, max_step = max_step, max_delta = max_delta)
}

#' Internal fit gllrm with added dif full refit helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param object Object dispatched to this helper.
#' @param item One-based item index.
#' @param background One-based exogenous-variable index.
#' @param max_step Maximum fitting iteration.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @return The internal `fit_gllrm_with_added_dif_full_refit()` computation result.
#' @keywords internal
#' @noRd
fit_gllrm_with_added_dif_full_refit <- function(object, item, background, max_step, max_delta) {
  model <- object$model %||% object$spec
  context <- object$fit$context
  candidate <- canonical_dif_term(
    vars = c(context$items$name[[item]], context$backgrounds$name[[background]]),
    items = context$items$name,
    exogenous = context$backgrounds$name,
    label = paste(context$items$name[[item]], context$backgrounds$name[[background]], sep = ":"),
    source = "user"
  )
  ld_terms <- source_order_ld_table(
    context$items,
    model$ld %||% empty_ld_terms()
  )
  dif_terms <- source_order_dif_table(
    context$items,
    context$backgrounds,
    rbind_fill(model$dif %||% empty_dif_terms(), candidate)
  )
  spec <- new_gRm_model_from_canonical_terms(
    model,
    ld = ld_terms,
    dif = dif_terms,
    call = match.call()
  )
  fit(spec, max_step = max_step, max_delta = max_delta)
}

#' Internal fit gllrm with added ld helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param object Object dispatched to this helper.
#' @param item1 Internal `item1` value used by this helper.
#' @param item2 Internal `item2` value used by this helper.
#' @param max_step Maximum fitting iteration.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @return The internal `fit_gllrm_with_added_ld()` computation result.
#' @keywords internal
#' @noRd
fit_gllrm_with_added_ld <- function(object, item1, item2, max_step, max_delta) {
  fit_gllrm_candidate_ld(
    object,
    item1 = item1,
    item2 = item2,
    max_step = max_step,
    max_delta = max_delta
  )
}

#' Internal fit gllrm with added dif helper
#'
#' Supports the gllrm fit implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param object Object dispatched to this helper.
#' @param item One-based item index.
#' @param background One-based exogenous-variable index.
#' @param max_step Maximum fitting iteration.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @return The internal `fit_gllrm_with_added_dif()` computation result.
#' @keywords internal
#' @noRd
fit_gllrm_with_added_dif <- function(object, item, background, max_step, max_delta) {
  fit_gllrm_candidate_dif(
    object,
    item = item,
    background = background,
    max_step = max_step,
    max_delta = max_delta
  )
}
