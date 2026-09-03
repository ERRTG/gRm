#' Source trace: source/digram_source_20260817/skunits/skbias22.pas::LD_Gamma_calculation and
#' source/digram_source_20260817/skunits/skbias22.pas::CalculateBiasedGammaValues2 evaluate score
#' polynomials component-wise when local-dependence terms join items. The R code
#' first identifies those connected item components from the included LD graph.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `gllrm_ld_components()` computation result.
#' @keywords internal
#' @noRd
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

#' Source trace: source/digram_source_20260817/skunits/skbias22.pas::InitializeParameters and
#' source/digram_source_20260817/skunits/skbias12b.pas::InitializeParameters initialize item score
#' gammas and the included IJ/IX parameters before the iterative updates. R stores
#' the same starting values in matrices/lists rather than Pascal arrays.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `initialize_gllrm_state()` computation result.
#' @keywords internal
#' @noRd
initialize_gllrm_state <- function(context) {
  item_gamma <- initial_gllrm_item_gamma(context)
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
    # Source trace: source/digram_source_20260817/skunits/skbias22.pas::GLLRM_estim initializes
    # Delta to Nvalid before the first Take_an_IPF_step copies it to
    # previous_delta.
    delta = as.numeric(context$counts$n_valid),
    report_delta = 0,
    n_step = 0L,
    converged = TRUE,
    stop_reason = NA_character_,
    delta_history = numeric(),
    previous_delta = NA_real_,
    initial_delta = NA_real_,
    min_delta = NA_real_,
    min_delta_step = NA_integer_,
    finish_count = 0L,
    recurring_delta_values = FALSE,
    convergence_before_final_acceptance = TRUE,
    attempted_n_step = NA_integer_,
    attempted_delta = NA_real_,
    attempted_converged = NA,
    attempted_stop_reason = NA_character_,
    reported_checkpoint_step = NA_integer_,
    report_value_source = NA_character_,
    log_likelihood = NA_real_
  )
}

#' Initialize a bootstrap GLLRM refit from the saved observed parameters
#'
#' `SKbias8.Start_random_Gllrm_with_exogene` calls `MakeIRTcopy` before the
#' simulation loop. Each later `Estimate_the_GLLRM` call passes
#' `ParametersInitialized = false`, but `skbias12b.Estimate_GLLRM` then takes
#' its `UseItemParameters` branch because that saved IRT copy exists. Thus
#' every bootstrap refit starts from the same observed-data item, IJ, and IX
#' parameters; it does not start from ones and it does not warm-start from the
#' preceding replicate.
#' Source: `source/digram_source_20260817/skunits/skbias12b.pas::UseItemParameters`.
#'
#' @param context Prepared context for one generated bootstrap sample.
#' @param parameters Saved observed-data GLLRM parameter state corresponding to
#'   Pascal's `ItemParmCopy`, `IJcopy`, and `IXcopy` structures.
#' @return A fresh iterative state with copied parameters and reset work and
#'   convergence fields.
#' @keywords internal
#' @noRd
initialize_gllrm_state_from_parameters <- function(context, parameters) {
  state <- initialize_gllrm_state(context)
  required <- c("item_gamma", "ld_parameters", "dif_parameters")
  if (!is.list(parameters) || !all(required %in% names(parameters))) {
    stop("Saved GLLRM bootstrap parameters are incomplete.", call. = FALSE)
  }
  if (!identical(dim(parameters$item_gamma), dim(state$item_gamma)) ||
      length(parameters$ld_parameters) != length(state$ld_parameters) ||
      length(parameters$dif_parameters) != length(state$dif_parameters)) {
    stop("Saved GLLRM bootstrap parameters do not match the sample model.",
         call. = FALSE)
  }

  # skbias12b.UseItemParameters first copies ItemParmCopy, then recomputes
  # Ifra/Itil from the generated sample's item margins. Gamma_calculation1
  # traverses only that inclusive score range. Zeroing values outside the same
  # range gives the full-width R matrices the identical effective support;
  # an unobserved interior score remains copied because Pascal still traverses
  # it whenever it lies between Ifra and Itil.
  item_gamma <- parameters$item_gamma
  observed <- context$counts$item_counts
  for (item in seq_len(context$n_items)) {
    scores <- seq.int(0L, context$item_raw_max[[item]] - 1L)
    present <- scores[as.numeric(observed[item, seq_along(scores)]) > 0]
    if (!length(present)) {
      item_gamma[item, ] <- 0
    } else {
      outside <- scores < min(present) | scores > max(present)
      item_gamma[item, which(outside)] <- 0
    }
  }
  state$item_gamma <- item_gamma

  copy_parameter_list <- function(saved, fresh, label) {
    for (index in seq_along(fresh)) {
      if (!identical(dim(saved[[index]]), dim(fresh[[index]]))) {
        stop("Saved GLLRM ", label, " parameters do not match the sample model.",
             call. = FALSE)
      }
    }
    saved
  }
  # skbias12b.UseItemParameters copies every matching IJ/IX cell from the
  # pre-bootstrap structures after allocating fresh per-sample information.
  state$ld_parameters <- copy_parameter_list(
    parameters$ld_parameters, state$ld_parameters, "IJ"
  )
  state$dif_parameters <- copy_parameter_list(
    parameters$dif_parameters, state$dif_parameters, "IX"
  )
  state$expected_items[,] <- 0
  state$expected_ld <- lapply(state$ld_parameters, function(value) value * 0)
  state$expected_dif <- lapply(state$dif_parameters, function(value) value * 0)
  state
}

#' Internal initial gllrm item gamma helper
#'
#' Supports the gllrm components implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `initial_gllrm_item_gamma()` computation result.
#' @keywords internal
#' @noRd
initial_gllrm_item_gamma <- function(context) {
  # Source trace: source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM initializes item
  # parameters only for item scores observed in the GLLRM estimation range.
  gamma <- initial_item_gamma(context$bundle)
  gamma[,] <- 0
  observed <- context$counts$item_counts
  for (item_index in seq_len(context$n_items)) {
    cols <- seq_len(context$item_raw_max[[item_index]])
    gamma[item_index, cols] <- as.numeric(observed[item_index, cols] > 0)
  }
  gamma
}

#' Internal gllrm dif map helper
#'
#' Supports the gllrm components implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `gllrm_dif_map()` computation result.
#' @keywords internal
#' @noRd
gllrm_dif_map <- function(context) {
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

#' Internal gllrm dif by item helper
#'
#' Supports the gllrm components implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `gllrm_dif_by_item()` computation result.
#' @keywords internal
#' @noRd
gllrm_dif_by_item <- function(context) {
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

#' Internal gllrm component key helper
#'
#' Supports the gllrm components implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param component_items Internal `component_items` value used by this helper.
#' @return The internal `gllrm_component_key()` computation result.
#' @keywords internal
#' @noRd
gllrm_component_key <- function(component_items) {
  paste(as.integer(component_items), collapse = ":")
}

#' Source trace: source/digram_source_20260817/skunits/skbias22.pas::LD_Gamma_calculation loops
#' over response configurations inside each LD-connected item component. This R
#' helper enumerates those configurations once so expected margins can reuse
#' them during each fit iteration.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param context Prepared GLLRM computation context.
#' @param component_items Internal `component_items` value used by this helper.
#' @return The internal `gllrm_component_configurations()` computation result.
#' @keywords internal
#' @noRd
gllrm_component_configurations <- function(context, component_items) {
  key <- gllrm_component_key(component_items)
  if (!is.null(context$component_configurations) &&
      !is.null(context$component_configurations[[key]])) {
    return(context$component_configurations[[key]])
  }
  gllrm_build_component_configurations(context, component_items)
}

#' Internal gllrm build component configurations helper
#'
#' Supports the gllrm components implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param context Prepared GLLRM computation context.
#' @param component_items Internal `component_items` value used by this helper.
#' @return The internal `gllrm_build_component_configurations()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal gllrm component config weight helper
#'
#' Supports the gllrm components implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param component_items Internal `component_items` value used by this helper.
#' @param item_values Internal `item_values` value used by this helper.
#' @param background_values Internal `background_values` value used by this helper.
#' @return The internal `gllrm_component_config_weight()` computation result.
#' @keywords internal
#' @noRd
gllrm_component_config_weight <- function(context, state, component_items, item_values, background_values) {
  weight <- 1
  dif_map <- context$dif_map %||% gllrm_dif_map(context)
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

#' Internal gllrm component config weight fast helper
#'
#' Supports the gllrm components implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param component_items Internal `component_items` value used by this helper.
#' @param item_values Internal `item_values` value used by this helper.
#' @param background_values Internal `background_values` value used by this helper.
#' @param key Internal `key` value used by this helper.
#' @return The internal `gllrm_component_config_weight_fast()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal gllrm component config weights fast helper
#'
#' Supports the gllrm components implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param component_items Internal `component_items` value used by this helper.
#' @param config_matrix Internal `config_matrix` value used by this helper.
#' @param background_values Internal `background_values` value used by this helper.
#' @param key Internal `key` value used by this helper.
#' @return The internal `gllrm_component_config_weights_fast()` computation result.
#' @keywords internal
#' @noRd
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

#' Source trace: source/digram_source_20260817/skunits/skbias22.pas::CalculateBiasedGammaValues2 multiplies
#' item score parameters, included IJ parameters, and included IX parameters before
#' aggregating by total score. This helper computes that same unnormalized
#' component contribution for one LD-connected component and one background
#' combination.
#' Source trace: `source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param component_items Internal `component_items` value used by this helper.
#' @param background_values Internal `background_values` value used by this helper.
#' @return The internal `gllrm_component_gamma()` computation result.
#' @keywords internal
#' @noRd
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

#' Source trace: source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation rescales
#' score polynomials for numerical stability without changing the fitted
#' probabilities. The R normalization is the same implementation-level
#' stabilization: only ratios are used downstream.
#' @param gamma Internal `gamma` value used by this helper.
#' @return The internal `gllrm_normalize_component_gamma()` computation result.
#' @keywords internal
#' @noRd
gllrm_normalize_component_gamma <- function(gamma) {
  scale <- max(gamma)
  if (scale <= 0) {
    scale <- 1
  }
  list(gamma = gamma / scale, scale = scale)
}
