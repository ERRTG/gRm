#' Internal gllrm uniform ld scoregroup tables helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param context Prepared GLLRM computation context.
#' @param groups Internal `groups` value used by this helper.
#' @param spec GLLRM model specification.
#' @param ld_index Internal `ld_index` value used by this helper.
#' @param probability_cache Internal `probability_cache` value used by this helper.
#' @return The internal `gllrm_uniform_ld_scoregroup_tables()` computation result.
#' @keywords internal
#' @noRd
gllrm_uniform_ld_scoregroup_tables <- function(context, groups, spec, ld_index, probability_cache) {
  gllrm_uniform_ld_scoregroup_tables_all(context, groups, probability_cache)[[ld_index]]
}

#' Internal gllrm uniform ld scoregroup tables all helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param context Prepared GLLRM computation context.
#' @param groups Internal `groups` value used by this helper.
#' @param probability_cache Internal `probability_cache` value used by this helper.
#' @return The internal `gllrm_uniform_ld_scoregroup_tables_all()` computation result.
#' @keywords internal
#' @noRd
gllrm_uniform_ld_scoregroup_tables_all <- function(context, groups, probability_cache) {
  tables_by_ld <- lapply(context$ld_specs, function(spec) {
    dimensions <- c(context$item_raw_max[[spec$item1]], context$item_raw_max[[spec$item2]], nrow(groups))
    list(
      observed = array(0, dim = dimensions),
      expected = array(0, dim = dimensions)
    )
  })
  if (length(tables_by_ld) == 0L) {
    return(tables_by_ld)
  }

  score_group_lookup <- global_homogeneity_uniform_score_group_lookup(groups, context$max_total_score)
  uniform_rows <- gllrm_uniform_complete_rows(context, score_group_lookup)
  for (row in uniform_rows) {
    group_index <- global_homogeneity_lookup_score(score_group_lookup, context$score[[row]])
    if (is.na(group_index)) {
      next
    }
    for (ld_index in seq_along(context$ld_specs)) {
      spec <- context$ld_specs[[ld_index]]
      score1 <- context$item_matrix[row, spec$item1] + 1L
      score2 <- context$item_matrix[row, spec$item2] + 1L
      tables_by_ld[[ld_index]]$observed[score1, score2, group_index] <-
        tables_by_ld[[ld_index]]$observed[score1, score2, group_index] + 1
    }
  }

  score_exo_groups <- gllrm_score_exo_groups(context, rows = uniform_rows)
  for (group_index in seq_len(nrow(score_exo_groups))) {
    group <- score_exo_groups[group_index, , drop = FALSE]
    score <- group$score[[1L]]
    homogeneity_group <- global_homogeneity_lookup_score(score_group_lookup, score)
    if (is.na(homogeneity_group)) {
      next
    }
    background_values <- gllrm_group_background_values(context, group)
    probabilities_by_ld <- gllrm_cached_ld_probabilities(
      probability_cache,
      total_score = score,
      background_values = background_values
    )
    for (ld_index in seq_along(context$ld_specs)) {
      tables_by_ld[[ld_index]]$expected[, , homogeneity_group] <-
        tables_by_ld[[ld_index]]$expected[, , homogeneity_group] +
        group$count[[1L]] * probabilities_by_ld[[ld_index]]
    }
  }

  tables_by_ld
}

#' Internal gllrm uniform summary stats helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param tables Internal `tables` value used by this helper.
#' @param n_groups Internal `n_groups` value used by this helper.
#' @return The internal `gllrm_uniform_summary_stats()` computation result.
#' @keywords internal
#' @noRd
gllrm_uniform_summary_stats <- function(tables, n_groups) {
  observed_gamma <- numeric(n_groups)
  expected_gamma <- numeric(n_groups)
  chi_square <- 0
  df <- 0L
  x_dim <- dim(tables$observed)[[1L]]
  y_dim <- dim(tables$observed)[[2L]]

  for (group_index in seq_len(n_groups)) {
    observed <- tables$observed[, , group_index, drop = FALSE][, , 1L]
    expected <- tables$expected[, , group_index, drop = FALSE][, , 1L]
    positive <- expected > 0
    if (any(positive)) {
      chi_square <- chi_square + sum((observed[positive] - expected[positive])^2 / expected[positive])
    }
    standardized <- gllrm_uniform_dif_standardize_expected(observed, expected)
    observed_gamma[[group_index]] <- source_gamma_from_table(observed)
    expected_gamma[[group_index]] <- source_gamma_from_table(standardized$table)
    df <- df + standardized$df
  }
  df <- df - (x_dim - 1L) * (y_dim - 1L)
  # Source-faithfulness guardrail: this uniform degrees-of-freedom floor
  # is intentionally left unchanged until source or validator evidence supports
  # a different convention. It was identified as mathematically unusual in
  # docs/GRM_R_PACKAGE_READ_ONLY_CODE_AUDIT_2026-06-12.md, so future edits
  # should not "simplify" it without a source-backed comparison; any change
  # requires source or validator evidence before changing the convention.
  if (df < 1L) {
    df <- 1L
  }

  list(
    observed_gamma = observed_gamma,
    expected_gamma = expected_gamma,
    chi_square = chi_square,
    df = as.integer(df),
    p_value = source_pfchi(df, chi_square)
  )
}

#' Internal gllrm uniform ld summary row helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param context Prepared GLLRM computation context.
#' @param groups Internal `groups` value used by this helper.
#' @param spec GLLRM model specification.
#' @param tables Internal `tables` value used by this helper.
#' @return The internal `gllrm_uniform_ld_summary_row()` computation result.
#' @keywords internal
#' @noRd
gllrm_uniform_ld_summary_row <- function(context, groups, spec, tables) {
  stats <- gllrm_uniform_summary_stats(tables, nrow(groups))

  data.frame(
    item1 = spec$item1,
    item2 = spec$item2,
    item1_label = context$items$label_code[[spec$item1]],
    item2_label = context$items$label_code[[spec$item2]],
    observed_gamma = I(list(stats$observed_gamma)),
    expected_gamma = I(list(stats$expected_gamma)),
    chi_square = stats$chi_square,
    df = stats$df,
    p_value = stats$p_value,
    stringsAsFactors = FALSE
  )
}

#' Internal gllrm group ld probabilities helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param total_score Internal `total_score` value used by this helper.
#' @param background_values Internal `background_values` value used by this helper.
#' @param components Internal `components` value used by this helper.
#' @return The internal `gllrm_group_ld_probabilities()` computation result.
#' @keywords internal
#' @noRd
gllrm_group_ld_probabilities <- function(context,
                                                state,
                                                total_score,
                                                background_values,
                                                components = NULL) {
  components <- components %||% context$ld_components_items %||% gllrm_ld_components(context)$items
  out <- lapply(context$ld_specs, function(spec) {
    matrix(0, nrow = context$item_raw_max[[spec$item1]], ncol = context$item_raw_max[[spec$item2]])
  })
  if (total_score < 0L || total_score > context$max_total_score) {
    return(out)
  }

  component_gamma <- lapply(components, function(component_items) {
    one <- gllrm_component_gamma(context, state, component_items, background_values)
    one$gamma * one$scale
  })
  convolutions <- gllrm_component_convolutions(component_gamma, context$max_total_score)
  denominator <- convolutions$full[[total_score + 1L]]
  if (denominator <= 0) {
    return(out)
  }

  for (component_index in seq_along(components)) {
    component_items <- components[[component_index]]
    key <- gllrm_component_key(component_items)
    ld_local <- context$component_ld_local_matrices[[key]] %||%
      context$component_ld_local_indices[[key]]
    if (is.null(ld_local) || nrow(ld_local) == 0L) {
      next
    }
    rest_gamma <- convolutions$rest[[component_index]]
    configs <- context$component_config_matrices[[key]]
    config_scores <- context$component_config_scores[[key]]
    weights <- gllrm_component_config_weights_fast(
      context,
      state,
      component_items,
      configs,
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

    for (ld_row in seq_len(nrow(ld_local))) {
      ld_index <- ld_local[ld_row, 1L]
      score1 <- configs[valid_index, ld_local[ld_row, 2L]] + 1L
      score2 <- configs[valid_index, ld_local[ld_row, 3L]] + 1L
      n_score1 <- nrow(out[[ld_index]])
      pair_index <- score1 + (score2 - 1L) * n_score1
      by_pair <- rowsum(probabilities, pair_index, reorder = FALSE)
      out[[ld_index]][as.integer(rownames(by_pair))] <-
        out[[ld_index]][as.integer(rownames(by_pair))] + as.numeric(by_pair[, 1L])
    }
  }

  out
}
