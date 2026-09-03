#' Internal gllrm global homogeneity uniform dif helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param fit Fitted gRm model.
#' @param groups Internal `groups` value used by this helper.
#' @param probability_cache Internal `probability_cache` value used by this helper.
#' @return The internal `gllrm_global_homogeneity_uniform_dif()` computation result.
#' @keywords internal
#' @noRd
gllrm_global_homogeneity_uniform_dif <- function(fit, groups, probability_cache = NULL) {
  context <- fit$fit$context
  state <- fit$fit
  if (length(context$dif_specs) == 0L) {
    return(data.frame())
  }

  rows <- vector("list", length(context$dif_specs))
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  probability_cache <- probability_cache %||%
    new_gllrm_probability_cache(context, state, components = components)
  tables_by_dif <- gllrm_uniform_dif_scoregroup_tables_all(context, groups, probability_cache)
  for (dif_index in seq_along(context$dif_specs)) {
    spec <- context$dif_specs[[dif_index]]
    rows[[dif_index]] <- gllrm_uniform_dif_summary_row(context, groups, spec, tables_by_dif[[dif_index]])
  }
  do.call(rbind, rows)
}

#' Internal gllrm uniform dif scoregroup tables helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param context Prepared GLLRM computation context.
#' @param groups Internal `groups` value used by this helper.
#' @param spec GLLRM model specification.
#' @param probability_cache Internal `probability_cache` value used by this helper.
#' @return The internal `gllrm_uniform_dif_scoregroup_tables()` computation result.
#' @keywords internal
#' @noRd
gllrm_uniform_dif_scoregroup_tables <- function(context, groups, spec, probability_cache) {
  gllrm_uniform_dif_scoregroup_tables_for_specs(
    context,
    groups,
    list(spec),
    probability_cache
  )[[1L]]
}

#' Internal gllrm uniform dif scoregroup tables all helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param context Prepared GLLRM computation context.
#' @param groups Internal `groups` value used by this helper.
#' @param probability_cache Internal `probability_cache` value used by this helper.
#' @return The internal `gllrm_uniform_dif_scoregroup_tables_all()` computation result.
#' @keywords internal
#' @noRd
gllrm_uniform_dif_scoregroup_tables_all <- function(context, groups, probability_cache) {
  gllrm_uniform_dif_scoregroup_tables_for_specs(
    context,
    groups,
    context$dif_specs,
    probability_cache
  )
}

#' Internal gllrm uniform dif scoregroup tables for specs helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param context Prepared GLLRM computation context.
#' @param groups Internal `groups` value used by this helper.
#' @param specs Internal `specs` value used by this helper.
#' @param probability_cache Internal `probability_cache` value used by this helper.
#' @return The internal `gllrm_uniform_dif_scoregroup_tables_for_specs()` computation result.
#' @keywords internal
#' @noRd
gllrm_uniform_dif_scoregroup_tables_for_specs <- function(context, groups, specs, probability_cache) {
  tables_by_dif <- lapply(specs, function(spec) {
    item_max <- context$item_raw_max[[spec$item]] - 1L
    background_max <- context$background_raw_max[[spec$background]]
    dimensions <- c(item_max + 1L, background_max, nrow(groups))
    list(
      observed = array(0, dim = dimensions),
      expected = array(0, dim = dimensions)
    )
  })
  if (length(tables_by_dif) == 0L) {
    return(tables_by_dif)
  }

  score_group_lookup <- global_homogeneity_uniform_score_group_lookup(groups, context$max_total_score)

  uniform_rows <- gllrm_uniform_complete_rows(context, score_group_lookup)
  for (row in uniform_rows) {
    group_index <- global_homogeneity_lookup_score(score_group_lookup, context$score[[row]])
    if (is.na(group_index)) {
      next
    }
    for (dif_index in seq_along(specs)) {
      spec <- specs[[dif_index]]
      item_score <- context$item_matrix[row, spec$item] + 1L
      background_value <- context$background_matrix[row, spec$background]
      tables_by_dif[[dif_index]]$observed[item_score, background_value, group_index] <-
        tables_by_dif[[dif_index]]$observed[item_score, background_value, group_index] + 1
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
    probabilities_by_item <- gllrm_cached_item_probabilities(
      probability_cache,
      total_score = score,
      background_values = background_values
    )
    for (dif_index in seq_along(specs)) {
      spec <- specs[[dif_index]]
      background_value <- background_values[[spec$background]]
      probabilities <- probabilities_by_item[[spec$item]]
      tables_by_dif[[dif_index]]$expected[seq_along(probabilities), background_value, homogeneity_group] <-
        tables_by_dif[[dif_index]]$expected[seq_along(probabilities), background_value, homogeneity_group] +
          group$count[[1L]] * probabilities
    }
  }

  tables_by_dif
}

#' Internal gllrm uniform dif summary row helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param context Prepared GLLRM computation context.
#' @param groups Internal `groups` value used by this helper.
#' @param spec GLLRM model specification.
#' @param tables Internal `tables` value used by this helper.
#' @return The internal `gllrm_uniform_dif_summary_row()` computation result.
#' @keywords internal
#' @noRd
gllrm_uniform_dif_summary_row <- function(context, groups, spec, tables) {
  stats <- gllrm_uniform_summary_stats(tables, nrow(groups))

  data.frame(
    item = spec$item,
    background = spec$background,
    item_label = context$items$label_code[[spec$item]],
    background_label = context$backgrounds$label_code[[spec$background]],
    observed_gamma = I(list(stats$observed_gamma)),
    expected_gamma = I(list(stats$expected_gamma)),
    chi_square = stats$chi_square,
    df = stats$df,
    p_value = stats$p_value,
    stringsAsFactors = FALSE
  )
}

#' Internal gllrm uniform dif standardize expected helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param observed Internal `observed` value used by this helper.
#' @param expected Internal `expected` value used by this helper.
#' @return The internal `gllrm_uniform_dif_standardize_expected()` computation result.
#' @keywords internal
#' @noRd
gllrm_uniform_dif_standardize_expected <- function(observed, expected) {
  # Source trace: source/digram_source_20260817/skunits/skfit2.pas::Standardize_ETAB2_to_TAB2_margins
  # derives observed margins, then delegates to
  # source/digram_source_20260817/skunits/skfit2.pas::Standardize_tab4. The R table stores only
  # body cells, so the shared helper applies the same fixed row/column scaling
  # to the expected table before the global-homogeneity gamma is calculated.
  row_margins <- rowSums(observed)
  col_margins <- colSums(observed)
  out <- source_standardize_table_margins(expected, row_margins, col_margins)
  list(
    table = out,
    df = max(0L, (sum(row_margins > 0) - 1L) * (sum(col_margins > 0) - 1L))
  )
}
