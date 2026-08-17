# M2/M3 value orchestration.

#' Internal m2 values helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param fit Fitted gRm model.
#' @param items Item selection or item metadata.
#' @param score_cuts Resolved total-score cut values.
#' @param bootstrap_control Internal `bootstrap_control` value used by this helper.
#' @return The internal `m2_values()` computation result.
#' @keywords internal
#' @noRd
m2_values <- function(fit,
                      items = NULL,
                      score_cuts = NULL,
                      bootstrap_control = list(enabled = FALSE)) {
  m2_m3_values(
    fit,
    items = items,
    score_cuts = score_cuts,
    include_three_way = FALSE,
    bootstrap_control = bootstrap_control
  )
}

#' Internal m3 values helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param fit Fitted gRm model.
#' @param items Item selection or item metadata.
#' @param score_cuts Resolved total-score cut values.
#' @param bootstrap_control Internal `bootstrap_control` value used by this helper.
#' @return The internal `m3_values()` computation result.
#' @keywords internal
#' @noRd
m3_values <- function(fit,
                      items = NULL,
                      score_cuts = NULL,
                      bootstrap_control = list(enabled = FALSE)) {
  m2_m3_values(
    fit,
    items = items,
    score_cuts = score_cuts,
    include_three_way = TRUE,
    bootstrap_control = bootstrap_control
  )
}

#' Internal m2 m3 values helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param fit Fitted gRm model.
#' @param items Item selection or item metadata.
#' @param score_cuts Resolved total-score cut values.
#' @param include_three_way Whether to include three-way CM3 margins.
#' @param bootstrap_control Internal `bootstrap_control` value used by this helper.
#' @return The internal `m2_m3_values()` computation result.
#' @keywords internal
#' @noRd
m2_m3_values <- function(fit,
                         items = NULL,
                         score_cuts = NULL,
                         include_three_way,
                         bootstrap_control = list(enabled = FALSE)) {
  # Source trace: source/PAS_skunits/skbias14.pas::CM3_analysis
  # (7289-7587) uses the current fitted GLLRM, prepares one margin list, and
  # runs CM3_tests once; R shares this engine for m2() and m3().
  fit <- as_public_gRm_fit(fit)
  analysis <- fit$analysis %||% fit$spec$analysis
  context_state <- m2_m3_fit_context_state(fit)
  context <- context_state$context
  state <- context_state$state
  score_cuts <- normalize_public_score_cuts(
    score_cuts,
    analysis$project,
    default = analysis$score_groups,
    bundle = fit$bundle %||% context$bundle %||% NULL
  )
  selected_indices <- m2_m3_selected_items(context, items)
  specs <- m2_m3_prepare_margins(context, selected_indices, include_three_way = include_three_way)
  score_groups <- global_homogeneity_score_groups(context$bundle, score_cuts)
  score_group_lookup <- m2_m3_score_group_lookup(context, score_cuts)
  probability_cache <- new_m2_m3_focus_probability_cache(context, state)
  tests <- if (length(specs)) {
    do.call(rbind, lapply(specs, function(spec) {
      m2_m3_analyze_margin(context, state, spec, score_group_lookup, probability_cache)
    }))
  } else {
    data.frame()
  }
  rownames(tests) <- NULL
  aggregates <- m2_m3_aggregate_values(tests)

  if (isTRUE(include_three_way)) {
    out <- list(
      tests = tests,
      m2 = aggregates$m2,
      m3 = aggregates$m3,
      item_trait = aggregates$item_trait,
      invariance = aggregates$invariance,
      margin_specs = specs,
      selected_items = m2_m3_selected_item_table(context, selected_indices),
      exogenous_names = as.character(m2_m3_context_backgrounds(context)$name),
      score_cuts = score_cuts,
      score_groups = score_groups,
      n_two_way_margins = sum(tests$is_m2 %||% logical()),
      n_three_way_margins = sum(!(tests$is_m2 %||% logical())),
      source_status = "m3_values"
    )
    out <- m2_m3_apply_parametric_bootstrap(
      fit,
      out,
      context,
      state,
      include_three_way = TRUE,
      control = bootstrap_control
    )
    class(out) <- c("gRm_m3_values", "list")
    out
  } else {
    out <- list(
      tests = tests,
      aggregate = aggregates$m2,
      item_trait = aggregates$item_trait,
      invariance = aggregates$invariance,
      margin_specs = specs,
      selected_items = m2_m3_selected_item_table(context, selected_indices),
      exogenous_names = as.character(m2_m3_context_backgrounds(context)$name),
      score_cuts = score_cuts,
      score_groups = score_groups,
      n_two_way_margins = sum(tests$is_m2 %||% logical()),
      n_three_way_margins = 0L,
      source_status = "m2_values"
    )
    out <- m2_m3_apply_parametric_bootstrap(
      fit,
      out,
      context,
      state,
      include_three_way = FALSE,
      control = bootstrap_control
    )
    class(out) <- c("gRm_m2_values", "list")
    out
  }
}

#' Internal m2 m3 fit context state helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param fit Fitted gRm model.
#' @return The internal `m2_m3_fit_context_state()` computation result.
#' @keywords internal
#' @noRd
m2_m3_fit_context_state <- function(fit) {
  if (!is.null(fit$fit$context)) {
    return(list(context = fit$fit$context, state = fit$fit))
  }

  bundle <- fit$bundle %||% build_item_parameters_bundle((fit$analysis %||% fit$spec$analysis)$project)
  context <- build_gllrm_context(fit$spec, bundle)
  state <- initialize_gllrm_state(context)
  state$item_gamma <- fit$fit$item_gamma
  state$expected_items <- fit$fit$expected_items %||% state$expected_items
  state$counts <- fit$fit$counts %||% context$counts
  state$converged <- fit$fit$converged %||% TRUE
  state$n_step <- fit$fit$n_step %||% NA_integer_
  state$delta <- fit$fit$delta %||% NA_real_
  list(context = context, state = state)
}

#' Internal m2 m3 selected item table helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param selected_indices Internal `selected_indices` value used by this helper.
#' @return The internal `m2_m3_selected_item_table()` computation result.
#' @keywords internal
#' @noRd
m2_m3_selected_item_table <- function(context, selected_indices) {
  items <- m2_m3_context_items(context)
  labels <- m2_m3_variable_labels(items)
  data.frame(
    item_index = as.integer(selected_indices),
    item_label = labels[selected_indices],
    item_name = as.character(items$name[selected_indices]),
    stringsAsFactors = FALSE
  )
}

#' Internal m2 m3 analyze margin helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param spec GLLRM model specification.
#' @param score_group_lookup Internal `score_group_lookup` value used by this helper.
#' @param probability_cache Internal `probability_cache` value used by this helper.
#' @return The internal `m2_m3_analyze_margin()` computation result.
#' @keywords internal
#' @noRd
m2_m3_analyze_margin <- function(context, state, spec, score_group_lookup, probability_cache) {
  # Source trace: source/PAS_skunits/skbias14.pas::Twoway_analysis and
  # ::Threeway_analysis dispatch observed and expected tables, then compute the
  # Pearson contribution for each prepared CM2/CM3 margin.
  observed <- m2_m3_count_observed(context, spec, score_group_lookup)
  expected <- m2_m3_expected_table(context, state, spec, score_group_lookup, probability_cache)
  dimensions <- dim(observed)
  df <- if (length(dimensions) == 2L) {
    m2_m3_df_two_way(dimensions)
  } else {
    m2_m3_df_three_way(dimensions)
  }
  stat <- m2_m3_pearson_stat(observed, expected, df = df)
  cbind(
    m2_m3_margin_metadata(context, spec),
    data.frame(
      chi_square = stat$chi_square,
      degrees_of_freedom = stat$degrees_of_freedom,
      p_value = stat$p_value,
      stringsAsFactors = FALSE
    )
  )
}

#' Internal m2 m3 pearson stat helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param observed Internal `observed` value used by this helper.
#' @param expected Internal `expected` value used by this helper.
#' @param df Internal `df` value used by this helper.
#' @return The internal `m2_m3_pearson_stat()` computation result.
#' @keywords internal
#' @noRd
m2_m3_pearson_stat <- function(observed, expected, df) {
  # Source trace: source/PAS_skunits/skbias14.pas::Twoway_analysis and
  # ::Threeway_analysis add (observed - expected)^2 / expected only for cells
  # with expected count greater than zero.
  keep <- expected > 0
  chi_square <- if (any(keep)) {
    sum((observed[keep] - expected[keep])^2 / expected[keep])
  } else {
    0
  }
  df <- as.integer(df)
  # Source trace: source/PAS_skunits/skbias14.pas calls PFCHI for the printed
  # CM2/CM3 row and aggregate p-values.
  p_value <- if (df > 0L) source_pfchi(df, chi_square) else 1
  list(
    chi_square = as.numeric(chi_square),
    degrees_of_freedom = df,
    p_value = as.numeric(p_value)
  )
}

#' Internal m2 m3 df two way helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param dimensions Internal `dimensions` value used by this helper.
#' @return The internal `m2_m3_df_two_way()` computation result.
#' @keywords internal
#' @noRd
m2_m3_df_two_way <- function(dimensions) {
  # Source trace: source/PAS_skunits/skbias14.pas::Twoway_analysis uses the
  # ordinary two-way margin df formula (columns - 1) * (rows - 1).
  dimensions <- as.integer(dimensions)
  as.integer((dimensions[[1L]] - 1L) * (dimensions[[2L]] - 1L))
}

#' Internal m2 m3 df three way helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param dimensions Internal `dimensions` value used by this helper.
#' @return The internal `m2_m3_df_three_way()` computation result.
#' @keywords internal
#' @noRd
m2_m3_df_three_way <- function(dimensions) {
  # Source trace: source/PAS_skunits/skbias14.pas::Threeway_analysis subtracts
  # the three one-way margins and the grand total from the full table cells.
  dimensions <- as.integer(dimensions)
  as.integer(prod(dimensions) - 1L - sum(dimensions - 1L))
}

#' Internal m2 m3 aggregate values helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param tests Diagnostic test rows.
#' @param context Prepared GLLRM computation context.
#' @return The internal `m2_m3_aggregate_values()` computation result.
#' @keywords internal
#' @noRd
m2_m3_aggregate_values <- function(tests, context = NULL) {
  # Source trace: source/PAS_skunits/skbias14.pas:7181-7211 stores the CM2 aggregate
  # after the first nM2tests rows, then initializes the CM3 aggregate with that
  # complete CM2 total before adding three-way rows.
  tests <- tests %||% data.frame()
  m2_rows <- if (is.data.frame(tests) && nrow(tests)) tests[tests$is_m2, , drop = FALSE] else data.frame()
  m3_rows <- tests
  item_trait_rows <- if (nrow(m2_rows)) {
    m2_rows[m2_rows$margin_type == "item_score_group", , drop = FALSE]
  } else {
    data.frame()
  }
  item_exogenous_rows <- if (nrow(m2_rows)) {
    m2_rows[m2_rows$margin_type == "item_exogenous", , drop = FALSE]
  } else {
    data.frame()
  }
  list(
    m2 = m2_m3_aggregate_row("M2", m2_rows),
    m3 = m2_m3_aggregate_row("M3", m3_rows),
    item_trait = m2_m3_aggregate_row("Item-trait", item_trait_rows),
    invariance = m2_m3_invariance_rows(item_exogenous_rows)
  )
}

#' Internal m2 m3 aggregate row helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param diagnostic Internal `diagnostic` value used by this helper.
#' @param rows Rows used by the computation.
#' @return The internal `m2_m3_aggregate_row()` computation result.
#' @keywords internal
#' @noRd
m2_m3_aggregate_row <- function(diagnostic, rows) {
  chi_square <- if (is.data.frame(rows) && nrow(rows)) sum(rows$chi_square) else 0
  degrees_of_freedom <- if (is.data.frame(rows) && nrow(rows)) sum(rows$degrees_of_freedom) else 0L
  p_value <- if (degrees_of_freedom > 0L) source_pfchi(degrees_of_freedom, chi_square) else 1
  data.frame(
    diagnostic = diagnostic,
    chi_square = as.numeric(chi_square),
    degrees_of_freedom = as.integer(degrees_of_freedom),
    p_value = as.numeric(p_value),
    stringsAsFactors = FALSE
  )
}

#' Internal m2 m3 invariance rows helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param rows Rows used by the computation.
#' @return The internal `m2_m3_invariance_rows()` computation result.
#' @keywords internal
#' @noRd
m2_m3_invariance_rows <- function(rows) {
  if (!is.data.frame(rows) || !nrow(rows)) {
    return(data.frame(
      background_label = character(),
      background_name = character(),
      chi_square = numeric(),
      degrees_of_freedom = integer(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  keys <- unique(rows$background_label)
  out <- lapply(keys, function(key) {
    group <- rows[rows$background_label == key, , drop = FALSE]
    aggregate <- m2_m3_aggregate_row(group$background_name[[1L]], group)
    data.frame(
      background_label = key,
      background_name = group$background_name[[1L]],
      chi_square = aggregate$chi_square,
      degrees_of_freedom = aggregate$degrees_of_freedom,
      p_value = aggregate$p_value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

#' Internal m2 m3 margin metadata helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @return The internal `m2_m3_margin_metadata()` computation result.
#' @keywords internal
#' @noRd
m2_m3_margin_metadata <- function(context, spec) {
  items <- m2_m3_context_items(context)
  backgrounds <- m2_m3_context_backgrounds(context)
  item_labels <- m2_m3_variable_labels(items)
  background_labels <- m2_m3_variable_labels(backgrounds)
  margin_variables <- m2_m3_margin_public_variables(context, spec)
  empty <- NA_character_
  row <- data.frame(
    margin = m2_m3_margin_name(context, spec),
    margin_label = m2_m3_margin_label(context, spec),
    margin_type = spec$kind,
    is_m2 = spec$kind %in% c("item_item", "item_exogenous", "item_score_group"),
    item_label = empty,
    item_name = empty,
    item1_label = empty,
    item1_name = empty,
    item2_label = empty,
    item2_name = empty,
    item3_label = empty,
    item3_name = empty,
    background_label = empty,
    background_name = empty,
    background1_label = empty,
    background1_name = empty,
    background2_label = empty,
    background2_name = empty,
    score_label = empty,
    score_name = empty,
    stringsAsFactors = FALSE
  )

  set_item <- function(prefix, item) {
    row[[paste0(prefix, "_label")]] <<- item_labels[[item]]
    row[[paste0(prefix, "_name")]] <<- as.character(items$name[[item]])
  }
  set_background <- function(prefix, background) {
    row[[paste0(prefix, "_label")]] <<- background_labels[[background]]
    row[[paste0(prefix, "_name")]] <<- as.character(backgrounds$name[[background]])
  }
  set_score <- function() {
    row$score_label <<- "#"
    row$score_name <<- "Score group"
  }

  switch(
    spec$kind,
    item_item = {
      set_item("item1", spec$item1)
      set_item("item2", spec$item2)
    },
    item_exogenous = {
      set_item("item", spec$item)
      set_item("item1", spec$item)
      set_background("background", spec$exogenous)
    },
    item_score_group = {
      set_item("item", spec$item)
      set_item("item1", spec$item)
      set_score()
    },
    item_item_item = {
      set_item("item1", spec$item1)
      set_item("item2", spec$item2)
      set_item("item3", spec$item3)
    },
    item_item_exogenous = {
      set_item("item1", spec$item1)
      set_item("item2", spec$item2)
      set_background("background", spec$exogenous)
    },
    item_item_score_group = {
      set_item("item1", spec$item1)
      set_item("item2", spec$item2)
      set_score()
    },
    item_exogenous_exogenous = {
      set_item("item", spec$item)
      set_item("item1", spec$item)
      set_background("background1", spec$exogenous1)
      set_background("background2", spec$exogenous2)
      row$background_label <- row$background1_label
      row$background_name <- row$background1_name
    },
    item_exogenous_score_group = {
      set_item("item", spec$item)
      set_item("item1", spec$item)
      set_background("background", spec$exogenous)
      set_score()
    },
    stop("Unknown M2/M3 margin kind.", call. = FALSE)
  )
  row$variable1 <- m2_m3_nth_character(margin_variables, 1L, empty)
  row$variable2 <- m2_m3_nth_character(margin_variables, 2L, empty)
  row$variable3 <- m2_m3_nth_character(margin_variables, 3L, empty)
  row
}

#' Internal m2 m3 public margin table helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param tests Diagnostic test rows.
#' @param context Prepared GLLRM computation context.
#' @return The internal `m2_m3_public_margin_table()` computation result.
#' @keywords internal
#' @noRd
m2_m3_public_margin_table <- function(tests, context = NULL) {
  tests
}

#' Internal m2 m3 nth character helper
#'
#' Supports the m2 m3 values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias14.pas::CM3_analysis`.
#' @param x Object or value to process.
#' @param index One-based internal index.
#' @param default Internal `default` value used by this helper.
#' @return The internal `m2_m3_nth_character()` computation result.
#' @keywords internal
#' @noRd
m2_m3_nth_character <- function(x, index, default = NA_character_) {
  if (length(x) < index || is.na(x[[index]])) {
    return(default)
  }
  as.character(x[[index]])
}
