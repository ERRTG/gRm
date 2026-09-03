# CM2/CM3 value orchestration.

#' Internal cm2 values helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param fit Fitted gRm model.
#' @param items Item selection or item metadata.
#' @param score_cuts Resolved total-score cut values.
#' @param bootstrap_control Internal `bootstrap_control` value used by this helper.
#' @return The internal `cm2_values()` computation result.
#' @keywords internal
#' @noRd
cm2_values <- function(fit,
                      items = NULL,
                      score_cuts = NULL,
                      bootstrap_control = list(enabled = FALSE)) {
  cm2_cm3_values(
    fit,
    items = items,
    score_cuts = score_cuts,
    include_three_way = FALSE,
    bootstrap_control = bootstrap_control
  )
}

#' Internal cm3 values helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param fit Fitted gRm model.
#' @param items Item selection or item metadata.
#' @param score_cuts Resolved total-score cut values.
#' @param bootstrap_control Internal `bootstrap_control` value used by this helper.
#' @return The internal `cm3_values()` computation result.
#' @keywords internal
#' @noRd
cm3_values <- function(fit,
                      items = NULL,
                      score_cuts = NULL,
                      bootstrap_control = list(enabled = FALSE)) {
  cm2_cm3_values(
    fit,
    items = items,
    score_cuts = score_cuts,
    include_three_way = TRUE,
    bootstrap_control = bootstrap_control
  )
}

#' Internal cm2 cm3 values helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param fit Fitted gRm model.
#' @param items Item selection or item metadata.
#' @param score_cuts Resolved total-score cut values.
#' @param include_three_way Whether to include three-way CM3 margins.
#' @param bootstrap_control Internal `bootstrap_control` value used by this helper.
#' @return The internal `cm2_cm3_values()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_values <- function(fit,
                         items = NULL,
                         score_cuts = NULL,
                         include_three_way,
                         bootstrap_control = list(enabled = FALSE)) {
  # Source trace: source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis
  # The canonical routine uses the current fitted GLLRM, prepares one margin list, and
  # runs CM3_tests once; R shares this engine for cm2() and cm3().
  fit <- as_public_gRm_fit(fit)
  analysis <- fit$analysis %||% fit$spec$analysis
  context_state <- cm2_cm3_fit_context_state(fit)
  context <- context_state$context
  state <- context_state$state
  score_cuts <- normalize_public_score_cuts(
    score_cuts,
    analysis$project,
    default = analysis$score_groups,
    bundle = fit$bundle %||% context$bundle %||% NULL
  )
  # SKbias2.SelectItems resolves the focal Boolean set once; CM3_analysis then
  # passes that set to Prepare_CM3tests without altering the fitted GLLRM.
  # `selection` additionally records the evaluated public request, while only
  # its source-ordered resolved indices enter the numerical margin loops.
  selection <- cm2_cm3_selected_items(context, items)
  selected_indices <- selection$resolved_items$item_index
  specs <- cm2_cm3_prepare_margins(context, selected_indices, include_three_way = include_three_way)
  score_groups <- global_homogeneity_score_groups(context$bundle, score_cuts)
  score_group_lookup <- cm2_cm3_score_group_lookup(context, score_cuts)
  probability_cache <- new_cm2_cm3_focus_probability_cache(context, state)
  tests <- if (length(specs)) {
    do.call(rbind, lapply(specs, function(spec) {
      cm2_cm3_analyze_margin(context, state, spec, score_group_lookup, probability_cache)
    }))
  } else {
    data.frame()
  }
  rownames(tests) <- NULL
  # CM3_tests initializes invariance over the complete fitted exogenous vector,
  # not merely the exogenous labels still represented by eligible DIF-filtered
  # rows, so aggregation must retain the full computation context.
  aggregates <- cm2_cm3_aggregate_values(tests, context)

  cm2_rows <- if (nrow(tests)) tests[tests$is_cm2, , drop = FALSE] else data.frame()
  cm3_rows <- if (nrow(tests)) tests[!tests$is_cm2, , drop = FALSE] else data.frame()
  cm2_bh <- cm2_cm3_bh_thresholds(cm2_rows, "CM2")
  cm3_bh <- cm2_cm3_bh_thresholds(cm3_rows, "CM3")

  if (isTRUE(include_three_way)) {
    out <- list(
      tests = tests,
      cm2 = aggregates$cm2,
      cm3 = aggregates$cm3,
      item_trait = aggregates$item_trait,
      invariance = aggregates$invariance,
      cm2_bh = cm2_bh,
      cm3_bh = cm3_bh,
      margin_specs = specs,
      selected_items = selection$resolved_items,
      selection = selection,
      exogenous_names = as.character(cm2_cm3_context_backgrounds(context)$name),
      score_cuts = score_cuts,
      score_groups = score_groups,
      n_two_way_margins = sum(tests$is_cm2 %||% logical()),
      n_three_way_margins = sum(!(tests$is_cm2 %||% logical())),
      source_status = "cm3_values"
    )
    out <- cm2_cm3_apply_parametric_bootstrap(
      fit,
      out,
      context,
      state,
      include_three_way = TRUE,
      control = bootstrap_control
    )
    class(out) <- c("gRm_cm3_values", "list")
    out
  } else {
    out <- list(
      tests = tests,
      aggregate = aggregates$cm2,
      item_trait = aggregates$item_trait,
      invariance = aggregates$invariance,
      cm2_bh = cm2_bh,
      margin_specs = specs,
      selected_items = selection$resolved_items,
      selection = selection,
      exogenous_names = as.character(cm2_cm3_context_backgrounds(context)$name),
      score_cuts = score_cuts,
      score_groups = score_groups,
      n_two_way_margins = sum(tests$is_cm2 %||% logical()),
      n_three_way_margins = 0L,
      source_status = "cm2_values"
    )
    out <- cm2_cm3_apply_parametric_bootstrap(
      fit,
      out,
      context,
      state,
      include_three_way = FALSE,
      control = bootstrap_control
    )
    class(out) <- c("gRm_cm2_values", "list")
    out
  }
}

#' Internal observed CM2/CM3 Benjamini-Hochberg thresholds
#'
#' Reproduces the two separate observed-data p-value blocks in canonical
#' DIGRAM 7.04. Source traces:
#' `source/digram_source_20260817/skunits/skbias14.pas::CM3_tests` and
#' `source/digram_source_20260817/skunits/SKmca.pas::BenjaminiHochberg1`.
#' @param rows Two-way or three-way diagnostic rows containing `p_value`.
#' @param diagnostic Diagnostic label, either `"CM2"` or `"CM3"`.
#' @return A data frame with the three source FDR critical levels and p-value
#'   capacity metadata.
#' @keywords internal
#' @noRd
cm2_cm3_bh_thresholds <- function(rows, diagnostic) {
  if (!is.data.frame(rows) || !nrow(rows)) {
    return(data.frame(
      diagnostic = character(),
      fdr = numeric(),
      critical_p = numeric(),
      n_p_values = integer(),
      source_limit_reached = logical(),
      stringsAsFactors = FALSE
    ))
  }

  # CM3_tests appends only while Npvalues < SKTypes.VECTORLENGTH. It resets
  # Npvalues between the two-way and three-way blocks, so each block has its
  # own 85-value capacity and BH calculation.
  source_limit <- 85L
  n_available <- nrow(rows)
  p_values <- as.numeric(rows$p_value[seq_len(min(n_available, source_limit))])
  fdr <- c(0.05, 0.01, 0.001)
  data.frame(
    diagnostic = rep(diagnostic, length(fdr)),
    fdr = fdr,
    critical_p = vapply(fdr, function(alpha) {
      source_bh_critical(p_values, alpha)
    }, numeric(1L)),
    n_p_values = rep.int(as.integer(length(p_values)), length(fdr)),
    source_limit_reached = rep.int(n_available > source_limit, length(fdr)),
    stringsAsFactors = FALSE
  )
}

#' Internal cm2 cm3 fit context state helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param fit Fitted gRm model.
#' @return The internal `cm2_cm3_fit_context_state()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_fit_context_state <- function(fit) {
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
  # skbias14.CM3_analysis consumes the stopping-step global delta left by the
  # fitted model. Preserve its R `report_delta` counterpart when this fallback
  # reconstructs a diagnostic state instead of silently retaining the
  # initializer's zero.
  state$report_delta <- fit$fit$report_delta %||%
    fit$convergence$report_delta %||% state$delta
  list(context = context, state = state)
}

#' Internal cm2 cm3 selected item table helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param selected_indices Internal `selected_indices` value used by this helper.
#' @return The internal `cm2_cm3_selected_item_table()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_selected_item_table <- function(context, selected_indices) {
  items <- cm2_cm3_context_items(context)
  labels <- cm2_cm3_variable_labels(items)
  data.frame(
    item_index = as.integer(selected_indices),
    item_label = labels[selected_indices],
    item_name = as.character(items$name[selected_indices]),
    stringsAsFactors = FALSE
  )
}

#' Internal cm2 cm3 analyze margin helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @param spec GLLRM model specification.
#' @param score_group_lookup Internal `score_group_lookup` value used by this helper.
#' @param probability_cache Internal `probability_cache` value used by this helper.
#' @return The internal `cm2_cm3_analyze_margin()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_analyze_margin <- function(context, state, spec, score_group_lookup, probability_cache) {
  # Source trace: source/digram_source_20260817/skunits/skbias14.pas::Twoway_analysis and
  # ::Threeway_analysis dispatch observed and expected tables, then compute the
  # Pearson contribution for each prepared CM2/CM3 margin.
  observed <- cm2_cm3_count_observed(context, spec, score_group_lookup)
  expected <- cm2_cm3_expected_table(context, state, spec, score_group_lookup, probability_cache)
  dimensions <- dim(observed)
  df <- if (length(dimensions) == 2L) {
    cm2_cm3_df_two_way(dimensions)
  } else {
    cm2_cm3_df_three_way(dimensions)
  }
  stat <- cm2_cm3_pearson_stat(observed, expected, df = df)
  cbind(
    cm2_cm3_margin_metadata(context, spec),
    data.frame(
      chi_square = stat$chi_square,
      degrees_of_freedom = stat$degrees_of_freedom,
      p_value = stat$p_value,
      stringsAsFactors = FALSE
    )
  )
}

#' Internal cm2 cm3 pearson stat helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param observed Internal `observed` value used by this helper.
#' @param expected Internal `expected` value used by this helper.
#' @param df Internal `df` value used by this helper.
#' @return The internal `cm2_cm3_pearson_stat()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_pearson_stat <- function(observed, expected, df) {
  # Source trace: source/digram_source_20260817/skunits/skbias14.pas::Twoway_analysis and
  # ::Threeway_analysis add (observed - expected)^2 / expected only for cells
  # with expected count greater than zero.
  keep <- expected > 0
  chi_square <- if (any(keep)) {
    sum((observed[keep] - expected[keep])^2 / expected[keep])
  } else {
    0
  }
  df <- as.integer(df)
  # Source trace: source/digram_source_20260817/skunits/skbias14.pas calls PFCHI for the printed
  # CM2/CM3 row and aggregate p-values.
  p_value <- if (df > 0L) source_pfchi(df, chi_square) else 1
  list(
    chi_square = as.numeric(chi_square),
    degrees_of_freedom = df,
    p_value = as.numeric(p_value)
  )
}

#' Internal cm2 cm3 df two way helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param dimensions Internal `dimensions` value used by this helper.
#' @return The internal `cm2_cm3_df_two_way()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_df_two_way <- function(dimensions) {
  # Source trace: source/digram_source_20260817/skunits/skbias14.pas::Twoway_analysis uses the
  # ordinary two-way margin df formula (columns - 1) * (rows - 1).
  dimensions <- as.integer(dimensions)
  as.integer((dimensions[[1L]] - 1L) * (dimensions[[2L]] - 1L))
}

#' Internal cm2 cm3 df three way helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param dimensions Internal `dimensions` value used by this helper.
#' @return The internal `cm2_cm3_df_three_way()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_df_three_way <- function(dimensions) {
  # Source trace: source/digram_source_20260817/skunits/skbias14.pas::Threeway_analysis subtracts
  # the three one-way margins and the grand total from the full table cells.
  dimensions <- as.integer(dimensions)
  as.integer(prod(dimensions) - 1L - sum(dimensions - 1L))
}

#' Internal cm2 cm3 aggregate values helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param tests Diagnostic test rows.
#' @param context Prepared GLLRM computation context.
#' @return The internal `cm2_cm3_aggregate_values()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_aggregate_values <- function(tests, context = NULL) {
  # Source trace: source/digram_source_20260817/skunits/skbias14.pas::CM3_tests stores the CM2 aggregate
  # after the first nM2tests rows, then initializes the CM3 aggregate with that
  # complete CM2 total before adding three-way rows.
  tests <- tests %||% data.frame()
  cm2_rows <- if (is.data.frame(tests) && nrow(tests)) tests[tests$is_cm2, , drop = FALSE] else data.frame()
  cm3_rows <- tests
  item_trait_rows <- if (nrow(cm2_rows)) {
    cm2_rows[cm2_rows$margin_type == "item_score_group", , drop = FALSE]
  } else {
    data.frame()
  }
  item_exogenous_rows <- if (nrow(cm2_rows)) {
    cm2_rows[cm2_rows$margin_type == "item_exogenous", , drop = FALSE]
  } else {
    data.frame()
  }
  list(
    cm2 = cm2_cm3_aggregate_row("CM2", cm2_rows),
    cm3 = cm2_cm3_aggregate_row("CM3", cm3_rows),
    item_trait = cm2_cm3_aggregate_row("Item-trait", item_trait_rows),
    invariance = cm2_cm3_invariance_rows(item_exogenous_rows, context)
  )
}

#' Internal cm2 cm3 aggregate row helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param diagnostic Internal `diagnostic` value used by this helper.
#' @param rows Rows used by the computation.
#' @return The internal `cm2_cm3_aggregate_row()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_aggregate_row <- function(diagnostic, rows) {
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

#' Internal cm2 cm3 invariance rows helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param rows Rows used by the computation.
#' @param context Optional prepared GLLRM computation context. When supplied,
#'   every fitted exogenous variable is retained even if no eligible CM2
#'   item-exogenous margin remains after DIF suppression.
#' @return The internal `cm2_cm3_invariance_rows()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_invariance_rows <- function(rows, context = NULL) {
  # Source trace: skbias14.CM3_tests initializes `invariance[i]`,
  # `invariancedf[i]`, and `NlargerInvariance[i]` to zero for every fitted
  # exogenous variable before it scans the prepared two-way rows. Consequently
  # a variable whose eligible item-exogenous family is empty remains a real
  # result row with chi-square 0, df 0, and asymptotic p-value 0. The final
  # report's `if invariancedf[i]>0 then ... else p:=0` guard is intentionally
  # different from the generic empty-aggregate convention used for CM2/CM3.
  if (!is.null(context)) {
    backgrounds <- cm2_cm3_context_backgrounds(context)
    labels <- cm2_cm3_variable_labels(backgrounds)
    keys <- data.frame(
      background_label = labels,
      background_name = as.character(backgrounds$name),
      stringsAsFactors = FALSE
    )
  } else if (is.data.frame(rows) && nrow(rows)) {
    keep <- !duplicated(rows$background_label)
    keys <- rows[keep, c("background_label", "background_name"), drop = FALSE]
  } else {
    keys <- data.frame(
      background_label = character(),
      background_name = character(),
      stringsAsFactors = FALSE
    )
  }

  if (!nrow(keys)) {
    return(data.frame(
      background_label = character(),
      background_name = character(),
      chi_square = numeric(),
      degrees_of_freedom = integer(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  out <- lapply(seq_len(nrow(keys)), function(index) {
    label <- keys$background_label[[index]]
    name <- keys$background_name[[index]]
    group <- if (is.data.frame(rows) && nrow(rows)) {
      rows[!is.na(rows$background_label) & rows$background_label == label, , drop = FALSE]
    } else {
      data.frame()
    }
    chi_square <- if (nrow(group)) sum(group$chi_square) else 0
    degrees_of_freedom <- if (nrow(group)) sum(group$degrees_of_freedom) else 0L
    p_value <- if (degrees_of_freedom > 0L) {
      source_pfchi(degrees_of_freedom, chi_square)
    } else {
      0
    }
    data.frame(
      background_label = label,
      background_name = name,
      chi_square = as.numeric(chi_square),
      degrees_of_freedom = as.integer(degrees_of_freedom),
      p_value = as.numeric(p_value),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, out)
  rownames(result) <- NULL
  result
}

#' Internal cm2 cm3 margin metadata helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @return The internal `cm2_cm3_margin_metadata()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_margin_metadata <- function(context, spec) {
  items <- cm2_cm3_context_items(context)
  backgrounds <- cm2_cm3_context_backgrounds(context)
  item_labels <- cm2_cm3_variable_labels(items)
  background_labels <- cm2_cm3_variable_labels(backgrounds)
  margin_variables <- cm2_cm3_margin_public_variables(context, spec)
  empty <- NA_character_
  row <- data.frame(
    margin = cm2_cm3_margin_name(context, spec),
    margin_label = cm2_cm3_margin_label(context, spec),
    margin_type = spec$kind,
    is_cm2 = spec$kind %in% c("item_item", "item_exogenous", "item_score_group"),
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
    stop("Unknown CM2/CM3 margin kind.", call. = FALSE)
  )
  row$variable1 <- cm2_cm3_nth_character(margin_variables, 1L, empty)
  row$variable2 <- cm2_cm3_nth_character(margin_variables, 2L, empty)
  row$variable3 <- cm2_cm3_nth_character(margin_variables, 3L, empty)
  row
}

#' Internal cm2 cm3 public margin table helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param tests Diagnostic test rows.
#' @param context Prepared GLLRM computation context.
#' @return The internal `cm2_cm3_public_margin_table()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_public_margin_table <- function(tests, context = NULL) {
  tests
}

#' Internal cm2 cm3 nth character helper
#'
#' Supports the cm2 cm3 values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param x Object or value to process.
#' @param index One-based internal index.
#' @param default Internal `default` value used by this helper.
#' @return The internal `cm2_cm3_nth_character()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_nth_character <- function(x, index, default = NA_character_) {
  if (length(x) < index || is.na(x[[index]])) {
    return(default)
  }
  as.character(x[[index]])
}
